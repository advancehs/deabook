"""Object-oriented core utilities for DEA models.

This module contains the shared, model-agnostic pieces used by the public
``deabook`` model classes: formula parsing, sample/reference data preparation,
direction-vector validation, solver execution and result checks.  The public
model modules keep their historical import paths and class names, while using
these small objects internally to avoid copy-pasted data and solver code.
"""

from __future__ import annotations

import ast
import math
import os
import re
from dataclasses import dataclass
from typing import Any, Iterable

import numpy as np
import pandas as pd
from pyomo.environ import SolverFactory, value
from pyomo.opt import SolverManagerFactory, SolverStatus, TerminationCondition

from .constant import OPT_DEFAULT, OPT_LOCAL, RTS_CRS, RTS_VRS1, RTS_VRS2
from .utils import tools


_FORMULA_SPLIT_RE = re.compile(r"\s+")
_QUERY_OPERATOR_RE = re.compile(r"(==|!=|<=|>=|<|>|\bin\b|\bnot\s+in\b)")
_EMAIL_RE = re.compile(r"([^@]+@[^@]+\.[a-zA-Z0-9]+)$")


@dataclass(frozen=True)
class FormulaSpec:
    """Parsed DEA variable declaration.

    Parameters
    ----------
    x : list[str]
        Input variable names.
    y : list[str]
        Desirable-output variable names.
    b : list[str]
        Undesirable-output variable names. Empty when the model has no
        undesirable outputs.
    raw : str
        Original formula string after separator normalization.
    input_side : str
        ``"left"`` for package-style ``"inputs=outputs:bad"`` declarations,
        or ``"right"`` for teaching-style ``"outputs:bad=inputs"`` formulas.
    """

    x: list[str]
    y: list[str]
    b: list[str]
    raw: str
    input_side: str = "left"

    @property
    def columns(self) -> list[str]:
        """Return all variables in input, desirable-output, bad-output order."""
        return [*self.x, *self.y, *self.b]

    @property
    def has_bad_outputs(self) -> bool:
        """Whether the declaration contains undesirable outputs."""
        return bool(self.b)


@dataclass(frozen=True)
class PreparedData:
    """Evaluation and reference data matrices with original index labels."""

    formula: FormulaSpec
    evaluated_index: list[Any]
    reference_index: list[Any]
    x: np.ndarray
    y: np.ndarray
    xref: np.ndarray
    yref: np.ndarray
    b: np.ndarray | None = None
    bref: np.ndarray | None = None


@dataclass(frozen=True)
class DirectionSpec:
    """Validated direction vectors and their implied orientation."""

    gx: list[float]
    gy: list[float]
    gb: list[float]

    @property
    def has_input_direction(self) -> bool:
        return sum(self.gx) > 0

    @property
    def has_output_direction(self) -> bool:
        return sum(self.gy) > 0

    @property
    def has_bad_direction(self) -> bool:
        return sum(self.gb) > 0

    @property
    def orientation(self) -> str:
        """Return ``input``, ``output``, ``bad`` or one of the hyper modes."""
        active = [
            name
            for name, enabled in (
                ("input", self.has_input_direction),
                ("output", self.has_output_direction),
                ("bad", self.has_bad_direction),
            )
            if enabled
        ]
        if not active:
            raise ValueError("At least one direction vector must contain a positive element.")
        if len(active) == 1:
            return active[0]
        return "hyper_" + "_".join(active)


@dataclass(frozen=True)
class SolverOutcome:
    """Normalized Pyomo solver outcome."""

    status: str
    termination_condition: str
    ok: bool
    raw: Any

    @property
    def legacy_status(self) -> str:
        """Status label historically used by deabook result tables."""
        return "ok" if self.ok else self.termination_condition or self.status


@dataclass(frozen=True)
class SolverConfig:
    """Solver settings shared by model classes.

    Parameters
    ----------
    solver : str or None
        Pyomo solver name. ``None`` keeps deabook's historical default of
        ``"mosek"`` for additive linear models.
    email : str
        ``OPT_LOCAL`` for local solving, otherwise a NEOS email address.
    tee : bool
        Whether Pyomo should print solver logs.
    """

    solver: str | None = OPT_DEFAULT
    email: str = OPT_LOCAL
    tee: bool = False

    @property
    def solver_name(self) -> str:
        """Concrete solver name used for local solving."""
        return "mosek" if self.solver is OPT_DEFAULT else str(self.solver)

    @property
    def use_neos(self) -> bool:
        """Whether to solve through NEOS."""
        return self.email != OPT_LOCAL


@dataclass(frozen=True)
class WeakDirectionFlags:
    """Boolean orientation flags for models with undesirable outputs.

    The historical weak-disposability modules use names such as
    ``hyper_orientedyb``.  Centralising those names keeps the public behaviour
    while avoiding repeated ad-hoc ``sum(g*)`` checks in each model class.
    """

    input_oriented: bool
    output_oriented: bool
    unoutput_oriented: bool
    hyper_orientedyx: bool
    hyper_orientedyb: bool
    hyper_orientedxb: bool
    hyper_orientedyxb: bool

    @property
    def is_single(self) -> bool:
        return self.input_oriented or self.output_oriented or self.unoutput_oriented

    @property
    def is_hyper(self) -> bool:
        return (
            self.hyper_orientedyx
            or self.hyper_orientedyb
            or self.hyper_orientedxb
            or self.hyper_orientedyxb
        )

    @property
    def active_efficiency_columns(self) -> list[str]:
        columns: list[str] = []
        if self.input_oriented or self.hyper_orientedyx or self.hyper_orientedxb or self.hyper_orientedyxb:
            columns.append("tei")
        if self.unoutput_oriented or self.hyper_orientedyb or self.hyper_orientedxb or self.hyper_orientedyxb:
            columns.append("teuo")
        if self.output_oriented or self.hyper_orientedyx or self.hyper_orientedyb or self.hyper_orientedyxb:
            columns.append("teo")
        return columns


def build_weak_direction_flags(gx: Iterable[float], gy: Iterable[float], gb: Iterable[float]) -> WeakDirectionFlags:
    """Return legacy-compatible orientation flags for weak DEA models."""
    sum_gx = sum(gx)
    sum_gy = sum(gy)
    sum_gb = sum(gb)
    flags = WeakDirectionFlags(
        input_oriented=sum_gx >= 1 and sum_gy == 0 and sum_gb == 0,
        output_oriented=sum_gy >= 1 and sum_gx == 0 and sum_gb == 0,
        unoutput_oriented=sum_gb >= 1 and sum_gx == 0 and sum_gy == 0,
        hyper_orientedyx=sum_gx >= 1 and sum_gy >= 1 and sum_gb == 0,
        hyper_orientedyb=sum_gb >= 1 and sum_gy >= 1 and sum_gx == 0,
        hyper_orientedxb=sum_gb >= 1 and sum_gx >= 1 and sum_gy == 0,
        hyper_orientedyxb=sum_gb >= 1 and sum_gx >= 1 and sum_gy >= 1,
    )
    if not (flags.is_single or flags.is_hyper):
        raise ValueError(
            "gx, gy and gb must activate at least one valid input, output, "
            "undesirable-output, or hyper orientation."
        )
    return flags


def normalize_formula(formula: str) -> str:
    """Normalize a DEA formula without changing variable-name case.

    Parameters
    ----------
    formula : str
        Raw formula string.

    Returns
    -------
    str
        Formula with full-width separators converted and repeated whitespace
        collapsed.
    """
    if not isinstance(formula, str):
        raise TypeError("formula/sent must be a string.")
    normalized = formula.replace("：", ":").replace("＝", "=").replace("＋", "+")
    normalized = _FORMULA_SPLIT_RE.sub(" ", normalized).strip()
    if not normalized:
        raise ValueError("formula/sent cannot be empty.")
    return normalized


def split_variables(text: str) -> list[str]:
    """Split a whitespace-delimited variable list."""
    return [part for part in _FORMULA_SPLIT_RE.split(text.strip()) if part]


def parse_sent(sent: str) -> FormulaSpec:
    """Parse package-style ``"inputs=outputs[:bad_outputs]"`` syntax."""
    return parse_dea_formula(sent, input_side="left")


def parse_formula(formula: str) -> FormulaSpec:
    """Parse teaching-style ``"outputs[:bad_outputs]=inputs"`` syntax."""
    return parse_dea_formula(formula, input_side="right")


def parse_dea_formula(formula: str, input_side: str = "left") -> FormulaSpec:
    """Parse DEA variable declarations.

    Parameters
    ----------
    formula : str
        DEA variable declaration.
    input_side : {"left", "right"}
        Whether inputs appear on the left or right side of ``=``.

    Returns
    -------
    FormulaSpec
        Parsed variable groups.
    """
    normalized = normalize_formula(formula)
    parts = normalized.split("=", 1)
    if len(parts) != 2 or not parts[0].strip() or not parts[1].strip():
        raise ValueError(
            f"Invalid formula '{formula}'. Expected exactly one '=' with variables on both sides."
        )

    left, right = parts[0].strip(), parts[1].strip()
    if input_side not in {"left", "right"}:
        raise ValueError("input_side must be either 'left' or 'right'.")

    input_part = left if input_side == "left" else right
    output_part = right if input_side == "left" else left
    output_parts = output_part.split(":", 1)

    x = split_variables(input_part)
    y = split_variables(output_parts[0])
    b = split_variables(output_parts[1]) if len(output_parts) == 2 else []

    if not x:
        raise ValueError("At least one input variable is required.")
    if not y:
        raise ValueError("At least one desirable output variable is required.")
    duplicates = _duplicates([*x, *y, *b])
    if duplicates:
        raise ValueError(f"Variables cannot be repeated in a formula: {duplicates}")
    return FormulaSpec(x=x, y=y, b=b, raw=normalized, input_side=input_side)


def prepare_dea_data(
    dataframe: pd.DataFrame,
    formula: FormulaSpec | str,
    eval_query: str | None = None,
    ref_query: str | None = None,
    input_side: str = "left",
    require_positive_x: bool = True,
    require_positive_y: bool = True,
    require_nonnegative_b: bool = True,
) -> PreparedData:
    """Prepare evaluation and reference matrices for DEA models.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        Source data. Rows are DMUs or DMU-period observations.
    formula : FormulaSpec or str
        Parsed formula, or a raw formula string.
    eval_query, ref_query : str or None
        Sample and reference selectors. Supports pandas ``query`` syntax and
        legacy deabook selectors such as ``"Year=[2010,2011]"``.
    input_side : {"left", "right"}
        Used only when ``formula`` is a string.

    Returns
    -------
    PreparedData
        Numeric matrices and original index labels.
    """
    if not isinstance(dataframe, pd.DataFrame):
        raise TypeError("dataframe/data must be a pandas.DataFrame.")
    spec = parse_dea_formula(formula, input_side=input_side) if isinstance(formula, str) else formula
    if not isinstance(spec, FormulaSpec):
        raise TypeError("formula must be a string or FormulaSpec.")

    missing = [column for column in spec.columns if column not in dataframe.columns]
    if missing:
        raise KeyError(
            "Missing DEA variable columns: "
            f"{missing}. Available columns: {list(dataframe.columns)}"
        )

    eval_df = select_rows(dataframe, eval_query, "eval_query/baseindex")
    ref_df = select_rows(dataframe, ref_query, "ref_query/refindex")
    if eval_df.empty:
        raise ValueError("The evaluation sample is empty.")
    if ref_df.empty:
        raise ValueError("The reference sample is empty.")

    x = numeric_matrix(eval_df, spec.x, "evaluated inputs", strictly_positive=require_positive_x)
    y = numeric_matrix(eval_df, spec.y, "evaluated desirable outputs", strictly_positive=require_positive_y)
    xref = numeric_matrix(ref_df, spec.x, "reference inputs", strictly_positive=require_positive_x)
    yref = numeric_matrix(ref_df, spec.y, "reference desirable outputs", strictly_positive=require_positive_y)

    b = bref = None
    if spec.b:
        b = numeric_matrix(
            eval_df,
            spec.b,
            "evaluated undesirable outputs",
            nonnegative=require_nonnegative_b,
        )
        bref = numeric_matrix(
            ref_df,
            spec.b,
            "reference undesirable outputs",
            nonnegative=require_nonnegative_b,
        )

    return PreparedData(
        formula=spec,
        evaluated_index=eval_df.index.tolist(),
        reference_index=ref_df.index.tolist(),
        x=x,
        y=y,
        xref=xref,
        yref=yref,
        b=b,
        bref=bref,
    )


def select_rows(dataframe: pd.DataFrame, selector: str | None, name: str) -> pd.DataFrame:
    """Select rows using either pandas query or legacy ``column=[values]`` syntax."""
    if selector is None:
        return dataframe.copy()
    if not isinstance(selector, str):
        raise TypeError(f"{name} must be a string or None.")
    selector = selector.strip()
    if not selector:
        raise ValueError(f"{name} cannot be an empty string.")

    legacy_match = re.fullmatch(r"\s*([A-Za-z_][\w]*)\s*=\s*(.+?)\s*", selector)
    if legacy_match and not _QUERY_OPERATOR_RE.search(selector):
        column, raw_value = legacy_match.groups()
        if column not in dataframe.columns:
            raise KeyError(
                f"Selector column '{column}' in {name} is not in the data. "
                f"Available columns: {list(dataframe.columns)}"
            )
        try:
            values = ast.literal_eval(raw_value)
        except (SyntaxError, ValueError) as exc:
            raise ValueError(f"Invalid value in {name}: {raw_value!r}") from exc
        if not isinstance(values, (list, tuple, set, pd.Index, np.ndarray)):
            values = [values]
        return dataframe.loc[dataframe[column].isin(values)].copy()

    try:
        return dataframe.query(selector).copy()
    except Exception as exc:
        raise ValueError(f"Invalid pandas query in {name}: {selector!r}") from exc


def numeric_matrix(
    dataframe: pd.DataFrame,
    columns: list[str],
    label: str,
    strictly_positive: bool = False,
    nonnegative: bool = False,
) -> np.ndarray:
    """Return a finite numeric matrix and enforce DEA sign restrictions."""
    try:
        matrix = dataframe.loc[:, columns].to_numpy(dtype=float)
    except Exception as exc:
        raise ValueError(f"{label} must be numeric. Columns: {columns}") from exc
    if matrix.ndim != 2:
        raise ValueError(f"{label} must be two-dimensional.")
    if matrix.shape[0] == 0 or matrix.shape[1] == 0:
        raise ValueError(f"{label} cannot be empty.")
    if not np.isfinite(matrix).all():
        raise ValueError(f"{label} contains NaN or infinite values.")
    if strictly_positive and not (matrix > 0).all():
        raise ValueError(f"{label} must be strictly positive.")
    if nonnegative and not (matrix >= 0).all():
        raise ValueError(f"{label} cannot contain negative values.")
    return matrix


def validate_direction_vector(vector: Iterable[float] | float, expected_length: int, name: str) -> list[float]:
    """Validate and normalize a DEA direction vector."""
    if expected_length < 0:
        raise ValueError("expected_length cannot be negative.")
    if expected_length == 0:
        values = [] if vector in (None, [], ()) else list(vector) if isinstance(vector, (list, tuple)) else [vector]
        if any(float(item) != 0 for item in values):
            raise ValueError(f"{name} must be empty when the corresponding variable group is empty.")
        return []

    if isinstance(vector, (int, float, np.integer, np.floating)):
        values = [float(vector)]
    else:
        try:
            values = [float(item) for item in vector]
        except TypeError as exc:
            raise TypeError(f"{name} must be a numeric scalar or iterable.") from exc

    if len(values) == 1 and expected_length > 1 and values[0] == 0:
        values = [0.0] * expected_length
    if len(values) != expected_length:
        raise ValueError(
            f"Length of {name} ({len(values)}) must match the number of "
            f"corresponding variables ({expected_length})."
        )
    if any(not math.isfinite(item) for item in values):
        raise ValueError(f"{name} contains NaN or infinite values.")
    if any(item < 0 for item in values):
        raise ValueError(f"{name} cannot contain negative values.")
    return values


def build_direction_spec(
    gx: Iterable[float] | float,
    gy: Iterable[float] | float,
    gb: Iterable[float] | float | None,
    num_inputs: int,
    num_outputs: int,
    num_bad_outputs: int = 0,
    allow_bad_direction: bool = False,
) -> DirectionSpec:
    """Validate ``gx``, ``gy`` and ``gb`` against parsed formula dimensions."""
    gx_values = validate_direction_vector(gx, num_inputs, "gx")
    gy_values = validate_direction_vector(gy, num_outputs, "gy")
    gb_values = validate_direction_vector(gb or [], num_bad_outputs, "gb")
    if gb_values and not allow_bad_direction:
        raise ValueError("This model does not support undesirable-output directions.")
    spec = DirectionSpec(gx=gx_values, gy=gy_values, gb=gb_values)
    _ = spec.orientation
    return spec


def validate_rts(rts: str, allowed: set[str] | None = None) -> str:
    """Validate returns-to-scale code."""
    allowed_values = allowed or {RTS_CRS, RTS_VRS1, RTS_VRS2}
    if rts not in allowed_values:
        raise ValueError(f"rts must be one of {sorted(allowed_values)}; got {rts!r}.")
    return rts


def solve_pyomo_model(model: Any, config: SolverConfig | None = None, strict: bool = True) -> SolverOutcome:
    """Solve a Pyomo model with availability and optimality checks.

    Parameters
    ----------
    model : pyomo.environ.ConcreteModel
        Built Pyomo model.
    config : SolverConfig or None
        Solver settings. ``None`` uses ``SolverConfig()``.
    strict : bool
        If ``True``, raise ``RuntimeError`` for unavailable solvers and
        non-optimal results.

    Returns
    -------
    SolverOutcome
        Normalized status information.
    """
    config = config or SolverConfig()
    solver_name = config.solver_name

    if config.use_neos:
        if not _EMAIL_RE.match(config.email):
            raise ValueError("Invalid NEOS email address.")
        os.environ["NEOS_EMAIL"] = config.email
        manager = SolverManagerFactory("neos")
        results = manager.solve(model, tee=config.tee, opt=solver_name)
    else:
        solver = SolverFactory(solver_name)
        available = False
        try:
            available = bool(solver.available(exception_flag=False))
        except Exception:
            available = False
        if not available:
            message = (
                f"Solver '{solver_name}' is not installed or not available. "
                "Install it or pass a different solver, for example solver='glpk' or solver='scip'."
            )
            if strict:
                raise RuntimeError(message)
            return SolverOutcome(status="unavailable", termination_condition="unavailable", ok=False, raw=None)
        results = solver.solve(model, tee=config.tee)

    status_obj = results.solver.status
    termination_obj = results.solver.termination_condition
    status = getattr(status_obj, "name", str(status_obj)).lower()
    termination = getattr(termination_obj, "name", str(termination_obj)).lower()
    ok = status_obj == SolverStatus.ok and termination_obj == TerminationCondition.optimal
    if strict and not ok:
        if termination_obj == TerminationCondition.infeasible:
            detail = "The linear program is infeasible; check variables, signs, and reference set."
        elif termination_obj == TerminationCondition.unbounded:
            detail = "The linear program is unbounded; check model orientation and normalization constraints."
        else:
            detail = f"Solver status={status}, termination_condition={termination}."
        raise RuntimeError(f"Solver did not reach an optimal solution. {detail}")
    return SolverOutcome(status=status, termination_condition=termination, ok=ok, raw=results)


def pyomo_value(component: Any) -> float:
    """Return a finite float value from a Pyomo component."""
    val = value(component)
    if val is None or not math.isfinite(float(val)):
        raise RuntimeError("The optimization result contains a missing or non-finite value.")
    return float(val)


def clean_small_values(values: Iterable[float], atol: float = 1e-9) -> np.ndarray:
    """Convert near-zero numerical noise to exact zero."""
    array = np.asarray(list(values), dtype=float)
    array[np.isclose(array, 0.0, atol=atol)] = 0.0
    return array


def indexed_component_to_frame(component: Any, value_accessor: str = "value") -> pd.DataFrame:
    """Convert a one- or two-dimensional Pyomo indexed component to a DataFrame.

    This helper mirrors deabook's repeated ``list(var[:, :].value)`` plus
    ``pivot`` pattern while handling missing values explicitly.
    """
    rows: list[dict[str, Any]] = []
    for key in component:
        raw_key = key if isinstance(key, tuple) else (key,)
        obj = component[key]
        val = getattr(obj, value_accessor, None)
        if val is None:
            val = value(obj)
        rows.append({"row": raw_key[0], "column": raw_key[1] if len(raw_key) > 1 else 0, "value": val})
    if not rows:
        return pd.DataFrame()
    frame = pd.DataFrame(rows)
    return frame.pivot(index="row", columns="column", values="value")


def indexed_component_to_numpy(component: Any) -> np.ndarray:
    """Return an indexed Pyomo component as a NumPy array."""
    return indexed_component_to_frame(component).to_numpy()


def scalar_component_to_series(component: Any) -> pd.Series:
    """Return a one-dimensional Pyomo component as a Series."""
    return indexed_component_to_frame(component).iloc[:, 0]


def format_ref_indexed_results(
    values: dict[Any, dict[Any, float] | None],
    reference_index: Iterable[Any],
) -> pd.DataFrame:
    """Format per-DMU dictionaries keyed by reference DMU into a DataFrame."""
    reference_index = list(reference_index)
    rows = []
    for evaluated_index, data_for_dmu in values.items():
        if data_for_dmu is None:
            rows.append(pd.Series(np.nan, index=reference_index, name=evaluated_index))
        else:
            rows.append(pd.Series(data_for_dmu, name=evaluated_index))
    if not rows:
        return pd.DataFrame(columns=reference_index)
    return pd.concat(rows, axis=1).T


def format_var_indexed_results(
    values: dict[Any, dict[Any, float] | list[float] | np.ndarray | None],
    columns: Iterable[Any] | None = None,
) -> pd.DataFrame:
    """Format per-DMU dictionaries or vectors keyed by variable index."""
    normalized: dict[Any, Any] = {}
    for evaluated_index, data_for_dmu in values.items():
        if data_for_dmu is None:
            normalized[evaluated_index] = np.nan
        else:
            normalized[evaluated_index] = data_for_dmu
    frame = pd.DataFrame(normalized).T
    if columns is not None and len(frame.columns) == len(list(columns)):
        frame.columns = list(columns)
    return frame


class ModelDictionaryMixin:
    """Mixin for classes that hold one Pyomo model per DMU."""

    _models_attr = "_models"

    def _model_dict(self) -> dict[Any, Any]:
        if hasattr(self, self._models_attr):
            return getattr(self, self._models_attr)
        for name in dir(self):
            if name.endswith("__modeldict"):
                return getattr(self, name)
        raise AttributeError("This object does not expose a model dictionary.")

    def iter_models(self):
        """Yield ``(dmu, model)`` pairs."""
        yield from self._model_dict().items()

    def info(self, dmu: Any = "all"):
        """Print Pyomo model details for selected DMUs."""
        models = self._model_dict()
        if dmu == "all":
            for ind, problem in models.items():
                print(f"\n--- Model for DMU: {ind} ---")
                problem.pprint()
            return None
        dmu_list = [dmu] if isinstance(dmu, (str, int, float)) else list(dmu)
        for ind in dmu_list:
            key = ind
            if key not in models:
                try:
                    key = int(ind)
                except (TypeError, ValueError):
                    pass
            if key in models:
                print(f"\n--- Model for DMU: {key} ---")
                models[key].pprint()
            else:
                print(f"DMU '{ind}' not found in the evaluated set.")
        return None


class PyomoResultMixin:
    """Mixin with reusable Pyomo variable display and extraction methods."""

    _model_attr = "_model"

    def _single_model(self):
        if hasattr(self, self._model_attr):
            return getattr(self, self._model_attr)
        for cls in self.__class__.mro():
            mangled = f"_{cls.__name__}__model__"
            if hasattr(self, mangled):
                return getattr(self, mangled)
        for name in dir(self):
            if name.endswith("__model__"):
                return getattr(self, name)
        raise AttributeError("This object does not expose a Pyomo model.")

    def _display_component(self, name: str):
        component = getattr(self._single_model(), name)
        component.display()

    def _component_array(self, name: str) -> np.ndarray:
        return indexed_component_to_numpy(getattr(self._single_model(), name))

    def _component_series(self, name: str) -> np.ndarray:
        component = getattr(self._single_model(), name)
        return np.asarray(list(component[:].value))


class CNLSResultMixin(PyomoResultMixin):
    """Common public result accessors for CNLS-style models."""

    def _optimize_cnls_model(self, email=OPT_LOCAL, solver=OPT_DEFAULT, cet=None):
        """Solve the underlying CNLS Pyomo model using legacy status fields."""
        self.problem_status, self.optimization_status = tools.optimize_model(
            self._single_model(), email, cet if cet is not None else getattr(self, "cet", None), solver
        )

    def display_status(self):
        tools.assert_optimized(self.optimization_status)
        print(self.optimization_status)

    def display_alpha(self):
        tools.assert_optimized(self.optimization_status)
        tools.assert_various_return_to_scale(self.rts)
        self._display_component("alpha")

    def display_delta(self):
        tools.assert_optimized(self.optimization_status)
        self._display_component("delta")

    def display_gamma(self):
        tools.assert_optimized(self.optimization_status)
        self._display_component("gamma")

    def display_kappa(self):
        tools.assert_optimized(self.optimization_status)
        self._display_component("kappa")

    def display_omega(self):
        tools.assert_optimized(self.optimization_status)
        tools.assert_contextual_variable(self.z)
        self._display_component("omega")

    def display_residual(self):
        tools.assert_optimized(self.optimization_status)
        self._display_component("epsilon")

    def get_status(self):
        return self.optimization_status

    def get_alpha(self):
        tools.assert_optimized(self.optimization_status)
        tools.assert_various_return_to_scale(self.rts)
        return self._component_series("alpha")

    def get_delta(self):
        tools.assert_optimized(self.optimization_status)
        return self._component_array("delta")

    def get_gamma(self):
        tools.assert_optimized(self.optimization_status)
        return self._component_array("gamma")

    def get_kappa(self):
        tools.assert_optimized(self.optimization_status)
        return self._component_array("kappa")

    def get_omega(self):
        tools.assert_optimized(self.optimization_status)
        tools.assert_contextual_variable(self.z)
        return self._component_series("omega")

    def get_adjusted_residual(self):
        tools.assert_optimized(self.optimization_status)
        return self.get_residual() - np.amax(self.get_residual())

    def get_adjusted_alpha(self):
        tools.assert_optimized(self.optimization_status)
        return self.get_alpha() + np.amax(self.get_residual())


def _duplicates(items: list[str]) -> list[str]:
    seen: set[str] = set()
    repeated: list[str] = []
    for item in items:
        if item in seen and item not in repeated:
            repeated.append(item)
        seen.add(item)
    return repeated
