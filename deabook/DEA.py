# import dependencies
from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, maximize, Constraint, value
import numpy as np
import pandas as pd
from .constant import (
    CONTEMPORARY,
    GLOBAL,
    OPT_DEFAULT,
    OPT_LOCAL,
    RTS_CRS,
    RTS_VRS,
    RTS_VRS1,
    SEQUENTIAL,
    WINDOW,
)


# ---------------------------------------------------------------------------
# OOP refactor for the public per-DMU models used throughout the book.
# ---------------------------------------------------------------------------

from .core import (
    SolverConfig,
    build_direction_spec,
    clean_small_values,
    parse_sent,
    prepare_dea_data,
    pyomo_value,
    solve_pyomo_model,
    validate_rts,
)

# Malmquist engine (leaf module; importing it here cannot create a cycle).
from .mqcore import (
    MalmquistModel,
    _shift,
    build_distance_frame,
    fgnz_decompose,
    global_decompose,
    normalize_tech,
    normalize_window,
    pl_global_decompose,
    rd_decompose,
    two_factor,
)

# Modern SBM implementation lives below as the private ``_SBM`` class;
# ``sbm`` / ``sbmu`` delegate to it.


class _PerDmuEnvelopmentModel:
    """Base class for one-model-per-DMU DEA estimators.

    Subclasses implement only the model-specific variable, objective and
    constraint rules.  Data parsing, sample/reference splitting, solver calls,
    result tables and ``info``/getter methods are shared here.
    """

    allowed_rts = {RTS_CRS, RTS_VRS, RTS_VRS1}
    result_value_name = "rho"

    def __init__(self, data, sent, gy, gx, rts=RTS_CRS, evalindex=None, refindex=None):
        self.data = data
        self.sent = sent
        self.rts = validate_rts(rts, self.allowed_rts)
        self.evalindex = evalindex
        self.refindex = refindex
        self.prepared = prepare_dea_data(
            data,
            sent,
            eval_query=evalindex,
            ref_query=refindex,
            input_side="left",
        )
        self.formula = self.prepared.formula
        self.xcol = self.formula.x
        self.ycol = self.formula.y
        self.y = self.prepared.y
        self.x = self.prepared.x
        self.yref = self.prepared.yref
        self.xref = self.prepared.xref
        self.evaluated_data_index = self.prepared.evaluated_index
        self.reference_data_index = self.prepared.reference_index
        self.direction = build_direction_spec(
            gx,
            gy,
            None,
            num_inputs=len(self.xcol),
            num_outputs=len(self.ycol),
        )
        self.gx = self.direction.gx
        self.gy = self.direction.gy
        self.input_oriented = self.direction.orientation == "input"
        self.output_oriented = self.direction.orientation == "output"
        self.hyper_oriented = self.direction.orientation == "hyper_input_output"
        self.evaluated_indices_range = range(len(self.evaluated_data_index))
        self.reference_indices_range = range(len(self.reference_data_index))
        self.num_inputs = len(self.xcol)
        self.num_outputs = len(self.ycol)
        self._models = {}
        self.results = {}
        for position, actual_index in enumerate(self.evaluated_data_index):
            self._models[actual_index] = self._build_model(position, actual_index)

    def _build_model(self, current_dmu_range_index, actual_index):
        model = ConcreteModel()
        model.R = Set(initialize=self.reference_indices_range)
        model.J = Set(initialize=range(self.num_inputs))
        model.K = Set(initialize=range(self.num_outputs))
        self._add_decision_variables(model, actual_index)
        self._add_objective(model)
        model.input = Constraint(
            model.J,
            rule=self._create_input_rule(current_dmu_range_index),
            doc="input constraint",
        )
        model.output = Constraint(
            model.K,
            rule=self._create_output_rule(current_dmu_range_index),
            doc="output constraint",
        )
        if self.rts in (RTS_VRS, RTS_VRS1):
            model.vrs = Constraint(rule=self._create_vrs_rule(), doc="variable return to scale rule")
        return model

    def _add_decision_variables(self, model, actual_index):
        raise NotImplementedError

    def _add_objective(self, model):
        raise NotImplementedError

    def _create_input_rule(self, current_dmu_range_index):
        raise NotImplementedError

    def _create_output_rule(self, current_dmu_range_index):
        raise NotImplementedError

    def _create_vrs_rule(self):
        def vrs_rule(model):
            return sum(model.lamda[r] for r in model.R) == 1
        return vrs_rule

    def _efficiency_columns(self, result_df):
        return result_df

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Solve every evaluated DMU and return a result table.

        Parameters
        ----------
        email : str
            ``OPT_LOCAL`` for local solving; otherwise a NEOS email address.
        solver : str or None
            Pyomo solver name. ``None`` keeps the historical default solver.

        Returns
        -------
        pandas.DataFrame
            Rows use the original evaluated-data index and include
            ``optimization_status``, ``rho`` and model-specific efficiency
            columns.
        """
        status_values = {}
        objective_values = {}
        scalar_values = {}
        lamda_values = {}
        config = SolverConfig(solver=solver, email=email, tee=False)
        for actual_index, model in self._models.items():
            try:
                outcome = solve_pyomo_model(model, config=config, strict=False)
                status_values[actual_index] = outcome.legacy_status
                if outcome.ok:
                    scalar_values[actual_index] = pyomo_value(getattr(model, self.result_value_name))
                    objective_values[actual_index] = pyomo_value(model.objective)
                    lamda_values[actual_index] = {
                        self.reference_data_index[r]: pyomo_value(model.lamda[r])
                        for r in model.R
                    }
                else:
                    scalar_values[actual_index] = np.nan
                    objective_values[actual_index] = np.nan
                    lamda_values[actual_index] = None
            except Exception as exc:
                status_values[actual_index] = f"error: {exc}"
                scalar_values[actual_index] = np.nan
                objective_values[actual_index] = np.nan
                lamda_values[actual_index] = None

        self.results = {
            "optimization_status": status_values,
            self.result_value_name: scalar_values,
            "objective_value": objective_values,
            "lamda_df": self._format_lamda(lamda_values),
        }
        result_df = pd.DataFrame(
            {
                "optimization_status": pd.Series(status_values),
                self.result_value_name: pd.Series(scalar_values, dtype="float64"),
                "objective_value": pd.Series(objective_values, dtype="float64"),
            }
        )
        result_df = self._efficiency_columns(result_df)
        return result_df

    def _format_lamda(self, lamda_values):
        rows = []
        for actual_index, values in lamda_values.items():
            if values is None:
                rows.append(pd.Series(np.nan, index=self.reference_data_index, name=actual_index))
            else:
                rows.append(pd.Series(values, name=actual_index))
        if not rows:
            return pd.DataFrame(index=self.evaluated_data_index, columns=self.reference_data_index)
        return pd.concat(rows, axis=1).T

    def display_status(self):
        """Display solver status for each DMU."""
        print(pd.Series(self.get_status()))

    def display_rho(self):
        """Display rho values for each DMU."""
        print(self.get_rho())

    def display_lamda(self):
        """Display intensity variables for each DMU."""
        print(self.get_lamda())

    def get_status(self):
        """Return solver status by original evaluated-data index."""
        return self.results.get("optimization_status", {})

    def get_rho(self):
        """Return rho as a pandas Series indexed by original DMU labels."""
        if not self.results:
            raise RuntimeError("The model has not been optimized yet.")
        return pd.Series(self.results.get(self.result_value_name, {}), dtype="float64")

    def get_lamda(self):
        """Return intensity variables as a DataFrame."""
        if not self.results:
            raise RuntimeError("The model has not been optimized yet.")
        return self.results.get("lamda_df")

    def info(self, dmu="all"):
        """Print Pyomo model details for selected DMUs."""
        if dmu == "all":
            dmu_list = list(self._models)
        elif isinstance(dmu, (str, int, float)):
            dmu_list = [dmu]
        else:
            dmu_list = list(dmu)
        for actual_index in dmu_list:
            if actual_index not in self._models:
                print(f"DMU '{actual_index}' not found in the evaluated set.")
                continue
            print(f"\n--- Model for DMU: {actual_index} ---")
            self._models[actual_index].pprint()


class _DdfRho(_PerDmuEnvelopmentModel):
    """Directional distance function solved separately for each evaluated DMU.

    The public :func:`ddf` factory returns this class.  Chapter 3 derives the
    efficiency from ``rho`` at the call site because it depends on the chosen
    direction vector (``te = 1 - rho`` on an input direction, ``te = 1 / (1 +
    rho)`` on an output direction).
    """

    def _add_decision_variables(self, model, actual_index):
        model.rho = Var(bounds=(None, None), doc=f"directional distance for DMU {actual_index}")
        model.lamda = Var(model.R, bounds=(0.0, None), doc="intensity variables")

    def _add_objective(self, model):
        model.objective = Objective(rule=lambda m: m.rho, sense=maximize, doc="objective function")

    def _create_input_rule(self, current_dmu_range_index):
        def input_rule(model, j):
            return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <= (
                self.x[current_dmu_range_index][j]
                - model.rho * self.gx[j] * self.x[current_dmu_range_index][j]
            )
        return input_rule

    def _create_output_rule(self, current_dmu_range_index):
        def output_rule(model, k):
            return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= (
                self.y[current_dmu_range_index][k]
                + model.rho * self.gy[k] * self.y[current_dmu_range_index][k]
            )
        return output_rule

    def _efficiency_columns(self, result_df):
        return result_df



# ---------------------------------------------------------------------------
# Decoupled per-interface models.
# (hyperbolic), so callers do not need to specify ``gy``/``gx``.
# ---------------------------------------------------------------------------

_ORIENT_ALIASES = {
    "oo": "output",
    "output": "output",
    "io": "input",
    "input": "input",
}


def _normalize_orient(orient):
    direction = _ORIENT_ALIASES.get(orient)
    if direction is None:
        raise ValueError("orient must be 'oo'/'output' or 'io'/'input'.")
    return direction


def _radial_direction(sent, orient):
    """Return ``gx``/``gy`` for a radial model from an orientation string."""
    spec = parse_sent(sent)
    n_inputs, n_outputs = len(spec.x), len(spec.y)
    if _normalize_orient(orient) == "input":
        return [1] * n_inputs, [0] * n_outputs
    return [0] * n_inputs, [1] * n_outputs


# ---------------------------------------------------------------------------
# Decoupled per-interface models.  Each subclasses ``_PerDmuEnvelopmentModel``
# directly (like ``DDF``) and implements exactly one linear program, so the
# convenience interfaces (``ccr``/``bcc``/``hyperbolic``/``ite``) no longer
# route through the monolithic ``DEA`` class.  ``DEA`` itself is kept
# unchanged for ``DEAweak`` / ``deabook.DEAModel`` / ``MQDEA``.
# ---------------------------------------------------------------------------


class _RadialModel(_PerDmuEnvelopmentModel):
    """Radial envelopment model shared by CCR, BCC and ITE.

    Implements the classic radial LP for one orientation; the direction
    vectors ``gx``/``gy`` (chosen by the subclass) decide whether it is
    input- or output-oriented.  There is no hyperbolic branch.
    """

    def _add_decision_variables(self, model, actual_index):
        bounds = (0.0, None)  # Cross-reference distances may exceed in-sample bounds.
        model.rho = Var(bounds=bounds, doc=f"efficiency for DMU {actual_index}")
        model.lamda = Var(model.R, bounds=(0.0, None), doc="intensity variables")

    def _add_objective(self, model):
        sense = minimize if self.input_oriented else maximize
        model.objective = Objective(rule=lambda m: m.rho, sense=sense, doc="objective function")

    def _create_input_rule(self, current_dmu_range_index):
        def input_rule(model, j):
            lhs = sum(model.lamda[r] * self.xref[r][j] for r in model.R)
            x0 = self.x[current_dmu_range_index][j]
            if self.input_oriented:
                return lhs <= (model.rho * x0 if self.gx[j] > 0 else x0)
            return lhs <= x0
        return input_rule

    def _create_output_rule(self, current_dmu_range_index):
        def output_rule(model, k):
            lhs = sum(model.lamda[r] * self.yref[r][k] for r in model.R)
            y0 = self.y[current_dmu_range_index][k]
            if self.output_oriented:
                return lhs >= (model.rho * y0 if self.gy[k] > 0 else y0)
            return lhs >= y0
        return output_rule

    def _efficiency_columns(self, result_df):
        rho = result_df["rho"]
        if self.input_oriented:
            result_df["te"] = rho
        else:
            result_df["te"] = rho.apply(lambda item: 1 / item if pd.notna(item) and item != 0 else np.nan)
        return result_df


class CCR(_RadialModel):
    """CCR radial DEA (constant returns to scale)."""

    allowed_rts = {RTS_CRS}

    def __init__(self, data, sent, orient="oo", evalindex=None, refindex=None):
        gx, gy = _radial_direction(sent, orient)
        super().__init__(data, sent, gy=gy, gx=gx, rts=RTS_CRS, evalindex=evalindex, refindex=refindex)


class BCC(_RadialModel):
    """BCC radial DEA (variable returns to scale)."""

    allowed_rts = {RTS_VRS}

    def __init__(self, data, sent, orient="oo", evalindex=None, refindex=None):
        gx, gy = _radial_direction(sent, orient)
        super().__init__(data, sent, gy=gy, gx=gx, rts=RTS_VRS, evalindex=evalindex, refindex=refindex)


class ITE(_RadialModel):
    """Single-input technical efficiency (Shephard input distance on one input).

    Only the **first** input listed in ``sent`` is adjusted; all other inputs
    and all outputs are held fixed.  Always input-oriented.
    """

    def __init__(self, data, sent, rts=RTS_CRS, evalindex=None, refindex=None):
        spec = parse_sent(sent)
        gx = [1] + [0] * (len(spec.x) - 1)
        gy = [0] * len(spec.y)
        super().__init__(data, sent, gy=gy, gx=gx, rts=rts, evalindex=evalindex, refindex=refindex)


class Hyperbolic(_PerDmuEnvelopmentModel):
    """Hyperbolic distance function: shrink all inputs and expand all outputs.

    The hyperbolic measure has no partial-input variant: every input is shrunk
    and every output expanded, so this class accepts neither ``gx``/``gy`` nor
    an ``orient`` option (``gx``/``gy`` are fixed to all ones internally).
    """

    def __init__(self, data, sent, rts=RTS_CRS, evalindex=None, refindex=None):
        spec = parse_sent(sent)
        gx = [1] * len(spec.x)
        gy = [1] * len(spec.y)
        super().__init__(data, sent, gy=gy, gx=gx, rts=rts, evalindex=evalindex, refindex=refindex)

    def _add_decision_variables(self, model, actual_index):
        bounds = (0.0, 1.0) if self.rts == RTS_CRS else (0.0, None)
        model.rho = Var(bounds=bounds, doc=f"hyperbolic efficiency for DMU {actual_index}")
        model.lamda = Var(model.R, bounds=(0.0, None), doc="intensity variables")

    def _add_objective(self, model):
        sense = minimize if self.rts == RTS_CRS else maximize
        model.objective = Objective(rule=lambda m: m.rho, sense=sense, doc="objective function")

    def _create_input_rule(self, current_dmu_range_index):
        # gx is fixed to all ones: every input is shrunk by rho (CRS) or by (1 - rho) (VRS).
        def input_rule(model, j):
            lhs = sum(model.lamda[r] * self.xref[r][j] for r in model.R)
            x0 = self.x[current_dmu_range_index][j]
            if self.rts == RTS_CRS:
                return lhs <= model.rho * x0
            return lhs <= x0 - model.rho * x0
        return input_rule

    def _create_output_rule(self, current_dmu_range_index):
        # gy is fixed to all ones: every output is expanded by (1 + rho) under VRS.
        def output_rule(model, k):
            lhs = sum(model.lamda[r] * self.yref[r][k] for r in model.R)
            y0 = self.y[current_dmu_range_index][k]
            if self.rts == RTS_CRS:
                return lhs >= y0
            return lhs >= y0 + model.rho * y0
        return output_rule

    # No ``_efficiency_columns`` override: ``hyperbolic`` reports the raw LP
    # value ``rho`` only.  Chapter 3 derives the efficiency from ``rho``
    # according to the model -- CRS: te = sqrt(rho); VRS: tei = 1 - rho on the
    # input side and teo = 1 / (1 + rho) on the output side -- because a single
    # efficiency column cannot express both sides of the VRS hyperbolic model.


def ccr(data, sent, orient="oo", evalindex=None, refindex=None):
    """CCR radial DEA (constant returns to scale).

    Returns :class:`CCR`, a self-contained radial model with ``rts=RTS_CRS``;
    the direction vectors ``gy``/``gx`` are derived from ``orient``.

    Parameters
    ----------
    data : pandas.DataFrame
        DEA data, one row per DMU.
    sent : str
        ``"inputs=outputs"`` declaration, e.g. ``"K L E=Y"``.
    orient : {"oo", "io"}, default "oo"
        Output ("oo") or input ("io") orientation.
    evalindex, refindex : str or None
        Optional evaluated/reference subset queries.
    """
    return CCR(data, sent, orient=orient, evalindex=evalindex, refindex=refindex)


def bcc(data, sent, orient="oo", evalindex=None, refindex=None):
    """BCC radial DEA (variable returns to scale).

    Returns :class:`BCC`, a self-contained radial model with ``rts=RTS_VRS``;
    the direction vectors ``gy``/``gx`` are derived from ``orient``.

    Parameters
    ----------
    data : pandas.DataFrame
        DEA data, one row per DMU.
    sent : str
        ``"inputs=outputs"`` declaration, e.g. ``"K L E=Y"``.
    orient : {"oo", "io"}, default "oo"
        Output ("oo") or input ("io") orientation.
    evalindex, refindex : str or None
        Optional evaluated/reference subset queries.
    """
    return BCC(data, sent, orient=orient, evalindex=evalindex, refindex=refindex)


def hyperbolic(data, sent, rts=RTS_CRS, evalindex=None, refindex=None):
    """Hyperbolic distance function: shrink all inputs and expand all outputs.

    Returns :class:`Hyperbolic`, a self-contained model with ``gy`` and ``gx``
    set to the all-inputs/all-outputs direction; there is no ``orient`` option.

    Parameters
    ----------
    data : pandas.DataFrame
        DEA data, one row per DMU.
    sent : str
        ``"inputs=outputs"`` declaration, e.g. ``"K L E=Y"``.
    rts : RTS_CRS or RTS_VRS, default RTS_CRS
        Returns-to-scale assumption.
    evalindex, refindex : str or None
        Optional evaluated/reference subset queries.
    """
    return Hyperbolic(data, sent, rts=rts, evalindex=evalindex, refindex=refindex)


def ite(data, sent, rts=RTS_CRS, evalindex=None, refindex=None):
    """Single-input technical efficiency: shrink only the first listed input.

    Returns :class:`ITE`, a self-contained radial model for the Shephard input
    distance function defined on a single input (e.g. energy or carbon); the direction
    vector ``gx`` marks only the **first** input in ``sent`` for radial
    shrinkage, while all other inputs and all outputs are held fixed.  There
    is no ``orient`` option (always input-oriented).

    Parameters
    ----------
    data : pandas.DataFrame
        DEA data, one row per DMU.
    sent : str
        ``"inputs=outputs"`` declaration, e.g. ``"E K L=Y"``.  The **first**
        input on the left of ``=`` is the one adjusted (here ``E``); the
        remaining inputs are held fixed.
    rts : RTS_CRS or RTS_VRS, default RTS_CRS
        Returns-to-scale assumption.
    evalindex, refindex : str or None
        Optional evaluated/reference subset queries.
    """
    return ITE(data, sent, rts=rts, evalindex=evalindex, refindex=refindex)


def ddf(data, sent, gy=[1], gx=[0], rts=RTS_CRS, evalindex=None, refindex=None):
    """Directional distance function: move along an explicit direction vector.

    Returns :class:`_DdfRho`, a self-contained directional-distance model that
    reports the raw LP value ``rho`` rather than an efficiency column; the
    efficiency follows from the chosen direction (see Chapter 3).  This is the
    only interface of the five that needs the direction vectors ``gy``/``gx``;
    ``ccr``/``bcc``/``hyperbolic``/``ite`` derive their direction internally
    and therefore do not accept them.

    Parameters
    ----------
    data : pandas.DataFrame
        DEA data, one row per DMU.
    sent : str
        ``"inputs=outputs"`` declaration, e.g. ``"K L E=Y"``.
    gy : list[int], default [1]
        Output direction vector; ``[0]`` leaves outputs unchanged.
    gx : list[int], default [0]
        Input direction vector; ``[0]`` leaves inputs unchanged.
    rts : RTS_CRS or RTS_VRS, default RTS_CRS
        Returns-to-scale assumption.
    evalindex, refindex : str or None
        Optional evaluated/reference subset queries.
    """
    return _DdfRho(data, sent, gy=gy, gx=gx, rts=rts, evalindex=evalindex, refindex=refindex)


# ---------------------------------------------------------------------------
# Malmquist interfaces (Chapter 4).  These are the ONLY places allowed to call
# the Chapter-3 factories ``ccr`` / ``bcc`` / ``ddf``; the shared engine lives
# in ``.mqcore`` (a leaf module) so that ``DEA.py`` never imports ``MQDEA.py``.
# ---------------------------------------------------------------------------


def _radial_distance(data, sent, orient, rts, evalindex, refindex, email, solver):
    """Efficiency value ``E`` (0 < E <= 1) for one (eval, ref) pair, via ccr/bcc.

    ``te`` is stored as-is for BOTH orientations: for ``orient="oo"`` it is the
    Shephard output distance ``D_o`` (eq-4-1); for ``orient="io"`` it is the
    Farrell input efficiency ``1/D_i`` required by eq-4-9/eq-4-10.  (The legacy
    ``MQDEA`` returned ``1/te`` for ``"io"``, which inverted the index.)
    """
    model = (
        ccr(data, sent, orient=orient, evalindex=evalindex, refindex=refindex)
        if rts == RTS_CRS
        else bcc(data, sent, orient=orient, evalindex=evalindex, refindex=refindex)
    )
    res = model.optimize(email, solver)
    out = res[["te"]].astype(float)
    out.columns = ["E"]
    return out


def _ddf_distance(data, sent, gy, gx, rts, evalindex, refindex, email, solver):
    """Directional-distance efficiency bundle(s), via ddf().

    ``ddf`` returns only ``rho``; the conversion to efficiency is done here at
    the call site, per Chapter 3 tbl-3-6:

    - input direction  -> ``E``  = ``1 - rho``
    - output direction -> ``E``  = ``1 / (1 + rho)``
    - both (hyperbolic) -> ``E_tei = 1 - rho`` and ``E_teo = 1 / (1 + rho)``
    """
    res = ddf(data, sent, gy=gy, gx=gx, rts=rts, evalindex=evalindex, refindex=refindex).optimize(email, solver)
    rho = res["rho"].astype(float)
    has_x = any(v > 0 for v in gx)
    has_y = any(v > 0 for v in gy)
    if not has_x and not has_y:
        raise ValueError("ddf direction vectors gx/gy cannot both be zero.")
    out = pd.DataFrame(index=res.index)
    if has_x and has_y:
        out["E_tei"] = 1.0 - rho
        out["E_teo"] = 1.0 / (1.0 + rho)
    elif has_x:
        out["E"] = 1.0 - rho
    else:
        out["E"] = 1.0 / (1.0 + rho)
    return out


def _ite_distance(data, sent, rts, evalindex, refindex, email, solver):
    """Energy-specific Shephard efficiency bundle from :func:`ite`.

    `ite` shrinks only the **first** input in `sent` (the energy) and holds the
    remaining inputs and all outputs fixed; its `te` column is the Farrell-type
    efficiency score theta.  Returns a one-column DataFrame (column ``"E"``)
    indexed by the evaluated rows, matching the contract expected by
    :func:`build_distance_frame`.
    """
    res = ite(data, sent, rts=rts, evalindex=evalindex, refindex=refindex).optimize(email, solver)
    out = pd.DataFrame(index=res.index)
    out["E"] = res["te"].astype(float)
    return out


def _make_model(data, id, year, tech, window, build, kind, orient, verbose):
    return MalmquistModel(
        data, id, year, tech=tech, window=window, build=build, kind=kind, orient=orient, verbose=verbose
    )


def mq(data, id, year, sent="inputvar=outputvar", orient="oo",
       rts=RTS_CRS, tech=CONTEMPORARY, window=3, verbose=False):
    """Radial Malmquist index (Shephard distance), built on :func:`ccr`/:func:`bcc`.

    Returns a :class:`~deabook.mqcore.MalmquistModel`; call ``.optimize(email,
    solver)`` to solve.  ``tech`` accepts ``CONTEMPORARY`` (default), ``WINDOW``,
    ``SEQUENTIAL`` and ``GLOBAL`` (Chapter 3 §3.4).
    """
    def build(work, years, t, w, year_col, email, solver):
        solve = lambda e, r, em, sv: _radial_distance(work, sent, orient, rts, e, r, em, sv)
        if t == GLOBAL:
            # Pastor-Lovell (2005) global Malmquist: mqpi = MEFFCH x MTECHCH.
            # MEFFCH is the *contemporary* efficiency change D^t/D^{t-1};
            # MTECHCH is the best-practice-gap change
            #   BPG_t / BPG_{t-1},  BPG_t = D^G(x^t)/D^t(x^t).
            glob = build_distance_frame(solve, work, years, GLOBAL, w, year_col, id_col=id, email=email, solver=solver, bundle="E")
            cont = build_distance_frame(solve, work, years, CONTEMPORARY, w, year_col, id_col=id, email=email, solver=solver, bundle="E")
            bpg_frame = pd.DataFrame({id: work[id], "BPG": glob["D11"] / cont["D11"]}, index=work.index)
            bpg_frame.attrs = cont.attrs.copy()
            result = pd.DataFrame(index=work.index)
            result["D11"] = glob["D11"]
            result["mqpi"] = glob["D11"] / _shift(glob, id, "D11")
            result["MEFFCH"] = cont["D11"] / _shift(cont, id, "D11")
            result["MTECHCH"] = bpg_frame["BPG"] / _shift(bpg_frame, id, "BPG")
            return result

        frame = build_distance_frame(solve, work, years, t, w, year_col, id_col=id, email=email, solver=solver, bundle="E")
        return pd.concat([frame[["D11", "D12", "D21"]], two_factor(frame, id, suffix="")], axis=1)

    return _make_model(data, id, year, tech, window, build, "mq", orient, verbose)


def mqddf(data, id, year, sent="inputvar=outputvar", gy=[1], gx=[0],
          rts=RTS_CRS, tech=CONTEMPORARY, window=3, verbose=False):
    """DDF Malmquist index, built on :func:`ddf` with the ``rho`` -> efficiency
    conversion of Chapter 3 tbl-3-6.  Same ``tech`` values as :func:`mq`.
    """
    def build(work, years, t, w, year_col, email, solver):
        solve = lambda e, r, em, sv: _ddf_distance(work, sent, gy, gx, rts, e, r, em, sv)
        hyperbolic = any(v > 0 for v in gx) and any(v > 0 for v in gy)
        if not hyperbolic:
            if t == GLOBAL:
                glob = build_distance_frame(solve, work, years, GLOBAL, w, year_col, id_col=id, email=email, solver=solver, bundle="E")
                cont = build_distance_frame(solve, work, years, CONTEMPORARY, w, year_col, id_col=id, email=email, solver=solver, bundle="E")
                bpg_frame = pd.DataFrame({id: work[id], "BPG": glob["D11"] / cont["D11"]}, index=work.index)
                bpg_frame.attrs = cont.attrs.copy()
                result = pd.DataFrame(index=work.index)
                result["D11"] = glob["D11"]
                result["mqpi"] = glob["D11"] / _shift(glob, id, "D11")
                result["MEFFCH"] = cont["D11"] / _shift(cont, id, "D11")
                result["MTECHCH"] = bpg_frame["BPG"] / _shift(bpg_frame, id, "BPG")
                return result
            frame = build_distance_frame(solve, work, years, t, w, year_col, id_col=id, email=email, solver=solver, bundle="E")
            return pd.concat([frame[["D11", "D12", "D21"]], two_factor(frame, id, suffix="")], axis=1)
        out = pd.DataFrame(index=work.index)
        for bundle, suffix in (("E_tei", "_tei"), ("E_teo", "_teo")):
            frame = build_distance_frame(solve, work, years, t, w, year_col, id_col=id, email=email, solver=solver, bundle=bundle)
            dcols = frame[["D11", "D12", "D21"]].rename(columns=lambda c: c + suffix)
            result = two_factor(frame, id, suffix=suffix)
            out = pd.concat([out, dcols, result], axis=1)
        return out

    return _make_model(data, id, year, tech, window, build, "mqddf", None, verbose)


def ep(data, id, year, sent="E K L=Y", rts=RTS_CRS,
       tech=CONTEMPORARY, window=3, verbose=False):
    """能源专项谢泼德口径 Malmquist 指数，基于 :func:`ite`（单投入径向距离函数）。

    与 :func:`mq` 结构相同，但底层距离函数由 `ccr`/`bcc` 径向距离替换为
    :func:`ite`：只缩减 `sent` 中的第一个投入（能源），其余投入与产出固定。
    返回 :class:`~deabook.mqcore.MalmquistModel`，调用 ``.optimize(email,
    solver)`` 求解；``tech`` 取值同 :func:`mq`（`CONTEMPORARY`、``WINDOW``、
    ``SEQUENTIAL``、`GLOBAL`）。
    """
    def build(work, years, t, w, year_col, email, solver):
        solve = lambda e, r, em, sv: _ite_distance(work, sent, rts, e, r, em, sv)
        if t == GLOBAL:
            glob = build_distance_frame(solve, work, years, GLOBAL, w, year_col, id_col=id, email=email, solver=solver, bundle="E")
            cont = build_distance_frame(solve, work, years, CONTEMPORARY, w, year_col, id_col=id, email=email, solver=solver, bundle="E")
            bpg_frame = pd.DataFrame({id: work[id], "BPG": glob["D11"] / cont["D11"]}, index=work.index)
            bpg_frame.attrs = cont.attrs.copy()
            result = pd.DataFrame(index=work.index)
            result["D11"] = glob["D11"]
            result["mqpi"] = glob["D11"] / _shift(glob, id, "D11")
            result["MEFFCH"] = cont["D11"] / _shift(cont, id, "D11")
            result["MTECHCH"] = bpg_frame["BPG"] / _shift(bpg_frame, id, "BPG")
            return result
        frame = build_distance_frame(solve, work, years, t, w, year_col, id_col=id, email=email, solver=solver, bundle="E")
        return pd.concat([frame[["D11", "D12", "D21"]], two_factor(frame, id, suffix="")], axis=1)

    return _make_model(data, id, year, tech, window, build, "ep", None, verbose)


def epddf(data, id, year, sent="K L E=Y", gy=[0], gx=[0, 0, 1],
          rts=RTS_CRS, tech=CONTEMPORARY, window=3, verbose=False):
    """能源专项方向距离口径 Malmquist 指数，基于 :func:`ddf`（能源投入方向）。

    与 :func:`mqddf` 结构相同，但方向向量固定为能源投入方向 `gx=[0,0,1]`、
    `gy=[0]`（只缩减能源、其余不变），底层距离函数由 :func:`ddf` 给出并按
    $EE=1-\\rho$ 转换为效率。返回
    :class:`~deabook.mqcore.MalmquistModel`，调用 ``.optimize(email, solver)``
    求解。
    """
    def build(work, years, t, w, year_col, email, solver):
        solve = lambda e, r, em, sv: _ddf_distance(work, sent, gy, gx, rts, e, r, em, sv)
        if t == GLOBAL:
            glob = build_distance_frame(solve, work, years, GLOBAL, w, year_col, id_col=id, email=email, solver=solver, bundle="E")
            cont = build_distance_frame(solve, work, years, CONTEMPORARY, w, year_col, id_col=id, email=email, solver=solver, bundle="E")
            bpg_frame = pd.DataFrame({id: work[id], "BPG": glob["D11"] / cont["D11"]}, index=work.index)
            bpg_frame.attrs = cont.attrs.copy()
            result = pd.DataFrame(index=work.index)
            result["D11"] = glob["D11"]
            result["mqpi"] = glob["D11"] / _shift(glob, id, "D11")
            result["MEFFCH"] = cont["D11"] / _shift(cont, id, "D11")
            result["MTECHCH"] = bpg_frame["BPG"] / _shift(bpg_frame, id, "BPG")
            return result
        frame = build_distance_frame(solve, work, years, t, w, year_col, id_col=id, email=email, solver=solver, bundle="E")
        return pd.concat([frame[["D11", "D12", "D21"]], two_factor(frame, id, suffix="")], axis=1)

    return _make_model(data, id, year, tech, window, build, "epddf", None, verbose)


# `DEA.ep.ddf` 作为 `ep` 的 DDF 变体（能源投入方向），与 `gmqdecom.ddf` 等命名保持一致。
ep.ddf = epddf


def _two_tech_build(sent, orient, id, decompose):
    """Build function for the CRS+VRS decomposition interfaces."""

    def build(work, years, t, w, year_col, email, solver):
        def frame_for(rts):
            solve = lambda e, r, em, sv: _radial_distance(work, sent, orient, rts, e, r, em, sv)
            return build_distance_frame(solve, work, years, t, w, year_col, id_col=id, email=email, solver=solver, bundle="E")

        return decompose(frame_for(RTS_VRS), frame_for(RTS_CRS), id)

    return build


def mqfgnz(data, id, year, sent="inputvar=outputvar", orient="oo",
           tech=CONTEMPORARY, window=3, verbose=False):
    """FGNZ (1994) three-factor Malmquist decomposition, via ccr AND bcc.

    Runs the radial distances under both CRS and VRS.  Same ``tech`` values as
    :func:`mq`.
    """
    build = _two_tech_build(sent, orient, id, fgnz_decompose)
    return _make_model(data, id, year, tech, window, build, "mqfgnz", orient, verbose)


def mqrd(data, id, year, sent="inputvar=outputvar", orient="oo",
         tech=CONTEMPORARY, window=3, verbose=False):
    """Ray & Desli (1997) decomposition, via ccr AND bcc.

    ``EFFCH`` is VRS efficiency change and ``TECHCH_RD`` uses VRS cross-period
    distances. ``SCALECH_RD`` completes the three-factor CRS Malmquist split:
    ``MQ = EFFCH * TECHCH_RD * SCALECH_RD``.
    """
    build = _two_tech_build(sent, orient, id, rd_decompose)
    return _make_model(data, id, year, tech, window, build, "mqrd", orient, verbose)


def _pl_global_build(sent, orient, id, distance):
    """Build the proper Pastor-Lovell three-factor decomposition under GLOBAL.

    ``distance`` is a callable ``(work, rts, evalindex, refindex, email, solver)``
    returning the efficiency frame; three distance sets are solved:
    contemporary VRS, contemporary CRS and global CRS.
    """

    def build(work, years, t, w, year_col, email, solver):
        def frame_for(rts, tech_value):
            solve = lambda e, r, em, sv: distance(work, rts, e, r, em, sv)
            return build_distance_frame(solve, work, years, tech_value, w, year_col, id_col=id, email=email, solver=solver, bundle="E")

        return pl_global_decompose(
            frame_for(RTS_VRS, CONTEMPORARY),
            frame_for(RTS_CRS, CONTEMPORARY),
            frame_for(RTS_CRS, GLOBAL),
            id,
        )

    return build


def gmqdecom(data, id, year, sent="inputvar=outputvar", orient="oo",
             tech=GLOBAL, window=3, verbose=False):
    """Pastor & Lovell (2005) global-frontier three-factor decomposition, via ccr AND bcc.

    Under ``tech=GLOBAL`` (the default) the global Malmquist index is decomposed
    the Pastor-Lovell way: ``PEFFCH`` (contemporary VRS efficiency change) x
    ``SEFFCH`` (contemporary scale-efficiency change) x ``TECHCH`` (best-practice-
    gap change), with ``MQ = mqpi = PEFFCH * SEFFCH * TECHCH``.  The other three
    ``tech`` values give the FGNZ-style three-factor split against that
    reference set.
    """

    def radial_distance(work, rts, e, r, em, sv):
        return _radial_distance(work, sent, orient, rts, e, r, em, sv)

    def build(work, years, t, w, year_col, email, solver):
        if t == GLOBAL:
            return _pl_global_build(sent, orient, id, radial_distance)(work, years, t, w, year_col, email, solver)
        return _two_tech_build(sent, orient, id, global_decompose)(work, years, t, w, year_col, email, solver)

    return _make_model(data, id, year, tech, window, build, "gmqdecom", orient, verbose)


def _ddf_two_tech_build(sent, gy, gx, id, decompose):
    """Build function for the DDF-based CRS+VRS decomposition interfaces."""

    def build(work, years, t, w, year_col, email, solver):
        def frame_for(rts):
            solve = lambda e, r, em, sv: _ddf_distance(work, sent, gy, gx, rts, e, r, em, sv)
            return build_distance_frame(solve, work, years, t, w, year_col, id_col=id, email=email, solver=solver, bundle="E")

        return decompose(frame_for(RTS_VRS), frame_for(RTS_CRS), id)

    return build


def _check_single_direction(gy, gx, name):
    if any(v > 0 for v in gx) and any(v > 0 for v in gy):
        raise ValueError(
            f"{name} supports a single direction (gx or gy active, not both); "
            "hyperbolic directions are not part of the Malmquist decomposition."
        )


def _gmqdecom_ddf(data, id, year, sent="inputvar=outputvar", gy=[1], gx=[0],
                  tech=GLOBAL, window=3, verbose=False):
    """DDF-based global Malmquist three-factor decomposition.

    The DDF variant of :func:`gmqdecom`: it decomposes the global Malmquist
    index computed from **directional distance functions** (via :func:`ddf`,
    with the ``rho`` -> efficiency conversion of Chapter 3 tbl-3-6) into the
    same three factors -- ``PEFFCH`` (contemporary VRS pure efficiency change),
    ``SEFFCH`` (scale efficiency change), ``TECHCH`` (best-practice-gap change)
    -- plus ``MQ``/``MEFFCH``/``MTECHCH`` and ``mqpi``.

    Only a single direction is supported (``gx`` or ``gy`` active, not both).
    Defaults to ``tech=GLOBAL`` (the Pastor-Lovell frontier).
    """
    _check_single_direction(gy, gx, "gmqdecom.ddf")

    def ddf_distance(work, rts, e, r, em, sv):
        return _ddf_distance(work, sent, gy, gx, rts, e, r, em, sv)

    def build(work, years, t, w, year_col, email, solver):
        if t == GLOBAL:
            return _pl_global_build(sent, None, id, ddf_distance)(work, years, t, w, year_col, email, solver)
        return _ddf_two_tech_build(sent, gy, gx, id, global_decompose)(work, years, t, w, year_col, email, solver)

    return _make_model(data, id, year, tech, window, build, "gmqdecom.ddf", None, verbose)


# ``DEA.gmqdecom.ddf(...)`` -- the DDF variant of gmqdecom.
gmqdecom.ddf = _gmqdecom_ddf


def _mqfgnz_ddf(data, id, year, sent="inputvar=outputvar", gy=[1], gx=[0],
                tech=CONTEMPORARY, window=3, verbose=False):
    """DDF-based FGNZ (1994) three-factor Malmquist decomposition.

    The DDF variant of :func:`mqfgnz`: the CRS and VRS distances are computed
    from directional distance functions (via :func:`ddf`, with the ``rho`` ->
    efficiency conversion of Chapter 3 tbl-3-6); the decomposition itself is
    identical.  Only a single direction is supported.
    """
    _check_single_direction(gy, gx, "mqfgnz.ddf")
    build = _ddf_two_tech_build(sent, gy, gx, id, fgnz_decompose)
    return _make_model(data, id, year, tech, window, build, "mqfgnz.ddf", None, verbose)


# ``DEA.mqfgnz.ddf(...)`` -- the DDF variant of mqfgnz.
mqfgnz.ddf = _mqfgnz_ddf


def _mqrd_ddf(data, id, year, sent="inputvar=outputvar", gy=[1], gx=[0],
              tech=CONTEMPORARY, window=3, verbose=False):
    """DDF-based Ray & Desli (1997) decomposition.

    The DDF variant of :func:`mqrd`: ``EFFCH`` equals ``mqfgnz.ddf``'s
    ``PEFFCH``; ``TECHCH_RD = MQ / EFFCH`` so the Malmquist index equals
    ``mqfgnz.ddf``'s ``MQ``.  Only a single direction is supported.
    """
    _check_single_direction(gy, gx, "mqrd.ddf")
    build = _ddf_two_tech_build(sent, gy, gx, id, rd_decompose)
    return _make_model(data, id, year, tech, window, build, "mqrd.ddf", None, verbose)


# ``DEA.mqrd.ddf(...)`` -- the DDF variant of mqrd.
mqrd.ddf = _mqrd_ddf


# --- static slacks-based measure (SBM) factory ---------------------------------
def _sbm_sent_has_bad(sent):
    """Return ``True`` if ``sent`` declares an undesirable-output component.

    The package ``sent`` convention is ``"inputs = outputs"`` or
    ``"inputs = outputs : undesirable_outputs"``; an undesirable part exists iff
    the sub-string after the first ``=`` contains a ``:``.
    """
    return len(parse_sent(sent).b) > 0


def sbm(data, sent="inputvar=outputvar", gx=None, gy=None, gb=None,
        rts=RTS_VRS1, evalindex=None, refindex=None):
    """Static slacks-based measure (SBM) *without* undesirable outputs (§8.3).

    Builds the Tone (2001) SBM efficiency model on the private
    :class:`~deabook.DEA._SBM` engine and returns the (unsolved) model object;
    call ``.optimize(email, solver)`` to solve and obtain the efficiency score
    ``rho`` together with the per-variable ``slack_x_*`` / ``slack_y_*`` columns.
    ``sent`` must declare only inputs and desirable outputs (no ``:``-separated
    undesirable part) -- use :func:`DEAweak.sbmu` for the undesirable-output
    version.

    Args:
        data (pandas.DataFrame): Input-output data.
        sent (str): ``"inputs = outputs"`` variable expression, e.g. ``"K L=Y"``.
        gx (list, optional): Input slack direction (1 = slack enters objective).
        gy (list, optional): Desirable-output slack direction.
        gb (list, optional): Undesirable-output slack direction (ignored here,
            since ``sbm`` models no undesirable output).
        rts (str): ``RTS_CRS`` or ``RTS_VRS1``.
        evalindex (str, optional): Evaluated-sample filter.
        refindex (str, optional): Reference-sample filter.

    Returns:
        deabook.DEA._SBM: An unsolved SBM model; call ``.optimize()`` on it.

    Example:
        >>> m = sbm(data, sent="K L=Y", gy=[1], gx=[1, 1], rts=RTS_VRS1)
        >>> res = m.optimize(OPT_LOCAL, "highs")
    """
    if _sbm_sent_has_bad(sent):
        raise ValueError(
            "DEA.sbm is for SBM *without* undesirable outputs, but the sent "
            f"{sent!r} declares an undesirable-output component. Use DEAweak.sbmu "
            "for the undesirable-output SBM (§8.4)."
        )
    kwargs = {}
    if gx is not None:
        kwargs["gx"] = gx
    if gy is not None:
        kwargs["gy"] = gy
    if gb is not None:
        kwargs["gb"] = gb
    return _SBM(data, sent, rts=rts, evalindex=evalindex, refindex=refindex, **kwargs)


def sddf(data, sent="inputvar=outputvar", gx=None, gy=None,
         rts=RTS_VRS1, evalindex=None, refindex=None):
    """Slack-based directional distance function (SD), Färe & Grosskopf (2010).

    Maximizes the *sum of absolute slacks* -- an independent slack per input and
    per desirable output (§8.5.1) -- which is the translation-invariant
    generalization of the Tone (2001) SBM measure and is strictly more
    discriminating than the directional, hyperbolic and Farrell measures.

    Returns the (unsolved) :class:`~deabook.DEA._SDDF` model; call
    ``.optimize(email, solver)`` to obtain ``objective_value`` (the SD
    inefficiency = the sum of slacks) together with the per-variable
    ``slack_x_*`` / ``slack_y_*`` columns.  ``sent`` must declare only inputs
    and desirable outputs (no ``:``-separated undesirable part) -- use
    :func:`DEAweak.sddf` for the undesirable-output extension of Pham & Zelenyuk
    (2018, §8.5.2).

    Args:
        data (pandas.DataFrame): Input-output data.
        sent (str): ``"inputs = outputs"`` variable expression, e.g. ``"K L=Y"``.
        gx (list, optional): Input slack direction; a positive element marks the
            input as eligible for the input-excess slack.  Defaults to all inputs.
        gy (list, optional): Desirable-output slack direction.  Defaults to all
            outputs.
        rts (str): ``RTS_CRS`` or ``RTS_VRS1``.
        evalindex (str, optional): Evaluated-sample filter.
        refindex (str, optional): Reference-sample filter.

    Returns:
        deabook.DEA._SDDF: An unsolved SDDF model; call ``.optimize()`` on it.

    Example:
        >>> m = sddf(data, sent="K L=Y", rts=RTS_VRS1)
        >>> res = m.optimize(OPT_LOCAL, "highs")
    """
    if _sbm_sent_has_bad(sent):
        raise ValueError(
            "DEA.sddf is the Färe-Grosskopf (2010) SD model *without* "
            f"undesirable outputs, but the sent {sent!r} declares an "
            "undesirable-output component. Use DEAweak.sddf for the "
            "Pham-Zelenyuk (2018) bad-output extension (§8.5.2)."
        )
    kwargs = {}
    if gx is not None:
        kwargs["gx"] = gx
    if gy is not None:
        kwargs["gy"] = gy
    return _SDDF(data, sent, rts=rts, evalindex=evalindex, refindex=refindex, **kwargs)


# ---------------------------------------------------------------------------
# Private SBM engine.
#
# Folded in from the now-deleted ``deabook.SBM`` module so that the public API
# is reached only through :func:`sbm` / :func:`DEAweak.sbmu` (the standalone
# ``SBM`` class is intentionally no longer importable).
# ---------------------------------------------------------------------------

class _SBM:
    """Slacks-based measure DEA model.

    Args:
        data (pandas.DataFrame): Input-output data.
        sent (str): Variable expression using the package convention
            ``"inputs = outputs"`` or ``"inputs = outputs : undesirable_outputs"``.
            Examples: ``"K L = Y"`` and ``"K L E = Y : CO2"``.
        gy (list): Desirable-output direction vector. A positive element means
            that the corresponding desirable output is adjusted upward.
        gx (list): Input direction vector. A positive element means that the
            corresponding input is adjusted downward.
        gb (list): Undesirable-output direction vector. A positive element
            means that the corresponding undesirable output is adjusted downward.
        rts (str): ``RTS_CRS`` or ``RTS_VRS1``.
        evalindex (str, optional): Evaluated-sample filter such as
            ``"Year=[2010]"``. Defaults to all observations.
        refindex (str, optional): Reference-sample filter such as
            ``"Year=[2010]"``. Defaults to all observations.
    """

    def __init__(
        self,
        data,
        sent,
        gy=[1],
        gx=[0],
        gb=[1],
        rts=RTS_VRS1,
        evalindex=None,
        refindex=None,
    ):
        self.data = data
        self.sent = sent
        self.rts = rts
        self.evalindex = evalindex
        self.refindex = refindex

        if self.rts not in {RTS_CRS, RTS_VRS1}:
            raise ValueError("rts must be RTS_CRS or RTS_VRS1.")

        self.formula = parse_sent(sent)
        self.xcol, self.ycol, self.bcol = self.formula.x, self.formula.y, self.formula.b
        self.direction_spec = build_direction_spec(
            gx,
            gy,
            gb if self.bcol else [],
            num_inputs=len(self.xcol),
            num_outputs=len(self.ycol),
            num_bad_outputs=len(self.bcol),
            allow_bad_direction=True,
        )
        self.gx = self.direction_spec.gx
        self.gy = self.direction_spec.gy
        self.gb = self.direction_spec.gb

        self.input_oriented = sum(self.gx) >= 1 and sum(self.gy) == 0 and sum(self.gb) == 0
        self.output_oriented = sum(self.gy) >= 1 and sum(self.gx) == 0 and sum(self.gb) == 0
        self.unoutput_oriented = sum(self.gb) >= 1 and sum(self.gx) == 0 and sum(self.gy) == 0
        self.hyper_orientedyx = sum(self.gx) >= 1 and sum(self.gy) >= 1 and sum(self.gb) == 0
        self.hyper_orientedyb = sum(self.gb) >= 1 and sum(self.gy) >= 1 and sum(self.gx) == 0
        self.hyper_orientedxb = sum(self.gb) >= 1 and sum(self.gx) >= 1 and sum(self.gy) == 0
        self.hyper_orientedyxb = sum(self.gb) >= 1 and sum(self.gx) >= 1 and sum(self.gy) >= 1
        self.direction = _direction_name(self.gx, self.gy, self.gb)

        prepared = prepare_dea_data(
            data,
            self.formula,
            eval_query=evalindex,
            ref_query=refindex,
            input_side="left",
            require_nonnegative_b=False,
        )
        self.evaluated_data_index = prepared.evaluated_index
        self.reference_data_index = prepared.reference_index
        self.x = prepared.x
        self.y = prepared.y
        self.xref = prepared.xref
        self.yref = prepared.yref
        self.b = prepared.b
        self.bref = prepared.bref
        if self.bcol and (np.any(self.b <= 0) or np.any(self.bref <= 0)):
            raise ValueError("Undesirable outputs must be strictly positive for SBM ratios.")

        self.evaluated_indices_range = range(len(self.evaluated_data_index))
        self.reference_indices_range = range(len(self.reference_data_index))
        self.num_inputs = len(self.xcol)
        self.num_outputs = len(self.ycol)
        self.num_unoutputs = len(self.bcol)
        self.x_active = [j for j, item in enumerate(self.gx) if item > 0]
        self.y_active = [k for k, item in enumerate(self.gy) if item > 0]
        self.b_active = [h for h, item in enumerate(self.gb) if item > 0]

        self.__modeldict = {}
        self.results = None
        self.lamda = None

        for i_range in self.evaluated_indices_range:
            actual_index = self.evaluated_data_index[i_range]
            model = self.__create_model(i_range, actual_index)
            self.__modeldict[actual_index] = model

    def __create_model(self, current_dmu_range_index, actual_index):
        model = ConcreteModel()
        model.R = Set(initialize=self.reference_indices_range)
        model.J = Set(initialize=range(self.num_inputs))
        model.K = Set(initialize=range(self.num_outputs))
        model.L = Set(initialize=range(self.num_unoutputs))
        model.AX = Set(initialize=self.x_active)
        model.AY = Set(initialize=self.y_active)
        model.AB = Set(initialize=self.b_active)
        model.lamda = Var(model.R, bounds=(0.0, None), doc="intensity variables")

        x0 = self.x[current_dmu_range_index]
        y0 = self.y[current_dmu_range_index]
        b0 = None if self.b is None else self.b[current_dmu_range_index]

        if self.x_active and (self.y_active or self.b_active):
            model.t = Var(bounds=(1e-9, None), doc=f"Charnes-Cooper variable for DMU {actual_index}")
            model.sx = Var(model.AX, bounds=(0.0, None), doc="transformed input slacks")
            model.sy = Var(model.AY, bounds=(0.0, None), doc="transformed desirable-output slacks")
            model.sb = Var(model.AB, bounds=(0.0, None), doc="transformed undesirable-output slacks")
            model.objective = Objective(rule=self.__ratio_objective_rule(x0), sense=minimize)
            model.cctrans = Constraint(rule=self.__cctrans_rule(y0, b0))
            model.input = Constraint(model.J, rule=self.__ratio_input_rule(current_dmu_range_index))
            model.output = Constraint(model.K, rule=self.__ratio_output_rule(current_dmu_range_index))
            model.unoutput = Constraint(model.L, rule=self.__ratio_unoutput_rule(current_dmu_range_index))
            if self.rts == RTS_VRS1:
                model.vrs = Constraint(rule=self.__ratio_vrs_rule())
            model._sbm_model_type = "ratio"
        elif self.x_active:
            model.sx = Var(model.AX, bounds=(0.0, None), doc="input slacks")
            model.objective = Objective(rule=self.__input_only_objective_rule(x0), sense=minimize)
            model.input = Constraint(model.J, rule=self.__input_only_input_rule(current_dmu_range_index))
            model.output = Constraint(model.K, rule=self.__input_only_output_rule(current_dmu_range_index))
            model.unoutput = Constraint(model.L, rule=self.__input_only_unoutput_rule(current_dmu_range_index))
            if self.rts == RTS_VRS1:
                model.vrs = Constraint(rule=self.__standard_vrs_rule())
            model._sbm_model_type = "input_only"
        else:
            model.sy = Var(model.AY, bounds=(0.0, None), doc="desirable-output slacks")
            model.sb = Var(model.AB, bounds=(0.0, None), doc="undesirable-output slacks")
            model.objective = Objective(rule=self.__output_bad_only_objective_rule(y0, b0), sense=maximize)
            model.input = Constraint(model.J, rule=self.__output_bad_only_input_rule(current_dmu_range_index))
            model.output = Constraint(model.K, rule=self.__output_bad_only_output_rule(current_dmu_range_index))
            model.unoutput = Constraint(model.L, rule=self.__output_bad_only_unoutput_rule(current_dmu_range_index))
            if self.rts == RTS_VRS1:
                model.vrs = Constraint(rule=self.__standard_vrs_rule())
            model._sbm_model_type = "output_bad_only"

        return model

    def __ratio_objective_rule(self, x0):
        def objective_rule(model):
            return model.t - sum(model.sx[j] / x0[j] for j in model.AX) / len(self.x_active)
        return objective_rule

    def __cctrans_rule(self, y0, b0):
        denom_count = len(self.y_active) + len(self.b_active)

        def cctrans_rule(model):
            y_term = sum(model.sy[k] / y0[k] for k in model.AY)
            b_term = sum(model.sb[h] / b0[h] for h in model.AB) if self.bcol else 0.0
            return model.t + (y_term + b_term) / denom_count == 1
        return cctrans_rule

    def __ratio_input_rule(self, current_dmu_range_index):
        def input_rule(model, j):
            lhs = sum(model.lamda[r] * self.xref[r][j] for r in model.R)
            rhs = model.t * self.x[current_dmu_range_index][j]
            if j in self.x_active:
                return lhs + model.sx[j] == rhs
            return lhs <= rhs
        return input_rule

    def __ratio_output_rule(self, current_dmu_range_index):
        def output_rule(model, k):
            lhs = sum(model.lamda[r] * self.yref[r][k] for r in model.R)
            rhs = model.t * self.y[current_dmu_range_index][k]
            if k in self.y_active:
                return lhs - model.sy[k] == rhs
            return lhs >= rhs
        return output_rule

    def __ratio_unoutput_rule(self, current_dmu_range_index):
        def unoutput_rule(model, h):
            lhs = sum(model.lamda[r] * self.bref[r][h] for r in model.R)
            rhs = model.t * self.b[current_dmu_range_index][h]
            if h in self.b_active:
                return lhs + model.sb[h] == rhs
            return lhs <= rhs
        return unoutput_rule

    def __ratio_vrs_rule(self):
        def vrs_rule(model):
            return sum(model.lamda[r] for r in model.R) == model.t
        return vrs_rule

    def __input_only_objective_rule(self, x0):
        def objective_rule(model):
            return 1 - sum(model.sx[j] / x0[j] for j in model.AX) / len(self.x_active)
        return objective_rule

    def __input_only_input_rule(self, current_dmu_range_index):
        def input_rule(model, j):
            lhs = sum(model.lamda[r] * self.xref[r][j] for r in model.R)
            rhs = self.x[current_dmu_range_index][j]
            if j in self.x_active:
                return lhs + model.sx[j] == rhs
            return lhs <= rhs
        return input_rule

    def __input_only_output_rule(self, current_dmu_range_index):
        def output_rule(model, k):
            return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= self.y[current_dmu_range_index][k]
        return output_rule

    def __input_only_unoutput_rule(self, current_dmu_range_index):
        def unoutput_rule(model, h):
            return sum(model.lamda[r] * self.bref[r][h] for r in model.R) <= self.b[current_dmu_range_index][h]
        return unoutput_rule

    def __output_bad_only_objective_rule(self, y0, b0):
        denom_count = len(self.y_active) + len(self.b_active)

        def objective_rule(model):
            y_term = sum(model.sy[k] / y0[k] for k in model.AY)
            b_term = sum(model.sb[h] / b0[h] for h in model.AB) if self.bcol else 0.0
            return 1 + (y_term + b_term) / denom_count
        return objective_rule

    def __output_bad_only_input_rule(self, current_dmu_range_index):
        def input_rule(model, j):
            return sum(model.lamda[r] * self.xref[r][j] for r in model.R) <= self.x[current_dmu_range_index][j]
        return input_rule

    def __output_bad_only_output_rule(self, current_dmu_range_index):
        def output_rule(model, k):
            lhs = sum(model.lamda[r] * self.yref[r][k] for r in model.R)
            rhs = self.y[current_dmu_range_index][k]
            if k in self.y_active:
                return lhs - model.sy[k] == rhs
            return lhs >= rhs
        return output_rule

    def __output_bad_only_unoutput_rule(self, current_dmu_range_index):
        def unoutput_rule(model, h):
            lhs = sum(model.lamda[r] * self.bref[r][h] for r in model.R)
            rhs = self.b[current_dmu_range_index][h]
            if h in self.b_active:
                return lhs + model.sb[h] == rhs
            return lhs <= rhs
        return unoutput_rule

    def __standard_vrs_rule(self):
        def vrs_rule(model):
            return sum(model.lamda[r] for r in model.R) == 1
        return vrs_rule

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Solve each DMU's Pyomo model.

        Args:
            email (str): NEOS email. Use ``OPT_LOCAL`` for local solving.
            solver (str): Pyomo solver name. Defaults to the package default,
                which is handled by the shared solver configuration.

        Returns:
            pandas.DataFrame: Solver status, direction name, SBM efficiency
            score ``rho`` and variable-specific slack columns.
        """

        rows = []
        lamda_rows = []
        config = SolverConfig(solver=solver, email=email, tee=False)

        for actual_index, model in self.__modeldict.items():
            try:
                outcome = solve_pyomo_model(model, config=config, strict=False)
                optimization_status = outcome.legacy_status
            except Exception as exc:
                print(f"Error optimizing DMU {actual_index}: {exc}")
                optimization_status = "solve_error"

            if optimization_status in {"ok", "optimal"}:
                try:
                    row, lamda = self.__extract_solution(model, optimization_status)
                except Exception as exc:
                    print(f"Warning: could not retrieve SBM results for DMU {actual_index}: {exc}")
                    row = self.__empty_result_row("result_error")
                    lamda = np.full(len(self.reference_data_index), np.nan)
            else:
                row = self.__empty_result_row(optimization_status)
                lamda = np.full(len(self.reference_data_index), np.nan)

            rows.append(pd.Series(row, name=actual_index))
            lamda_rows.append(pd.Series(lamda, index=self.reference_data_index, name=actual_index))

        self.results = pd.DataFrame(rows)
        self.lamda = pd.DataFrame(lamda_rows)
        return self.results

    def __extract_solution(self, model, optimization_status):
        sx = np.zeros(len(self.xcol))
        sy = np.zeros(len(self.ycol))
        sb = np.zeros(len(self.bcol))
        lamda = np.asarray([value(model.lamda[r]) for r in model.R], dtype=float)
        objective_value = value(model.objective)
        row = {
            "optimization_status": optimization_status,
            "direction": self.direction,
            "objective_value": objective_value,
        }

        if model._sbm_model_type == "ratio":
            t_value = value(model.t)
            lamda = clean_small_values(lamda / t_value)
            for j in model.AX:
                sx[j] = value(model.sx[j]) / t_value
            for k in model.AY:
                sy[k] = value(model.sy[k]) / t_value
            for h in model.AB:
                sb[h] = value(model.sb[h]) / t_value
            row["rho"] = _clean_efficiency(objective_value)
            row["t"] = t_value
        elif model._sbm_model_type == "input_only":
            lamda = clean_small_values(lamda)
            for j in model.AX:
                sx[j] = value(model.sx[j])
            row["rho"] = _clean_efficiency(objective_value)
        elif model._sbm_model_type == "output_bad_only":
            lamda = clean_small_values(lamda)
            for k in model.AY:
                sy[k] = value(model.sy[k])
            for h in model.AB:
                sb[h] = value(model.sb[h])
            phi = objective_value
            row["rho"] = _clean_efficiency(np.nan if phi == 0 else 1 / phi)
            row["phi"] = phi

        sx = clean_small_values(sx)
        sy = clean_small_values(sy)
        sb = clean_small_values(sb)
        for j, name in enumerate(self.xcol):
            row[f"slack_x_{name}"] = sx[j]
        for k, name in enumerate(self.ycol):
            row[f"slack_y_{name}"] = sy[k]
        for h, name in enumerate(self.bcol):
            row[f"slack_b_{name}"] = sb[h]
        return row, lamda

    def __empty_result_row(self, status):
        row = {
            "optimization_status": status,
            "direction": self.direction,
            "objective_value": np.nan,
            "rho": np.nan,
        }
        if self.x_active and (self.y_active or self.b_active):
            row["t"] = np.nan
        if not self.x_active:
            row["phi"] = np.nan
        for name in self.xcol:
            row[f"slack_x_{name}"] = np.nan
        for name in self.ycol:
            row[f"slack_y_{name}"] = np.nan
        for name in self.bcol:
            row[f"slack_b_{name}"] = np.nan
        return row

    def get_status(self):
        """Return the solver status for each evaluated DMU."""
        _assert_optimized(self.results)
        return self.results["optimization_status"]

    def get_rho(self):
        """Return SBM efficiency scores."""
        _assert_optimized(self.results)
        return self.results["rho"]

    def get_slacks(self):
        """Return input, output and undesirable-output slacks."""
        _assert_optimized(self.results)
        slack_cols = [col for col in self.results.columns if col.startswith("slack_")]
        return self.results.loc[:, slack_cols]

    def get_lamda(self):
        """Return intensity variables with evaluated DMUs in rows."""
        _assert_optimized(self.lamda)
        return self.lamda

    def display_status(self):
        """Display the solver status for each evaluated DMU."""
        _assert_optimized(self.results)
        print(self.results["optimization_status"])

    def display_rho(self):
        """Display SBM efficiency scores."""
        _assert_optimized(self.results)
        print(self.results["rho"])

    def display_lamda(self):
        """Display intensity variables."""
        _assert_optimized(self.lamda)
        print(self.lamda)

    def info(self, dmu="all"):
        """Show Pyomo model information for selected DMU(s)."""
        if dmu == "all":
            for actual_index, model in self.__modeldict.items():
                print(f"\n--- SBM model for DMU: {actual_index} ---")
                model.pprint()
            return

        dmu_list = [dmu] if isinstance(dmu, (str, int, float)) else list(dmu)
        for actual_index in dmu_list:
            if actual_index not in self.__modeldict:
                print(f"DMU '{actual_index}' not found in the evaluated set.")
                continue
            print(f"\n--- SBM model for DMU: {actual_index} ---")
            self.__modeldict[actual_index].pprint()


def _direction_name(gx, gy, gb):
    has_x = sum(gx) >= 1
    has_y = sum(gy) >= 1
    has_b = sum(gb) >= 1

    if has_x and not has_y and not has_b:
        return "input_oriented"
    if has_y and not has_x and not has_b:
        return "output_oriented"
    if has_b and not has_x and not has_y:
        return "unoutput_oriented"
    if has_x and has_y and not has_b:
        return "hyper_orientedyx"
    if has_y and has_b and not has_x:
        return "hyper_orientedyb"
    if has_x and has_b and not has_y:
        return "hyper_orientedxb"
    if has_x and has_y and has_b:
        return "hyper_orientedyxb"
    raise ValueError(
        "gx, gy and gb must represent input_oriented, output_oriented, "
        "unoutput_oriented, hyper_orientedyx, hyper_orientedyb, "
        "hyper_orientedxb or hyper_orientedyxb direction."
    )


def _clean_efficiency(item):
    item = float(item)
    if abs(item) < 1e-10:
        return 0.0
    if 1.0 < item < 1.0 + 1e-8:
        return 1.0
    if -1e-8 < item < 0.0:
        return 0.0
    return item


def _assert_optimized(result):
    if result is None:
        raise RuntimeError("The model has not been optimized. Call optimize() first.")


# ---------------------------------------------------------------------------
# Slack-based directional distance function (SD).
#
# Färe & Grosskopf (2010, EJOR 200, 320-322) propose the SD measure: a
# slack-based directional distance function that maximizes the *sum of absolute
# slacks* with one independent slack variable per input and per output, instead
# of the SBM ratio of weighted-average relative slacks.  Because every variable
# carries its own non-negative slack, the model is translation-invariant (zeros
# and negatives are admissible in the reference data) and strictly more
# discriminating than the directional, hyperbolic and Farrell measures.  A DMU is
# SD-efficient iff the SD objective equals zero.
#
# Pham & Zelenyuk (2018, Empirical Economics 54, 153-187) extend SD to
# undesirable outputs under joint weak disposability of good and bad outputs, by
# adding an independent bad-output slack to the SD objective (§8.5.2).  The two
# public factories are :func:`sddf` (no undesirable outputs) and
# :func:`DEAweak.sddf` (undesirable outputs); both build this private engine.
# ---------------------------------------------------------------------------

class _SDDF:
    """Slack-based directional distance function (SD) efficiency model.

    Args:
        data (pandas.DataFrame): Input-output data.
        sent (str): Variable expression ``"inputs = outputs"`` or
            ``"inputs = outputs : undesirable_outputs"``.  Undesirable outputs
            are required by :func:`DEAweak.sddf` and rejected by :func:`sddf`.
        gx (list, optional): Input slack direction; a positive element marks the
            input as eligible for the input-excess slack.  Defaults to all ones.
        gy (list, optional): Desirable-output slack direction.  Defaults to all
            ones.
        gb (list, optional): Undesirable-output slack direction.  Defaults to all
            ones (only used when ``sent`` declares undesirable outputs).
        rts (str): ``RTS_CRS`` or ``RTS_VRS1``.
        evalindex, refindex (str, optional): Evaluated / reference sample filters.
    """

    def __init__(self, data, sent, gx=None, gy=None, gb=None,
                 rts=RTS_VRS1, evalindex=None, refindex=None):
        self.data = data
        self.sent = sent
        self.rts = rts
        self.evalindex = evalindex
        self.refindex = refindex

        if self.rts not in {RTS_CRS, RTS_VRS1}:
            raise ValueError("rts must be RTS_CRS or RTS_VRS1.")

        self.formula = parse_sent(sent)
        self.xcol = self.formula.x
        self.ycol = self.formula.y
        self.bcol = self.formula.b
        n_x, n_y, n_b = len(self.xcol), len(self.ycol), len(self.bcol)

        gb_default = [1] * n_b if n_b else []
        self.direction_spec = build_direction_spec(
            gx if gx is not None else [1] * n_x,
            gy if gy is not None else [1] * n_y,
            gb if gb is not None else gb_default,
            num_inputs=n_x,
            num_outputs=n_y,
            num_bad_outputs=n_b,
            allow_bad_direction=True,
        )
        self.gx = self.direction_spec.gx
        self.gy = self.direction_spec.gy
        self.gb = self.direction_spec.gb
        self.x_active = [j for j, v in enumerate(self.gx) if v > 0]
        self.y_active = [k for k, v in enumerate(self.gy) if v > 0]
        self.b_active = [h for h, v in enumerate(self.gb) if v > 0]

        prepared = prepare_dea_data(
            data, self.formula, eval_query=evalindex, ref_query=refindex,
            input_side="left", require_nonnegative_b=False,
        )
        self.evaluated_data_index = prepared.evaluated_index
        self.reference_data_index = prepared.reference_index
        self.x = prepared.x
        self.y = prepared.y
        self.xref = prepared.xref
        self.yref = prepared.yref
        self.b = prepared.b
        self.bref = prepared.bref
        if self.bcol and (np.any(self.b <= 0) or np.any(self.bref <= 0)):
            raise ValueError("Undesirable outputs must be strictly positive for SDDF.")

        self.evaluated_indices_range = range(len(self.evaluated_data_index))
        self.reference_indices_range = range(len(self.reference_data_index))
        self.num_inputs = n_x
        self.num_outputs = n_y
        self.num_unoutputs = n_b

        self.__modeldict = {}
        self.results = None
        self.lamda = None
        for i_range in self.evaluated_indices_range:
            actual_index = self.evaluated_data_index[i_range]
            self.__modeldict[actual_index] = self.__create_model(i_range, actual_index)

    def __create_model(self, current_dmu_range_index, actual_index):
        model = ConcreteModel()
        model.R = Set(initialize=self.reference_indices_range)
        model.J = Set(initialize=range(self.num_inputs))
        model.K = Set(initialize=range(self.num_outputs))
        model.L = Set(initialize=range(self.num_unoutputs))
        model.AX = Set(initialize=self.x_active)
        model.AY = Set(initialize=self.y_active)
        model.AB = Set(initialize=self.b_active)
        model.lamda = Var(model.R, bounds=(0.0, None), doc="intensity variables")
        model.ax = Var(model.AX, bounds=(0.0, None), doc="input excess slack")
        model.ay = Var(model.AY, bounds=(0.0, None), doc="output shortfall slack")
        model.ab = Var(model.AB, bounds=(0.0, None), doc="undesirable-output slack")

        model.objective = Objective(
            expr=(
                sum(model.ax[j] for j in model.AX)
                + sum(model.ay[k] for k in model.AY)
                + sum(model.ab[h] for h in model.AB)
            ),
            sense=maximize,
            doc="sum of absolute slacks (SD inefficiency)",
        )
        model.input = Constraint(model.J, rule=self.__input_rule(current_dmu_range_index))
        model.output = Constraint(model.K, rule=self.__output_rule(current_dmu_range_index))
        if self.num_unoutputs:
            model.unoutput = Constraint(model.L, rule=self.__unoutput_rule(current_dmu_range_index))
        if self.rts == RTS_VRS1:
            model.vrs = Constraint(rule=lambda m: sum(m.lamda[r] for r in m.R) == 1)
        return model

    def __input_rule(self, current_dmu_range_index):
        def input_rule(model, j):
            lhs = sum(model.lamda[r] * self.xref[r][j] for r in model.R)
            if j in self.x_active:
                return lhs + model.ax[j] == self.x[current_dmu_range_index][j]
            return lhs <= self.x[current_dmu_range_index][j]
        return input_rule

    def __output_rule(self, current_dmu_range_index):
        def output_rule(model, k):
            lhs = sum(model.lamda[r] * self.yref[r][k] for r in model.R)
            if k in self.y_active:
                return lhs - model.ay[k] == self.y[current_dmu_range_index][k]
            return lhs >= self.y[current_dmu_range_index][k]
        return output_rule

    def __unoutput_rule(self, current_dmu_range_index):
        def unoutput_rule(model, h):
            lhs = sum(model.lamda[r] * self.bref[r][h] for r in model.R)
            if h in self.b_active:
                return lhs + model.ab[h] == self.b[current_dmu_range_index][h]
            return lhs <= self.b[current_dmu_range_index][h]
        return unoutput_rule

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Solve each DMU's SDDF Pyomo model.

        Returns:
            pandas.DataFrame: Solver status, the SD inefficiency
            (``objective_value`` = sum of absolute slacks), a normalized
            efficiency ``te = 1 / (1 + objective_value)`` (1 for an efficient
            DMU, decreasing with the slack sum), and the per-variable slack
            columns ``slack_x_*`` / ``slack_y_*`` / ``slack_b_*``.
        """
        rows = []
        lamda_rows = []
        config = SolverConfig(solver=solver, email=email, tee=False)
        for actual_index, model in self.__modeldict.items():
            try:
                outcome = solve_pyomo_model(model, config=config, strict=False)
                optimization_status = outcome.legacy_status
            except Exception as exc:
                print(f"Error optimizing DMU {actual_index}: {exc}")
                optimization_status = "solve_error"

            if optimization_status in {"ok", "optimal"}:
                try:
                    row, lamda = self.__extract_solution(model, optimization_status)
                except Exception as exc:
                    print(f"Warning: could not retrieve SDDF results for DMU {actual_index}: {exc}")
                    row = self.__empty_result_row("result_error")
                    lamda = np.full(len(self.reference_data_index), np.nan)
            else:
                row = self.__empty_result_row(optimization_status)
                lamda = np.full(len(self.reference_data_index), np.nan)
            rows.append(pd.Series(row, name=actual_index))
            lamda_rows.append(pd.Series(lamda, index=self.reference_data_index, name=actual_index))

        self.results = pd.DataFrame(rows)
        self.lamda = pd.DataFrame(lamda_rows)
        return self.results

    def __extract_solution(self, model, optimization_status):
        ax = np.zeros(self.num_inputs)
        ay = np.zeros(self.num_outputs)
        ab = np.zeros(self.num_unoutputs)
        lamda = np.asarray([value(model.lamda[r]) for r in model.R], dtype=float)
        lamda = clean_small_values(lamda)
        for j in model.AX:
            ax[j] = value(model.ax[j])
        for k in model.AY:
            ay[k] = value(model.ay[k])
        for h in model.AB:
            ab[h] = value(model.ab[h])
        ax = clean_small_values(ax)
        ay = clean_small_values(ay)
        ab = clean_small_values(ab)
        objective_value = value(model.objective)
        obj_scalar = float(objective_value) if objective_value is not None else np.nan
        row = {
            "optimization_status": optimization_status,
            "objective_value": obj_scalar,
        }
        row["te"] = _clean_efficiency(
            1.0 / (1.0 + obj_scalar) if pd.notna(obj_scalar) and obj_scalar > -1 else np.nan
        )
        for j, name in enumerate(self.xcol):
            row[f"slack_x_{name}"] = ax[j]
        for k, name in enumerate(self.ycol):
            row[f"slack_y_{name}"] = ay[k]
        for h, name in enumerate(self.bcol):
            row[f"slack_b_{name}"] = ab[h]
        return row, lamda

    def __empty_result_row(self, status):
        row = {
            "optimization_status": status,
            "objective_value": np.nan,
            "te": np.nan,
        }
        for name in self.xcol:
            row[f"slack_x_{name}"] = np.nan
        for name in self.ycol:
            row[f"slack_y_{name}"] = np.nan
        for name in self.bcol:
            row[f"slack_b_{name}"] = np.nan
        return row

    def get_status(self):
        """Return the solver status for each evaluated DMU."""
        if self.results is None:
            raise RuntimeError("The model has not been optimized. Call optimize() first.")
        return self.results["optimization_status"]

    def get_objective(self):
        """Return the SD inefficiency (sum of absolute slacks) per DMU."""
        if self.results is None:
            raise RuntimeError("The model has not been optimized. Call optimize() first.")
        return self.results["objective_value"]

    def get_slacks(self):
        """Return input, output and undesirable-output slacks."""
        if self.results is None:
            raise RuntimeError("The model has not been optimized. Call optimize() first.")
        slack_cols = [col for col in self.results.columns if col.startswith("slack_")]
        return self.results.loc[:, slack_cols]

    def get_lamda(self):
        """Return intensity variables with evaluated DMUs in rows."""
        if self.lamda is None:
            raise RuntimeError("The model has not been optimized. Call optimize() first.")
        return self.lamda

    def display_status(self):
        """Display the solver status for each evaluated DMU."""
        print(self.get_status())

    def display_objective(self):
        """Display the SD inefficiency for each evaluated DMU."""
        print(self.get_objective())

    def display_lamda(self):
        """Display intensity variables."""
        print(self.get_lamda())

    def info(self, dmu="all"):
        """Show Pyomo model information for selected DMU(s)."""
        if dmu == "all":
            for actual_index, model in self.__modeldict.items():
                print(f"\n--- SDDF model for DMU: {actual_index} ---")
                model.pprint()
            return

        dmu_list = [dmu] if isinstance(dmu, (str, int, float)) else list(dmu)
        for actual_index in dmu_list:
            if actual_index not in self.__modeldict:
                print(f"DMU '{actual_index}' not found in the evaluated set.")
                continue
            print(f"\n--- SDDF model for DMU: {actual_index} ---")
            self.__modeldict[actual_index].pprint()
