"""Slacks-based measure (SBM) DEA models implemented with Pyomo.

The SBM class follows the direction-vector convention used by ``DDFweak2``:
``gx`` marks input directions, ``gy`` marks desirable-output directions, and
``gb`` marks undesirable-output directions.  The active direction name is stored
in ``self.direction`` as one of ``input_oriented``, ``output_oriented``,
``unoutput_oriented``, ``hyper_orientedyx``, ``hyper_orientedyb``,
``hyper_orientedxb`` or ``hyper_orientedyxb``.
"""

import numpy as np
import pandas as pd
from pyomo.environ import ConcreteModel, Constraint, Objective, Set, Var, maximize, minimize, value

from .constant import OPT_DEFAULT, OPT_LOCAL, RTS_CRS, RTS_VRS1
from .core import (
    SolverConfig,
    build_direction_spec,
    clean_small_values,
    parse_sent,
    prepare_dea_data,
    solve_pyomo_model,
)


class SBM:
    """Slacks-based measure DEA model.

    Args:
        data (pandas.DataFrame): Input-output data.
        sent (str): Variable expression using the package convention
            ``"inputs = outputs"`` or ``"inputs = outputs : undesirable_outputs"``.
            Examples: ``"K L = Y"`` and ``"K L E = Y : CO2"``.
        gy (list): Desirable-output direction vector. Same convention as
            ``DDFweak2``. A positive element means that the corresponding
            desirable output is adjusted upward.
        gx (list): Input direction vector. A positive element means that the
            corresponding input is adjusted downward.
        gb (list): Undesirable-output direction vector. A positive element
            means that the corresponding undesirable output is adjusted downward.
        rts (str): ``RTS_CRS`` or ``RTS_VRS1``.
        baseindex (str, optional): Evaluated-sample filter such as
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
        baseindex=None,
        refindex=None,
    ):
        self.data = data
        self.sent = sent
        self.rts = rts
        self.baseindex = baseindex
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
            eval_query=baseindex,
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
