"""Shared object-oriented bases for shadow-price dual DEA models."""

from __future__ import annotations

import ast
from typing import Any

import numpy as np
import pandas as pd
from pyomo.environ import ConcreteModel, Constraint, Objective, PositiveReals, Reals, Set, Var, maximize, minimize, value

from .constant import LEFT, RIGHT, OPT_DEFAULT, OPT_LOCAL, RTS_CRS, RTS_VRS
from .core import SolverConfig, solve_pyomo_model


class ShadowPriceDualBase:
    """One-model-per-DMU base for shadow-price dual models.

    Subclasses keep the economic formulas; this base centralises formula
    parsing, sample/reference selection, Pyomo set construction, solving,
    result-table formatting and ``info``.
    """

    objective_sense = minimize

    def __init__(self, data, sent, rts=RTS_VRS, baseindex=None, refindex=None):
        self.data = data
        self.sent = sent
        self.outputvars, self.inputvars, self.unoutputvars = self._parse_sent(sent)
        self.rts = rts
        self.baseindex = baseindex
        self.refindex = refindex
        self.y, self.x, self.b = self._select_data(baseindex)
        self.yref, self.xref, self.bref = self._select_data(refindex)
        self.xcol = self.x.columns
        self.ycol = self.y.columns
        self.bcol = self.b.columns
        self.I = self.x.index
        self._models = {}
        self.__modeldict = self._models
        for i in self.I:
            self.I0 = i
            self._models[i] = self._build_model()

    @staticmethod
    def _parse_sent(sent):
        left, right = sent.split("=", 1)
        output_part, bad_part = right.split(":", 1)
        inputvars = [item for item in left.strip().split() if item]
        outputvars = [item for item in output_part.strip().split() if item]
        unoutputvars = [item for item in bad_part.strip().split() if item]
        return outputvars, inputvars, unoutputvars

    def _select_data(self, selector):
        if selector is not None:
            column, raw_values = selector.split("=", 1)
            values = ast.literal_eval(raw_values)
            if not isinstance(values, (list, tuple, set, pd.Index, np.ndarray)):
                values = [values]
            frame = self.data.loc[self.data[column.strip()].isin(values)]
        else:
            frame = self.data
        return (
            frame.loc[:, self.outputvars],
            frame.loc[:, self.inputvars],
            frame.loc[:, self.unoutputvars],
        )

    def _build_model(self):
        model = ConcreteModel()
        model.I2 = Set(initialize=self.xref.index)
        model.K = Set(initialize=range(len(self.x.iloc[0])))
        model.L = Set(initialize=range(len(self.y.iloc[0])))
        model.J = Set(initialize=range(len(self.b.iloc[0])))
        self._add_shadow_price_variables(model)
        model.objective = Objective(rule=self._objective_rule(), sense=self._objective_sense())
        model.first = Constraint(model.I2, rule=self._technology_rule(), doc="first constraint")
        model.second = Constraint(model.I2, rule=self._input_shadow_rule(), doc="second constraint")
        self._add_model_specific_constraints(model)
        return model

    def _add_shadow_price_variables(self, model):
        model.spx = Var(model.K, initialize=1, bounds=(0.0, None), within=Reals, doc="shadow price of x")
        model.spy = Var(model.L, initialize=1, bounds=(0.0, None), within=Reals, doc="shadow price of y")
        model.spb = Var(model.J, bounds=(1e-6, None), within=PositiveReals, doc="shadow price of b")
        if self.rts == RTS_VRS:
            model.spalpha = Var(Set(initialize=range(1)), within=Reals, doc="shadow price of 1")

    def _alpha_term(self, model):
        return model.spalpha[0] if self.rts == RTS_VRS else 0

    def _technology_expression(self, model, row):
        return (
            sum(model.spx[k] * self.xref.loc[row, self.xcol[k]] for k in model.K)
            - sum(model.spy[l] * self.yref.loc[row, self.ycol[l]] for l in model.L)
            + sum(model.spb[j] * self.bref.loc[row, self.bcol[j]] for j in model.J)
            + self._alpha_term(model)
        )

    def _input_shadow_expression(self, model, row):
        return sum(model.spx[k] * self.xref.loc[row, self.xcol[k]] for k in model.K) + self._alpha_term(model)

    def _technology_rule(self):
        def first_rule(model, i2):
            return self._technology_expression(model, i2) >= 0
        return first_rule

    def _input_shadow_rule(self):
        def second_rule(model, i2):
            return self._input_shadow_expression(model, i2) >= 0
        return second_rule

    def _objective_sense(self):
        return self.objective_sense

    def _objective_rule(self):
        raise NotImplementedError

    def _add_model_specific_constraints(self, model):
        raise NotImplementedError

    def optimize(self, solver=OPT_DEFAULT):
        obj, spx, spy, spb = {}, {}, {}, {}
        config = SolverConfig(solver=solver, email=OPT_LOCAL, tee=False)
        for ind, problem in self._models.items():
            outcome = solve_pyomo_model(problem, config=config, strict=False)
            if outcome.ok:
                obj[ind] = value(problem.objective, exception=False)
                spx[ind] = np.asarray([value(problem.spx[k], exception=False) for k in problem.K])
                spy[ind] = np.asarray([value(problem.spy[l], exception=False) for l in problem.L])
                spb[ind] = np.asarray([value(problem.spb[j], exception=False) for j in problem.J])
            else:
                obj[ind] = np.nan
                spx[ind] = np.full(len(list(problem.K)), np.nan)
                spy[ind] = np.full(len(list(problem.L)), np.nan)
                spb[ind] = np.full(len(list(problem.J)), np.nan)
        obj_df = pd.DataFrame(obj, index=["obj"]).T
        spx_df = pd.DataFrame(spx).T
        spx_df.columns = spx_df.columns.map(lambda x: "Input" + str(x) + "'s shadow price")
        spy_df = pd.DataFrame(spy).T
        spy_df.columns = spy_df.columns.map(lambda y: "Output" + str(y) + "'s shadow price")
        spb_df = pd.DataFrame(spb).T
        spb_df.columns = spb_df.columns.map(lambda b: "Undesirable Output" + str(b) + "'s shadow price")
        return pd.concat([obj_df, spx_df, spy_df, spb_df], axis=1)

    def info(self, dmu="all"):
        if dmu == "all":
            for ind, problem in self._models.items():
                print(ind, "\n", problem.pprint())
            return None
        key = dmu
        if key not in self._models:
            try:
                key = int(dmu)
            except (TypeError, ValueError):
                pass
        if key in self._models:
            print(self._models[key].pprint())
        else:
            print(f"DMU '{dmu}' not found in the evaluated set.")
        return None


class DirectionalDistanceDualBase(ShadowPriceDualBase):
    """Dual model for the directional distance function."""

    def __init__(self, data, sent, gy, gx, gb, rts=RTS_VRS, baseindex=None, refindex=None):
        self._raw_gy, self._raw_gx, self._raw_gb = gy, gx, gb
        self.gy, self.gx, self.gb = gy, gx, gb
        super().__init__(data, sent, rts=rts, baseindex=baseindex, refindex=refindex)
        self.gy = self._direction_frame(self._raw_gy, self.ycol, self.y.index)
        self.gx = self._direction_frame(self._raw_gx, self.xcol, self.x.index)
        self.gb = self._direction_frame(self._raw_gb, self.bcol, self.b.index)

    @staticmethod
    def _direction_frame(values, columns, index):
        if isinstance(values, pd.DataFrame):
            return values.loc[:, columns]
        if isinstance(values, pd.Series):
            values = values.to_list()
        array = np.asarray(values, dtype=float).reshape(1, -1)
        if array.shape[1] != len(columns):
            raise ValueError("Direction vector length must match the number of variables.")
        return pd.DataFrame(np.repeat(array, len(index), axis=0), index=index, columns=columns)

    def _direction_at(self, direction, columns, row, idx):
        if isinstance(direction, pd.DataFrame):
            return direction.loc[row, columns[idx]]
        if isinstance(direction, pd.Series):
            if columns[idx] in direction.index:
                return direction.loc[columns[idx]]
            return direction.iloc[idx]
        return list(direction)[idx]

    def _objective_rule(self):
        def objective_rule(model):
            return (
                sum(model.spx[k] * self.x.loc[self.I0, self.xcol[k]] for k in model.K)
                - sum(model.spy[l] * self.y.loc[self.I0, self.ycol[l]] for l in model.L)
                + sum(model.spb[j] * self.b.loc[self.I0, self.bcol[j]] for j in model.J)
                + self._alpha_term(model)
            )
        return objective_rule

    def _add_model_specific_constraints(self, model):
        model.third = Constraint(rule=self._normalization_rule(), doc="third constraint")

    def _normalization_rule(self):
        def third_rule(model):
            return (
                sum(model.spx[k] * self._direction_at(self.gx, self.xcol, self.I0, k) for k in model.K)
                + sum(model.spy[l] * self._direction_at(self.gy, self.ycol, self.I0, l) for l in model.L)
                + sum(model.spb[j] * self._direction_at(self.gb, self.bcol, self.I0, j) for j in model.J)
                == 1
            )
        return third_rule


class DirectionalResponseDualBase(ShadowPriceDualBase):
    """Dual model for the directional response function."""

    def __init__(self, data, year, sent, fenmu, fenzi, side=LEFT, rts=RTS_VRS, baseindex=None, refindex=None):
        self.year = year
        self.tlt = pd.Series(year).drop_duplicates().sort_values().reset_index(drop=True)
        self.fenmu = fenmu
        self.fenzi = fenzi
        self.side = side
        self.outputvars, self.inputvars, self.unoutputvars = self._parse_sent(sent)
        self.obj_coeflt, self.rule4_coeflt, self.neg_obj = self._build_ratio_coefficients(fenmu, fenzi)
        self.xobj_coef = self.obj_coeflt["xobj_coef"]
        self.yobj_coef = self.obj_coeflt["yobj_coef"]
        self.bobj_coef = self.obj_coeflt["bobj_coef"]
        self.xrule4_coef = self.rule4_coeflt["xrule4_coef"]
        self.yrule4_coef = self.rule4_coeflt["yrule4_coef"]
        self.brule4_coef = self.rule4_coeflt["brule4_coef"]
        super().__init__(data, sent, rts=rts, baseindex=baseindex, refindex=refindex)

    def _build_ratio_coefficients(self, fenmu, fenzi):
        variables = self.inputvars + self.outputvars + self.unoutputvars
        if fenmu not in variables:
            raise ValueError("fenmu must be in sent.")
        if fenzi not in variables:
            raise ValueError("fenzi must be in sent.")
        obj = {
            "xobj_coef": [1 if fenzi == item else 0 for item in self.inputvars],
            "yobj_coef": [1 if fenzi == item else 0 for item in self.outputvars],
            "bobj_coef": [1 if fenzi == item else 0 for item in self.unoutputvars],
        }
        rule4 = {
            "xrule4_coef": [1 if fenmu == item else 0 for item in self.inputvars],
            "yrule4_coef": [1 if fenmu == item else 0 for item in self.outputvars],
            "brule4_coef": [1 if fenmu == item else 0 for item in self.unoutputvars],
        }
        neg_obj = fenmu in self.inputvars or fenmu in self.unoutputvars
        return obj, rule4, neg_obj

    def _add_shadow_price_variables(self, model):
        model.spx = Var(model.K, bounds=self._bounds_from_coefficients(self.xobj_coef), within=Reals, doc="shadow price of x")
        model.spy = Var(model.L, bounds=self._bounds_from_coefficients(self.yobj_coef), within=Reals, doc="shadow price of y")
        model.spb = Var(model.J, bounds=self._bounds_from_coefficients(self.bobj_coef), within=Reals, doc="shadow price of b")
        if self.rts == RTS_VRS:
            model.spalpha = Var(Set(initialize=range(1)), within=Reals, doc="shadow price of 1")

    @staticmethod
    def _bounds_from_coefficients(coefficients):
        lower = {idx: None if coef > 0 else 0 for idx, coef in enumerate(coefficients)}
        upper = {idx: None for idx, _ in enumerate(coefficients)}

        def bounds(_model, idx):
            return lower[idx], upper[idx]
        return bounds

    def _objective_sense(self):
        if self.side == RIGHT:
            return maximize if self.neg_obj else minimize
        if self.side == LEFT:
            return minimize if self.neg_obj else maximize
        raise ValueError("side must be LEFT or RIGHT.")

    def _objective_rule(self):
        sign = -1 if self.neg_obj else 1

        def objective_rule(model):
            return sign * (
                sum(model.spx[k] * self.xobj_coef[k] for k in model.K)
                - sum(model.spy[l] * self.yobj_coef[l] for l in model.L)
                + sum(model.spb[j] * self.bobj_coef[j] for j in model.J)
            )
        return objective_rule

    def _add_model_specific_constraints(self, model):
        model.third = Constraint(rule=self._zero_inefficiency_rule(), doc="inefficiency=0 constraint")
        model.forth = Constraint(rule=self._denominator_rule(), doc="fenmu=1 constraint")

    def _zero_inefficiency_rule(self):
        def third_rule(model):
            return (
                sum(model.spx[k] * self.x.loc[self.I0, self.xcol[k]] for k in model.K)
                - sum(model.spy[l] * self.y.loc[self.I0, self.ycol[l]] for l in model.L)
                + sum(model.spb[j] * self.b.loc[self.I0, self.bcol[j]] for j in model.J)
                + self._alpha_term(model)
                == 0
            )
        return third_rule

    def _denominator_rule(self):
        def forth_rule(model):
            return (
                sum(self.xrule4_coef[k] * model.spx[k] for k in model.K)
                + sum(self.yrule4_coef[l] * model.spy[l] for l in model.L)
                + sum(self.brule4_coef[j] * model.spb[j] for j in model.J)
                == 1
            )
        return forth_rule
