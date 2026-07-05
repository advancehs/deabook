"""Object-oriented backend for the basic material-balance DEA model.

The public :class:`deabook.MB.MB` entry point historically dispatches to four
internal classes named by whether non-polluting and polluting desirable outputs
are present: ``MB1111``, ``MB1110``, ``MB1101`` and ``MB1100``.  The old file
implemented those four cases by copy-pasting almost the same Pyomo model.  This
module keeps the old class names, constructor signatures and result columns, but
moves the shared model lifecycle into ``_BasicMBBase``.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from pyomo.environ import ConcreteModel, Constraint, Objective, Reals, Set, Var, minimize, value

from ..constant import OPT_DEFAULT, RTS_VRS1
from ..core import SolverConfig, select_rows, solve_pyomo_model


class _BasicMBBase:
    """Shared implementation for ``MB1111``/``MB1110``/``MB1101``/``MB1100``.

    Subclasses only normalize their historical positional arguments into the
    four optional variable groups.  The base class preserves the old Pyomo
    component names (``Knp``, ``Kp``, ``Lp``, ``Lnp``, ``objb`` and so on), the
    ``level``-controlled variable targets, and the compact result table returned
    by ``optimize``.
    """

    def _init_basic_mb(
        self,
        data,
        inputvars_np,
        inputvars_p,
        outputvars_np,
        outputvars_p,
        unoutputvars,
        sx,
        sy,
        rts,
        level,
        baseindex,
        refindex,
    ):
        self.data = data
        self.inputvars_np = list(inputvars_np)
        self.inputvars_p = list(inputvars_p)
        self.outputvars_np = None if outputvars_np is None else list(outputvars_np)
        self.outputvars_p = None if outputvars_p is None else list(outputvars_p)
        self.unoutputvars = list(unoutputvars)
        self.sx = sx
        self.sy = sy
        self.rts = rts
        self.level = level
        self.baseindex = baseindex
        self.refindex = refindex

        self.sx_np = np.asarray(self.sx)[:, : len(self.inputvars_np)]
        self.sx_p = np.asarray(self.sx)[:, len(self.inputvars_np) :]
        if self.outputvars_np is None:
            self.sy_np = None
            self.sy_p = np.asarray(self.sy) if self.outputvars_p is not None else None
        else:
            split_at = len(self.outputvars_np)
            self.sy_np = np.asarray(self.sy)[:, :split_at]
            self.sy_p = np.asarray(self.sy)[:, split_at:] if self.outputvars_p is not None else None

        eval_df = select_rows(self.data, self.baseindex, "baseindex")
        ref_df = select_rows(self.data, self.refindex, "refindex")
        self._assign_frames(eval_df, ref_df)
        self._build_all_models()

    @property
    def _has_y_np(self):
        return self.outputvars_np is not None

    @property
    def _has_y_p(self):
        return self.outputvars_p is not None

    @property
    def _output_np_variable_level(self):
        return 5 if self._has_y_p else 4

    def _assign_frames(self, eval_df, ref_df):
        self.x_p = eval_df.loc[:, self.inputvars_p]
        self.x_np = eval_df.loc[:, self.inputvars_np]
        self.b = eval_df.loc[:, self.unoutputvars]
        self.xref_p = ref_df.loc[:, self.inputvars_p]
        self.xref_np = ref_df.loc[:, self.inputvars_np]
        self.bref = ref_df.loc[:, self.unoutputvars]

        if self._has_y_np:
            self.y_np = eval_df.loc[:, self.outputvars_np]
            self.yref_np = ref_df.loc[:, self.outputvars_np]
            self.ycol_np = self.y_np.columns
        else:
            self.y_np = None
            self.yref_np = None
            self.ycol_np = None

        if self._has_y_p:
            self.y_p = eval_df.loc[:, self.outputvars_p]
            self.yref_p = ref_df.loc[:, self.outputvars_p]
            self.ycol_p = self.y_p.columns
        else:
            self.y_p = None
            self.yref_p = None
            self.ycol_p = None

        self.xcol_p = self.x_p.columns
        self.xcol_np = self.x_np.columns
        self.bcol = self.b.columns
        self.I = self.x_p.index

    def _build_all_models(self):
        self.__modeldict = {}
        self._models = self.__modeldict
        for i in self.I:
            self.I0 = i
            self.__modeldict[i] = self._build_model()

    def _build_model(self):
        model = ConcreteModel()
        model.I2 = Set(initialize=self.xref_p.index)
        model.Knp = Set(initialize=range(len(self.xcol_np)))
        model.Kp = Set(initialize=range(len(self.xcol_p)))
        if self._has_y_np:
            model.Lnp = Set(initialize=range(len(self.ycol_np)))
        if self._has_y_p:
            model.Lp = Set(initialize=range(len(self.ycol_p)))
        model.B = Set(initialize=range(len(self.bcol)))

        model.thetax_p = Var(model.Kp, bounds=(0.0, None), doc="slack x_p")
        model.thetax_np = Var(model.Knp, bounds=(0.0, None), doc="slack x_np")
        if self._has_y_p:
            model.thetay_p = Var(model.Lp, bounds=(0.0, None), doc="slack y_p")
        if self._has_y_np:
            model.thetay_np = Var(model.Lnp, bounds=(0.0, None), doc="slack y_np")
        model.thetab = Var(model.B, bounds=(0.0, None), doc="slack b")
        model.objb = Var(model.B, bounds=(0.0, None), within=Reals, doc="object b")

        if self.level >= 2:
            model.objx_p = Var(model.Kp, bounds=(0.0, None), within=Reals, doc="object x_p")
        if self.level >= 3:
            model.objx_np = Var(model.Knp, bounds=(0.0, None), within=Reals, doc="object x_np")
        if self._has_y_p and self.level >= 4:
            model.objy_p = Var(model.Lp, bounds=(0.0, None), within=Reals, doc="object y_p")
        if self._has_y_np and self.level >= self._output_np_variable_level:
            model.objy_np = Var(model.Lnp, bounds=(0.0, None), within=Reals, doc="object y_np")

        model.lamda = Var(model.I2, bounds=(0.0, None), within=Reals, doc="intensity variables")
        model.objective = Objective(rule=self._objective_rule(), sense=minimize, doc="objective function")
        model.input_np = Constraint(model.Knp, rule=self._input_np_rule(), doc="input_np constraint")
        model.input_p = Constraint(model.Kp, rule=self._input_p_rule(), doc="input_p constraint")
        if self._has_y_np:
            model.output_np = Constraint(model.Lnp, rule=self._output_np_rule(), doc="output_np constraint")
        if self._has_y_p:
            model.output_p = Constraint(model.Lp, rule=self._output_p_rule(), doc="output_p constraint")
        model.undesirable_output = Constraint(model.B, rule=self._undesirable_output_rule(), doc="undesirable output constraint")
        model.mb = Constraint(model.B, rule=self._mb_rule(), doc="material balance constraint")
        if self.rts == RTS_VRS1:
            model.vrs = Constraint(rule=self._vrs_rule(), doc="various return to scale rule")
        return model

    def _objective_rule(self):
        def objective_rule(model):
            return sum(model.objb[b] for b in model.B)

        return objective_rule

    def _input_p_rhs(self, model, kp):
        if self.level >= 2:
            return model.objx_p[kp]
        return self.x_p.loc[self.I0, self.xcol_p[kp]]

    def _input_np_rhs(self, model, knp):
        if self.level >= 3:
            return model.objx_np[knp]
        return self.x_np.loc[self.I0, self.xcol_np[knp]]

    def _output_p_rhs(self, model, lp):
        if self.level >= 4:
            return model.objy_p[lp]
        return self.y_p.loc[self.I0, self.ycol_p[lp]]

    def _output_np_rhs(self, model, lnp):
        if self.level >= self._output_np_variable_level:
            return model.objy_np[lnp]
        return self.y_np.loc[self.I0, self.ycol_np[lnp]]

    def _input_p_rule(self):
        def input_p_rule(model, kp):
            return sum(model.lamda[i2] * self.xref_p.loc[i2, self.xcol_p[kp]] for i2 in model.I2) + model.thetax_p[kp] == self._input_p_rhs(model, kp)

        return input_p_rule

    def _input_np_rule(self):
        def input_np_rule(model, knp):
            return sum(model.lamda[i2] * self.xref_np.loc[i2, self.xcol_np[knp]] for i2 in model.I2) + model.thetax_np[knp] == self._input_np_rhs(model, knp)

        return input_np_rule

    def _output_p_rule(self):
        def output_p_rule(model, lp):
            return sum(model.lamda[i2] * self.yref_p.loc[i2, self.ycol_p[lp]] for i2 in model.I2) - model.thetay_p[lp] == self._output_p_rhs(model, lp)

        return output_p_rule

    def _output_np_rule(self):
        def output_np_rule(model, lnp):
            return sum(model.lamda[i2] * self.yref_np.loc[i2, self.ycol_np[lnp]] for i2 in model.I2) - model.thetay_np[lnp] == self._output_np_rhs(model, lnp)

        return output_np_rule

    def _undesirable_output_rule(self):
        def undesirable_output_rule(model, b):
            return sum(model.lamda[i2] * self.bref.loc[i2, self.bcol[b]] for i2 in model.I2) + model.thetab[b] == model.objb[b]

        return undesirable_output_rule

    def _input_pollution_gap(self, model, kp):
        if self.level >= 2:
            return self.x_p.loc[self.I0, self.xcol_p[kp]] - model.objx_p[kp]
        return model.thetax_p[kp]

    def _output_pollution_gap(self, model, lp):
        if self.level >= 4:
            return model.objy_p[lp] - self.y_p.loc[self.I0, self.ycol_p[lp]]
        return model.thetay_p[lp]

    def _mb_lhs(self, model, b):
        input_term = sum(self.sx_p[b][kp] * self._input_pollution_gap(model, kp) for kp in model.Kp)
        output_term = 0 if not self._has_y_p else sum(self.sy_p[b][lp] * self._output_pollution_gap(model, lp) for lp in model.Lp)
        return input_term + output_term

    def _mb_rule(self):
        def mb_rule(model, b):
            return self._mb_lhs(model, b) == self.b.loc[self.I0, self.bcol[b]] - model.objb[b]

        return mb_rule

    def _vrs_rule(self):
        def vrs_rule(model):
            return sum(model.lamda[i2] for i2 in model.I2) == 1

        return vrs_rule

    def optimize(self, solver=OPT_DEFAULT):
        data2 = pd.DataFrame()
        obj = {}
        objb = {}
        config = SolverConfig(solver=solver)
        for ind, problem in self.__modeldict.items():
            outcome = solve_pyomo_model(problem, config=config, strict=False)
            data2.loc[ind, "optimization_status"] = outcome.legacy_status
            if outcome.ok:
                obj[ind] = value(problem.objective)
                objb[ind] = [value(problem.objb[b]) for b in problem.B]
            else:
                obj[ind] = np.nan
                objb[ind] = [np.nan for _ in problem.B]

        obj_df = pd.DataFrame(obj, index=["obj"]).T
        objb_df = pd.DataFrame(objb).T
        if len(objb_df.columns) == 1:
            objb_df.columns = ["best of Undesirable"]
        else:
            objb_df.columns = [f"best of Undesirable{b}" for b in objb_df.columns]
        return pd.concat([data2, obj_df, objb_df], axis=1)

    def info(self, dmu="all"):
        """Keep the historical quiet ``info`` behaviour for these internals."""
        return None


class MB1111(_BasicMBBase):
    """Basic MB model with ``x_np + x_p = y_np + y_p : b``."""

    def __init__(self, data, inputvars_np, inputvars_p, outputvars_np, outputvars_p, unoutputvars, sx, sy, rts, level, baseindex, refindex):
        self._init_basic_mb(data, inputvars_np, inputvars_p, outputvars_np, outputvars_p, unoutputvars, sx, sy, rts, level, baseindex, refindex)


class MB1110(_BasicMBBase):
    """Basic MB model with ``x_np + x_p = y_np : b``."""

    def __init__(self, data, inputvars_np, inputvars_p, outputvars_np, unoutputvars, sx, sy, rts, level, baseindex, refindex):
        self._init_basic_mb(data, inputvars_np, inputvars_p, outputvars_np, None, unoutputvars, sx, sy, rts, level, baseindex, refindex)


class MB1101(_BasicMBBase):
    """Basic MB model with ``x_np + x_p = y_p : b``."""

    def __init__(self, data, inputvars_np, inputvars_p, outputvars_p, unoutputvars, sx, sy, rts, level, baseindex, refindex):
        self._init_basic_mb(data, inputvars_np, inputvars_p, None, outputvars_p, unoutputvars, sx, sy, rts, level, baseindex, refindex)


class MB1100(_BasicMBBase):
    """Basic MB model with ``x_np + x_p = : b``."""

    def __init__(self, data, inputvars_np, inputvars_p, unoutputvars, sx, sy, rts, level, baseindex, refindex):
        self._init_basic_mb(data, inputvars_np, inputvars_p, None, None, unoutputvars, sx, sy, rts, level, baseindex, refindex)


__all__ = ["MB1111", "MB1110", "MB1101", "MB1100"]
