"""Main module."""
# import dependencies

from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, maximize, Constraint, Reals
import numpy as np
import pandas as pd
from .constant import RTS_VRS1, OPT_DEFAULT
from .utils import tools,MB_
from .core import SolverConfig, solve_pyomo_model
import ast


class MB():
    def __init__(self, data,sent = "inputvar_np + inputvar_p =outputvar_np + outputvar_p:unoutputvar",  \
                 sx=[[1,1,1],[1,1,1]], sy=[[1],[1]], level=5 ,rts=RTS_VRS1, baseindex=None,refindex=None):
        """DEA: Directional distance function

        Args:
            data (pandas.DataFrame): input pandas.
            sent (str): inputvars=outputvars[: unoutputvars]. e.g.: "K L+ E = Y : CO2" 使用“+”区分污染和废污染投入（或产出）
            sx (list): 投入包含污染物质系数. Defaults to [[1,1,1],[1,1,1]].
            sy (list, optional): 期望产出包含污染物质系数. Defaults to [[1],[1]].
            level(int, optional): 返回求的层级。1：只变量化b；2：再变量化x_p；3：再变量化x_np；\
                                                4：再变量化y_p（没有则变量化y_np）；5：再变量化y_np
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # Initialize DEA model
        self.data=data
        # self.year = year
        self.sent = sent
        # self.tlt=pd.Series(self.year).drop_duplicates().sort_values()
        self.rts = rts
        self.baseindex = baseindex
        self.refindex = refindex

        self.inputvars_np, self.inputvars_p,self.outputvars_np,\
            self.outputvars_p,self.unoutputvars, self.sx, self.sy,self.level = tools.split_MB(self.sent, sx, sy,level)

    def _backend_spec(self):
        """Return the legacy MB backend class and its variable-group arguments."""
        shared_args = [self.unoutputvars, self.sx, self.sy, self.rts,
                       self.level, self.baseindex, self.refindex]
        if self.outputvars_np is not None and self.outputvars_p is not None:
            return MB_.MB1111, [
                self.inputvars_np, self.inputvars_p,
                self.outputvars_np, self.outputvars_p, *shared_args
            ]
        if self.outputvars_np is not None:
            return MB_.MB1110, [
                self.inputvars_np, self.inputvars_p,
                self.outputvars_np, *shared_args
            ]
        if self.outputvars_p is not None:
            return MB_.MB1101, [
                self.inputvars_np, self.inputvars_p,
                self.outputvars_p, *shared_args
            ]
        return MB_.MB1100, [self.inputvars_np, self.inputvars_p, *shared_args]

    def _make_backend(self):
        backend_class, args = self._backend_spec()
        return backend_class(self.data, *args)

    def optimize(self, solver=OPT_DEFAULT, dmu="1"):
        backend = self._make_backend()
        return backend.optimize(solver=solver), backend.info(dmu=dmu)

    def info(self, dmu="all"):
        return self._make_backend().info(dmu=dmu)



class _PanelMBBase(MB):
    """Shared object-oriented base for panel material-balance models."""

    _has_objx = False
    _has_objy = False
    _has_theta = False
    _objective_sense = minimize
    _objective_component = "thetab"
    _include_var_undesirable = False
    _print_all_info = True

    def __init__(self, data, year, sent="inputvar=outputvar", sx=[[1, 1, 1], [1, 1, 1]],
                 sy=[[1], [1]], rts=RTS_VRS1, baseindex=None, refindex=None):
        self.data = data
        self.year = year
        self.sent = sent
        self.tlt = pd.Series(self.year).drop_duplicates().sort_values().reset_index(drop=True)
        self.inputvars, self.outputvars, self.unoutputvars = self._parse_sent(sent)
        self.sx, self.sy = sx, sy
        self.rts = rts
        self.baseindex = baseindex
        self.refindex = refindex
        self.y, self.x, self.b = self._select_data(baseindex)
        self.yref, self.xref, self.bref = self._select_data(refindex)
        self.xcol = self.x.columns
        self.ycol = self.y.columns
        self.bcol = self.b.columns if self.b is not None else None
        self.I = self.x.index
        self._models = {}
        self.__modeldict = self._models
        for i in self.I:
            self.I0 = i
            self._models[i] = self._build_model()

    @staticmethod
    def _parse_sent(sent):
        left, right = sent.split('=', 1)
        inputvars = [item for item in left.strip().split() if item]
        if ':' in right:
            output_part, unoutput_part = right.split(':', 1)
            outputvars = [item for item in output_part.strip().split() if item]
            unoutputvars = [item for item in unoutput_part.strip().split() if item]
        else:
            outputvars = [item for item in right.strip().split() if item]
            unoutputvars = None
        return inputvars, outputvars, unoutputvars or None

    def _select_data(self, selector):
        if selector is not None:
            varname, raw_values = selector.split('=', 1)
            values = ast.literal_eval(raw_values)
            if not isinstance(values, (list, tuple, set, pd.Index, np.ndarray)):
                values = [values]
            frame = self.data.loc[self.data[varname.strip()].isin(values)]
        else:
            frame = self.data
        y = frame.loc[:, self.outputvars]
        x = frame.loc[:, self.inputvars]
        b = frame.loc[:, self.unoutputvars] if self.unoutputvars is not None else None
        return y, x, b

    def _build_model(self):
        model = ConcreteModel()
        model.I2 = Set(initialize=self.xref.index)
        model.K = Set(initialize=range(len(self.x.iloc[0])))
        model.L = Set(initialize=range(len(self.y.iloc[0])))
        if self.b is not None:
            model.B = Set(initialize=range(len(self.b.iloc[0])))
        self._add_variables(model)
        model.objective = Objective(rule=self._objective_rule(), sense=self._objective_sense,
                                    doc='objective function')
        model.input = Constraint(model.K, rule=self._input_rule(), doc='input constraint')
        model.output = Constraint(model.L, rule=self._output_rule(), doc='output constraint')
        if self.b is not None:
            model.undesirable_output = Constraint(
                model.B, rule=self._undesirable_output_rule(), doc='undesirable output constraint'
            )
            model.mb = Constraint(model.B, rule=self._mb_rule(), doc='material balance constraint')
        if self.rts == RTS_VRS1:
            model.vrs = Constraint(rule=self._vrs_rule(), doc='various return to scale rule')
        return model

    def _add_variables(self, model):
        if self._has_objx:
            model.objx = Var(model.K, bounds=(0.0, None), within=Reals, doc='object x')
        if self._has_objy:
            model.objy = Var(model.L, bounds=(0.0, None), within=Reals, doc='object y')
        model.thetax = Var(model.K, bounds=(0.0, None), doc='slack x')
        model.thetay = Var(model.L, bounds=(0.0, None), doc='slack y')
        if self.b is not None:
            model.thetab = Var(model.B, bounds=(0.0, None), doc='slack b')
            if self._has_theta:
                model.theta = Var(model.B, bounds=(0.0, None), within=Reals, doc='object b')
        model.lamda = Var(model.I2, bounds=(0.0, None), within=Reals, doc='intensity variables')

    def _objective_rule(self):
        def objective_rule(model):
            component = getattr(model, self._objective_component)
            return sum(component[b] for b in model.B)
        return objective_rule

    def _input_rhs(self, model, k):
        return model.objx[k] if self._has_objx else self.x.loc[self.I0, self.xcol[k]]

    def _output_rhs(self, model, l):
        return model.objy[l] if self._has_objy else self.y.loc[self.I0, self.ycol[l]]

    def _undesirable_rhs(self, model, b):
        return model.theta[b] if self._has_theta else self.b.loc[self.I0, self.bcol[b]]

    def _input_rule(self):
        def input_rule(model, k):
            return sum(model.lamda[i2] * self.xref.loc[i2, self.xcol[k]] for i2 in model.I2) \
                + model.thetax[k] == self._input_rhs(model, k)
        return input_rule

    def _output_rule(self):
        def output_rule(model, l):
            return sum(model.lamda[i2] * self.yref.loc[i2, self.ycol[l]] for i2 in model.I2) \
                - model.thetay[l] == self._output_rhs(model, l)
        return output_rule

    def _undesirable_output_rule(self):
        def undesirable_output_rule(model, b):
            return sum(model.lamda[i2] * self.bref.loc[i2, self.bcol[b]] for i2 in model.I2) \
                + model.thetab[b] == self._undesirable_rhs(model, b)
        return undesirable_output_rule

    def _x_gap(self, model, k):
        if self._has_objx:
            return model.thetax[k] + self.x.loc[self.I0, self.xcol[k]] - model.objx[k]
        return model.thetax[k]

    def _y_gap(self, model, l):
        if self._has_objy:
            return model.thetay[l] + model.objy[l] - self.y.loc[self.I0, self.ycol[l]]
        return model.thetay[l]

    def _mb_rhs(self, model, b):
        return model.thetab[b]

    def _mb_rule(self):
        def mb_rule(model, b):
            return sum(self.sx[b][k] * self._x_gap(model, k) for k in model.K) \
                + sum(self.sy[b][l] * self._y_gap(model, l) for l in model.L) == self._mb_rhs(model, b)
        return mb_rule

    def _vrs_rule(self):
        def vrs_rule(model):
            return sum(model.lamda[i2] for i2 in model.I2) == 1
        return vrs_rule

    def optimize(self, solver=OPT_DEFAULT):
        data2 = pd.DataFrame()
        obj, theta, thetax, thetay, thetab, lamda = {}, {}, {}, {}, {}, {}
        config = SolverConfig(solver=solver)
        for ind, problem in self._models.items():
            outcome = solve_pyomo_model(problem, config=config, strict=False)
            data2.loc[ind, "optimization_status"] = outcome.legacy_status
            if outcome.ok:
                obj[ind] = problem.objective()
                if self._include_var_undesirable:
                    theta[ind] = np.asarray(list(problem.theta[:].value))
                thetax[ind] = np.asarray(list(problem.thetax[:].value))
                thetay[ind] = np.asarray(list(problem.thetay[:].value))
                if self.b is not None:
                    thetab[ind] = np.asarray(list(problem.thetab[:].value))
                lamda[ind] = np.asarray(list(problem.lamda[:].value))
            else:
                obj[ind] = np.nan
                if self._include_var_undesirable:
                    theta[ind] = np.full(len(list(problem.B)), np.nan)
                thetax[ind] = np.full(len(list(problem.K)), np.nan)
                thetay[ind] = np.full(len(list(problem.L)), np.nan)
                if self.b is not None:
                    thetab[ind] = np.full(len(list(problem.B)), np.nan)
                lamda[ind] = np.full(len(list(problem.I2)), np.nan)

        obj_df = pd.DataFrame(obj, index=["obj"]).T
        pieces = []
        if self._include_var_undesirable:
            theta_df = pd.DataFrame(theta).T
            theta_df.columns = theta_df.columns.map(lambda b: "var Undesirable" if len(theta_df.columns) == 1 else "var Undesirable" + str(b))
            if len(theta_df.columns) == 1:
                theta_df.columns = ["var Undesirable"]
            pieces.append(theta_df)
        thetax_df = pd.DataFrame(thetax).T
        thetax_df.columns = thetax_df.columns.map(lambda x: "Input" + str(x) + "'s slack")
        thetay_df = pd.DataFrame(thetay).T
        thetay_df.columns = thetay_df.columns.map(lambda y: "Output" + str(y) + "'s slack")
        pieces.extend([thetax_df, thetay_df])
        if self.b is not None:
            thetab_df = pd.DataFrame(thetab).T
            thetab_df.columns = thetab_df.columns.map(lambda b: "Undesirable Output" + str(b) + "'s slack")
            pieces.append(thetab_df)
        theta_all = pd.concat(pieces, axis=1)
        lamda_df = pd.DataFrame(lamda).T
        lamda_df.columns = lamda_df.columns.map(lambda x: "lamda" + str(x))
        data3 = pd.concat([data2, obj_df], axis=1)
        data3 = pd.concat([data3, theta_all], axis=1)
        return data3

    def info(self, dmu="all"):
        if dmu == "all":
            if self._print_all_info:
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


class _MBThetaMixin:
    _has_theta = True
    _objective_sense = minimize
    _objective_component = "theta"
    _include_var_undesirable = True

    def _undesirable_rhs(self, model, b):
        return model.theta[b]

    def _mb_rhs(self, model, b):
        return self.b.loc[self.I0, self.bcol[b]] - model.theta[b]


class MBx(_MBThetaMixin, _PanelMBBase):
    """Input-side material-balance panel model."""
    _has_objx = True
    _has_objy = False
    _print_all_info = False


class MBx2(MBx):
    """Backward-compatible alias of :class:`MBx`."""
    _print_all_info = True


class MBxy(_MBThetaMixin, _PanelMBBase):
    """Input-output material-balance panel model with variable bad output target."""
    _has_objx = True
    _has_objy = True

    def _mb_rhs(self, model, b):
        return model.thetab[b] + self.b.loc[self.I0, self.bcol[b]] - model.theta[b]


class MBxy2(_PanelMBBase):
    """Input-output material-balance panel model with fixed bad output."""
    _has_objx = True
    _has_objy = True
    _objective_sense = maximize
    _objective_component = "thetab"


class MB2(_PanelMBBase):
    """Simplified material-balance panel model with observed input/output targets."""
    _has_objx = False
    _has_objy = False
    _objective_sense = maximize
    _objective_component = "thetab"

