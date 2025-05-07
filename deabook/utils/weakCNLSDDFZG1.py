# import dependencies
from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, Constraint, log
from pyomo.core.expr.numvalue import NumericValue
import numpy as np
import pandas as pd
from ..constant import CET_ADDI, FUN_PROD, FUN_COST, RTS_CRS, RTS_VRS, OPT_DEFAULT, OPT_LOCAL
from .tools import optimize_model
from . import tools


class weakCNLSDDFZG1():
    """initial Group-VC-added weakCNLSDDFZ (weakCNLSDDFZ+G) model
    """

    def __init__(self, data, sent, z, cutactive,gy=[1], gx=[1], gb=[1], fun=FUN_PROD, rts=RTS_VRS, \
                 baseindex=None, refindex=None):
        """CNLSDDFZ+G model

        Args:
            data (pandas.DataFrame): input pandas.
            sent (str): inputvars=outputvars: unoutputvars. e.g.: "K L = Y : CO2"
            z (float): Contextual variable(s). Defaults to None.
            cutactive (DataFrame, optional): active concavity constraint.
            gy (list, optional): output directional vector. Defaults to [1].
            gx (list, optional): input directional vector. Defaults to [1].
            gb (list, optional): undesirable output directional vector. Defaults to [1].
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS.
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        self.outputvars,self.inputvars,self.unoutputvars,self.zvars,self.gy, self.gx, self.gb \
            = tools.assert_valid_yxbz(sent,z,gy,gx,gb)
        self.y, self.x, self.b, self.z, self.yref, self.xref, self.bref, self.zref,self.referenceflag\
            = tools.assert_valid_yxbz2(baseindex,refindex,data,\
                                       self.outputvars,self.inputvars,self.unoutputvars,self.zvars)

        self.xcol = self.x.columns
        self.ycol = self.y.columns
        self.bcol = self.b.columns
        self.zcol = self.z.columns if type(z) != type(None) else None

        self.fun = fun
        self.rts = rts

        self.cutactive = cutactive
        # print("sssaaa111111",self.cutactive)
        # Initialize the CNLS model
        self.__model__ = ConcreteModel()

        # Initialize the sets
        self.__model__.I = Set(initialize=self.x.index)  ## I 是 被评价决策单元的数量
        self.__model__.I.pprint()
        if self.referenceflag:
            self.__model__.I2 = Set(initialize=self.xref.index)  ## I2 是 参考决策单元的数量
        self.__model__.K = Set(initialize=range(len(self.x.iloc[0])))  ## K 是投入个数
        self.__model__.L = Set(initialize=range(len(self.y.iloc[0])))  ## L 是产出个数 被评价单元和参考单元的K，L一样
        self.__model__.J = Set(initialize=range(len(self.b.iloc[0])))  ## B 是 非期望产出个数


        # Initialize the variables
        self.__model__.alpha = Var(self.__model__.I, doc='alpha')
        self.__model__.beta = Var(
            self.__model__.I, self.__model__.K, bounds=(0.0, None), doc='beta')
        self.__model__.gamma = Var(
            self.__model__.I, self.__model__.L, bounds=(0.0, None), doc='gamma')
        self.__model__.delta = Var(
            self.__model__.I, self.__model__.J, bounds=(0.0, None), doc='delta')
        self.__model__.epsilon = Var(self.__model__.I, doc='residuals')

        if type(self.z) != type(None):
            # Initialize the set of z
            self.__model__.M = Set(initialize=range(len(self.z.iloc[0])))
            # Initialize the variables for z variable
            self.__model__.lamda = Var(self.__model__.M, doc='z coefficient')

        # Setup the objective function and constraints
        self.__model__.objective = Objective(rule=self.__objective_rule(),
                                             sense=minimize,
                                             doc='objective function')
        self.__model__.regression_rule = Constraint(self.__model__.I,
                                                    rule=self.__regression_rule(),
                                                    doc='regression equation')
        self.__model__.regression_rule.pprint()
        self.__model__.translation_rule = Constraint(self.__model__.I,
                                                     rule=self.__translation_property(),
                                                     doc='translation property')
        self.__model__.translation_rule.pprint()

        self.__model__.afriat_rule = Constraint(self.__model__.I,
                                                rule=self.__afriat_rule(),
                                                doc='elementary Afriat approach')
        self.__model__.afriat_rule.pprint()
        self.__model__.disposability_rule = Constraint(self.__model__.I,
                                                        rule=self.__disposability_rule(),
                                                        doc='elementary weak disposibility')
        self.__model__.disposability_rule.pprint()

        self.__model__.sweet_rule = Constraint(self.__model__.I,
                                               self.__model__.I,
                                               rule=self.__sweet_rule(),
                                               doc='sweet spot approach')
        self.__model__.sweet_rule.pprint()
        self.__model__.sweet_weak_rule = Constraint(self.__model__.I,
                                               self.__model__.I,
                                               rule=self.__sweet_weak_rule(),
                                               doc='sweet spot approach for weak dis')
        self.__model__.sweet_weak_rule.pprint()
        # if self.referenceflag:
        #     self.__model__.afriat_ref_rule = Constraint(self.__model__.I,
        #                                                self.__model__.I2,
        #                                             rule=self.__afriat_ref_rule(),
        #                                             doc='afriat reference inequality')
        #     # self.__model__.afriat_ref_rule.pprint()
        #
        #     self.__model__.disposability_ref_rule = Constraint(self.__model__.I,
        #                                                self.__model__.I2,
        #                                             rule=self.__disposability_ref_rule(),
        #                                             doc='weak disposability reference ')
        #     self.__model__.disposability_ref_rule.pprint()
        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Optimize the function by requested method

        Args:
            email (string): The email address for remote optimization. It will optimize locally if OPT_LOCAL is given.
            solver (string): The solver chosen for optimization. It will optimize with default solver if OPT_DEFAULT is given.
        """
        # TODO(error/warning handling): Check problem status after optimization
        self.problem_status, self.optimization_status = optimize_model(
            self.__model__, email, CET_ADDI, solver)

    def __objective_rule(self):
        """Return the proper objective function"""
        def objective_rule(model):
            return sum(model.epsilon[i] ** 2 for i in model.I)
        return objective_rule

    def __regression_rule(self):
        """Return the proper regression constraint"""
        if self.rts == RTS_VRS:
            if type(self.z) != type(None):
                def regression_rule(model, i):
                    return sum(model.gamma[i, l] * self.y.loc[i,self.ycol[l]] for l in model.L) \
                        == model.alpha[i] \
                        + sum(model.beta[i, k] * self.x.loc[i,self.xcol[k]]  for k in model.K) \
                        + sum(model.delta[i, j] * self.b.loc[i,self.bcol[j]] for j in model.J) \
                        - sum(model.lamda[m] * self.z.loc[i,self.zcol[m]] for m in model.M) \
                        - model.epsilon[i]
                return regression_rule

            elif type(self.z) == type(None):
                def regression_rule(model, i):
                    return sum(model.gamma[i, l] * self.y.loc[i,self.ycol[l]] for l in model.L) \
                        == model.alpha[i] \
                        + sum(model.beta[i, k] * self.x.loc[i,self.xcol[k]] for k in model.K) \
                        + sum(model.delta[i, j] * self.b.loc[i,self.bcol[j]] for j in model.J) \
                        - model.epsilon[i]
                return regression_rule

        elif self.rts == RTS_CRS:
            if type(self.z) != type(None):
                def regression_rule(model, i):
                    return sum(model.gamma[i, l] * self.y.loc[i,self.ycol[l]] for l in model.L) \
                        == sum(model.beta[i, k] * self.x.loc[i,self.xcol[k]] for k in model.K) \
                        + sum(model.delta[i, j] * self.b.loc[i,self.bcol[j]] for j in model.J) \
                        - sum(model.lamda[m] * self.z.loc[i,self.zcol[m]] for m in model.M) \
                        - model.epsilon[i]
                return regression_rule

            elif type(self.z) == type(None):
                def regression_rule(model, i):
                    return sum(model.gamma[i, l] * self.y.loc[i,self.ycol[l]] for l in model.L) \
                        == sum(model.beta[i, k] * self.x.loc[i,self.xcol[k]] for k in model.K) \
                        + sum(model.delta[i, j] * self.b.loc[i,self.bcol[j]] for j in model.J) \
                        - model.epsilon[i]
                return regression_rule

        raise ValueError("Undefined model parameters.")

    def __translation_property(self):
        """Return the proper translation property"""
        def translation_rule(model, i):
            return sum(model.beta[i, k] * self.x.loc[i,self.xcol[k]] for k in model.K) \
                + sum(model.gamma[i, l] * self.y.loc[i,self.ycol[l]] for l in model.L) \
                + sum(model.delta[i, j] * self.b.loc[i,self.bcol[j]] for j in model.J) == 1

        return translation_rule

    def __afriat_rule(self):
        """Return the proper elementary Afriat approach constraint"""
        if self.fun == FUN_PROD:
            __operator = NumericValue.__le__
        elif self.fun == FUN_COST:
            __operator = NumericValue.__ge__

        if self.rts == RTS_VRS:
            def afriat_rule(model, i):
                return __operator(
                    model.alpha[i] + sum(model.beta[i, k] *  self.x.loc[i,self.xcol[k]]  for k in model.K) \
                               + sum(model.delta[i, j] *  self.b.loc[i,self.bcol[j]] for j in model.J) \
                                - sum(model.gamma[i, l] * self.y.loc[i,self.ycol[l]] for l in model.L),
                    model.alpha[self.__model__.I.nextw(i)] \
           + sum(model.beta[self.__model__.I.nextw(i), k] * self.x.loc[i,self.xcol[k]]  for k in model.K) \
             + sum(model.delta[self.__model__.I.nextw(i), j] * self.b.loc[i,self.bcol[j]] for j in model.J)\
             - sum(model.gamma[self.__model__.I.nextw(i), l] * self.y.loc[i,self.ycol[l]] for l in model.L))
            return afriat_rule

        elif self.rts == RTS_CRS:
            def afriat_rule(model, i):
                return __operator( sum(model.beta[i, k] *  self.x.loc[i,self.xcol[k]]  for k in model.K) \
                               + sum(model.delta[i, j] *  self.b.loc[i,self.bcol[j]] for j in model.J) \
                                - sum(model.gamma[i, l] * self.y.loc[i,self.ycol[l]] for l in model.L),
            sum(model.beta[self.__model__.I.nextw(i), k] * self.x.loc[i,self.xcol[k]]  for k in model.K) \
             + sum(model.delta[self.__model__.I.nextw(i), j] * self.b.loc[i,self.bcol[j]] for j in model.J)\
             - sum(model.gamma[self.__model__.I.nextw(i), l] * self.y.loc[i,self.ycol[l]] for l in model.L))
            return afriat_rule

        raise ValueError("Undefined model parameters.")

    def __disposability_rule(self):
        """Return the proper elementary weak disposability constraint"""
        if self.rts == RTS_VRS:
            def disposability_rule(model, i):
                return model.alpha[self.__model__.I.nextw(i)] \
                    + sum(model.beta[self.__model__.I.nextw(i), k] * self.x.loc[i,self.xcol[k]] for k in model.K) >= 0
            return disposability_rule

        elif self.rts == RTS_CRS:
            def disposability_rule(model, i):
                return sum(model.beta[self.__model__.I.nextw(i), k] * self.x.loc[i,self.xcol[k]] for k in model.K) >= 0
            return disposability_rule
        raise ValueError("Undefined model parameters.")

    def __sweet_rule(self ):
        """Return the proper sweet spot approach constraint"""
        if self.fun == FUN_PROD:
            __operator = NumericValue.__le__
        elif self.fun == FUN_COST:
            __operator = NumericValue.__ge__

        if self.rts == RTS_VRS:
            def sweet_rule(model, i, h):
                if self.cutactive.loc[i,h]:
                    if i == h:
                        return Constraint.Skip
                    return __operator(model.alpha[i] \
                              + sum(model.beta[i, k] * self.x.loc[i,self.xcol[k]]  for k in model.K) \
                              + sum(model.delta[i, j] * self.b.loc[i,self.bcol[j]] for j in model.J) \
                              - sum(model.gamma[i, l] * self.y.loc[i,self.ycol[l]] for l in model.L),
                              model.alpha[h] \
                              + sum(model.beta[h, k] * self.x.loc[i,self.xcol[k]]  for k in model.K) \
                              + sum(model.delta[h, j] * self.b.loc[i,self.bcol[j]] for j in model.J) \
                              - sum(model.gamma[h, l] * self.y.loc[i,self.ycol[l]] for l in model.L))
                return Constraint.Skip
            return sweet_rule

        elif self.rts == RTS_CRS:
            def sweet_rule(model, i, h):
                if self.cutactive.loc[i,h]:
                    if i == h:
                        return Constraint.Skip
                    return __operator(sum(model.beta[i, k] * self.x.loc[i,self.xcol[k]]  for k in model.K) \
                              + sum(model.delta[i, j] * self.b.loc[i,self.bcol[j]] for j in model.J) \
                              - sum(model.gamma[i, l] *  self.y.loc[i,self.ycol[l]] for l in model.L),
                              sum(model.beta[h, k] * self.x.loc[i,self.xcol[k]]  for k in model.K) \
                              + sum(model.delta[h, j] * self.b.loc[i,self.bcol[j]] for j in model.J) \
                              - sum(model.gamma[h, l] * self.y.loc[i,self.ycol[l]] for l in model.L))
                return Constraint.Skip
            return sweet_rule

        raise ValueError("Undefined model parameters.")

    def __sweet_weak_rule(self):
        """Return the proper sweet spot approach constraint"""

        if self.rts == RTS_VRS:
            def sweet_weaak_rule(model, i, h):
                if self.cutactive.loc[i,h]:
                    if i == h:
                        return Constraint.Skip
                    return model.alpha[i] + sum(model.beta[i, k] * self.x.loc[h, self.xcol[k]] for k in model.K) >= 0
                return Constraint.Skip
            return sweet_weaak_rule

        elif self.rts == RTS_CRS:
            def sweet_weaak_rule(model, i, h):
                if self.cutactive.loc[i,h]:
                    if i == h:
                        return Constraint.Skip
                    return sum(model.beta[i, k] * self.x.loc[h, self.xcol[k]] for k in model.K) >= 0
                return Constraint.Skip
            return sweet_weaak_rule

        raise ValueError("Undefined model parameters.")

    def get_alpha(self):
        """Return alpha value by array"""
        if self.optimization_status == 0:
            self.optimize()
        alpha = pd.Series(self.__model__.alpha.extract_values())
        return alpha

    def get_beta(self):
        """Return beta value by array"""
        if self.optimization_status == 0:
            self.optimize()
        beta = pd.Series(self.__model__.beta.extract_values(),index=self.__model__.beta.extract_values().keys())
        # if the series is multi-indexed we need to unstack it...
        if type(beta.index[0]) == tuple:  # it is multi-indexed
            beta = beta.unstack(level=1)
        else:
            beta = pd.DataFrame(beta)  # force transition from Series -> df
        # multi-index the columns
        beta.columns = map(lambda x: "beta"+str(x) ,beta.columns)
        return beta

    def get_delta(self):
        """Return delta value by array"""
        if self.optimization_status == 0:
            self.optimize()
        delta = pd.Series(self.__model__.delta.extract_values(),index=self.__model__.delta.extract_values().keys())
        # if the series is multi-indexed we need to unstack it...
        if type(delta.index[0]) == tuple:  # it is multi-indexed
            delta = delta.unstack(level=1)
        else:
            delta = pd.DataFrame(delta)  # force transition from Series -> df
        # multi-index the columns
        delta.columns = map(lambda x: "delta"+str(x) ,delta.columns)
        return delta

    def get_gamma(self):
        """Return delta value by array"""
        if self.optimization_status == 0:
            self.optimize()
        gamma = pd.Series(self.__model__.gamma.extract_values(),index=self.__model__.gamma.extract_values().keys())
        # if the series is multi-indexed we need to unstack it...
        if type(gamma.index[0]) == tuple:  # it is multi-indexed
            gamma = gamma.unstack(level=1)
        else:
            gamma = pd.DataFrame(gamma)  # force transition from Series -> df
        # multi-index the columns
        gamma.columns = map(lambda x: "delta"+str(x) ,gamma.columns)
        return gamma
