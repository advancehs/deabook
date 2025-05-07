# import dependencies
from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, Constraint, log
from pyomo.core.expr.numvalue import NumericValue
import numpy as np
import pandas as pd
from ..constant import CET_ADDI, CET_MULT, FUN_PROD, FUN_COST, RTS_CRS, RTS_VRS, OPT_DEFAULT, OPT_LOCAL
from .tools import optimize_model, trans_list, to_2d_list
from . import weakCNLSbZG1
from . import tools


class weakCNLSbZG2(weakCNLSbZG1.weakCNLSbZG1):
    """weakCNLSbZ+G in iterative loop
    """

    def __init__(self, data, sent, z, cutactive, active,activeweak, cet=CET_ADDI, fun=FUN_PROD, rts=RTS_VRS,\
                 baseindex=None,refindex=None):
        """weakCNLSbZ+G model

        Args:
            sent (str): inputvars=outputvars: unoutputvars. e.g.: "K L = Y:CO2 "
            z (float, optional): Contextual variable(s). Defaults to None.
            cutactive (float or list): active concavity constraint.
            active (float or ndarray): violated concavity constraint.
            activeweak (float or ndarray): violated concavity constraint for weak disposibility.
            cet (String, optional): CET_ADDI (additive composite error term) or CET_MULT (multiplicative composite error term). Defaults to CET_ADDI.
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS.
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        self.outputvars, self.inputvars, self.unoutputvars, self.zvars = tools.assert_valid_yxbz_nog(sent, z)

        self.y, self.x, self.b, self.z, self.yref, self.xref, self.bref, self.zref \
            = tools.assert_valid_yxbz2(baseindex, refindex, data, \
                                       self.outputvars, self.inputvars, self.unoutputvars, self.zvars)
        self.xcol = self.x.columns
        self.ycol = self.y.columns
        self.bcol = self.b.columns
        self.zcol = self.z.columns if type(z) != type(None) else None

        self.cet = cet
        self.fun = fun
        self.rts = rts
        self.cutactive = cutactive
        self.active = to_2d_list(trans_list(active))
        self.activeweak = to_2d_list(trans_list(activeweak))

        # Initialize the CNLS model
        self.__model__ = ConcreteModel()

        # Initialize the sets
        self.__model__.I = Set(initialize=self.x.index)  ## I 是 被评价决策单元的数量
        self.__model__.I2 = Set(initialize=self.xref.index)  ## I2 是 参考决策单元的数量
        self.__model__.K = Set(initialize=range(len(self.x.iloc[0])))  ## K 是投入个数
        # self.__model__.J = Set(initialize=range(len(self.b.iloc[0])))  ## B 是 非期望产出个数
        self.__model__.L = Set(initialize=range(len(self.y.iloc[0])))  ## L 是产出个数 被评价单元和参考单元的K，L一样

        if type(self.z) != type(None):
            self.__model__.M = Set(initialize=range(len(self.z.iloc[0])))
            # Initialize the variables for z variable
            self.__model__.lamda = Var(self.__model__.M, doc='z coefficient')

        # Initialize the variables
        self.__model__.alpha = Var(self.__model__.I, doc='alpha')
        self.__model__.beta = Var(self.__model__.I,
                                  self.__model__.K,
                                  bounds=(0.0, None),
                                  doc='beta')
        self.__model__.gamma = Var(self.__model__.I,
                                   self.__model__.L,
                                   bounds=(0.0, None),
                                   doc='gamma')
        self.__model__.epsilon = Var(self.__model__.I, doc='residual')
        self.__model__.frontier = Var(self.__model__.I,
                                      bounds=(0.0, None),
                                      doc='estimated frontier')

        # Setup the objective function and constraints
        self.__model__.objective = Objective(rule=self._weakCNLSbZG1__objective_rule(),
                                             sense=minimize,
                                             doc='objective function')
        self.__model__.regression_rule = Constraint(self.__model__.I,
                                                    rule=self._weakCNLSbZG1__regression_rule(),
                                                    doc='regression equation')
        if self.cet == CET_MULT:
            self.__model__.log_rule = Constraint(self.__model__.I,
                                                 rule=self._weakCNLSbZG1__log_rule(),
                                                 doc='log-transformed regression equation')
        # self.__model__.afriat_rule = Constraint(self.__model__.I,
        #                                         rule=self._weakCNLSbZG1__afriat_rule(),
        #                                         doc='elementary Afriat approach')
        # self.__model__.disposability_rule = Constraint(self.__model__.I,
        #                                                 rule=self._weakCNLSbZG1__disposability_rule(),
        #                                                 doc='elementary weak disposibility')
        self.__model__.sweet_rule = Constraint(self.__model__.I,
                                               self.__model__.I2,
                                               rule=self._weakCNLSbZG1__sweet_rule(),
                                               doc='sweet spot approach')
        self.__model__.sweet_rule2 = Constraint(self.__model__.I,
                                               self.__model__.I2,
                                               rule=self.__sweet_rule2(),
                                               doc='sweet spot-2 approach')
        self.__model__.sweet_rule_weak2 = Constraint(self.__model__.I,
                                               self.__model__.I2,
                                               rule=self.__sweet_rule_weak2(),
                                               doc='sweet spot-2 approach for weak dis')
        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0


    def __sweet_rule2(self):
        """Return the proper sweet spot (step2) approach constraint"""
        if self.fun == FUN_PROD:
            __operator = NumericValue.__le__
        elif self.fun == FUN_COST:
            __operator = NumericValue.__ge__

        if self.rts == RTS_VRS:
            def sweet_rule2(model, i, h):
                if self.active[i][h]:
                    if i == h:
                        return Constraint.Skip
                    return __operator(model.alpha[i] \
                                  + sum(model.beta[i, k] * self.x.loc[i,self.xcol[k]] for k in model.K) \
                                  - sum(model.gamma[i, l] * self.y.loc[i,self.ycol[l]]  for l in model.L),
                                  model.alpha[h] \
                                  + sum(model.beta[h, k] * self.xref.loc[i,self.xcol[k]] for k in model.K) \
                                  - sum(model.gamma[h, l] * self.yref.loc[i,self.ycol[l]]  for l in model.L))
                return Constraint.Skip
            return sweet_rule2

        elif self.rts == RTS_CRS:
            def sweet_rule2(model, i, h):
                if self.active[i][h]:
                    if i == h:
                        return Constraint.Skip
                    return __operator(sum(model.beta[i, k] * self.x.loc[i,self.xcol[k]] for k in model.K) \
                                  - sum(model.gamma[i, l] * self.y.loc[i,self.ycol[l]]  for l in model.L),
                                 sum(model.beta[h, k] * self.xref.loc[i,self.xcol[k]] for k in model.K) \
                                  - sum(model.gamma[h, l] * self.yref.loc[i,self.ycol[l]]  for l in model.L))
                return Constraint.Skip
            return sweet_rule2

        raise ValueError("Undefined model parameters.")

    def __sweet_rule_weak2(self):
        """Return the proper sweet spot (step2) approach constraint for weak dis"""

        if self.rts == RTS_VRS:
            def sweet_rule_weak2(model, i, h):
                if self.activeweak[i][h]:
                    if i == h:
                        return Constraint.Skip
                    return model.alpha[i] \
                        + sum(model.beta[i, k] * self.xref.loc[h,self.xcol[k]] for k in model.K) >= 0
                return Constraint.Skip
            return sweet_rule_weak2

        elif self.rts == RTS_CRS:
            def sweet_rule_weak2(model, i, h):
                if self.activeweak[i][h]:
                    if i == h:
                        return Constraint.Skip
                    return sum(model.beta[i, k] * self.xref.loc[h,self.xcol[k]] for k in model.K) >= 0
                return Constraint.Skip
            return sweet_rule_weak2

        raise ValueError("Undefined model parameters.")
