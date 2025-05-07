# import dependencies
from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, Constraint, log
from ..constant import CET_ADDI, CET_MULT, FUN_PROD, FUN_COST, RTS_CRS, RTS_VRS, OPT_DEFAULT, OPT_LOCAL
from . import tools
from . import weakCNLSZG1
import numpy as np
import pandas as pd

class weakCNLSG1(weakCNLSZG1.weakCNLSZG1):
    """initial Group-VC-added weakCNLSZ (weakCNLS+G) model
    """

    def __init__(self, data,sent, cutactive, cet=CET_ADDI, fun=FUN_PROD, rts=RTS_VRS,baseindex=None,refindex=None):
        """CNLSZ+G model

        Args:
            sent (str): inputvars=outputvars:unoutputvars. e.g.: "K L = Y:CO2 "
            cutactive (float): active concavity constraint.
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

        # Initialize the CNLS model
        self.__model__ = ConcreteModel()

        # Initialize the sets
        self.__model__.I = Set(initialize=self.x.index)  ## I 是 被评价决策单元的数量
        self.__model__.I2 = Set(initialize=self.xref.index)  ## I2 是 参考决策单元的数量
        self.__model__.K = Set(initialize=range(len(self.x.iloc[0])))  ## K 是投入个数
        self.__model__.J = Set(initialize=range(len(self.b.iloc[0])))  ## B 是 非期望产出个数

        self.__model__.M = Set(initialize=range(len(self.z.iloc[0])))
        # Initialize the variables for z variable
        self.__model__.lamda = Var(self.__model__.M, doc='z coefficient')

        # Initialize the variables
        self.__model__.alpha = Var(self.__model__.I, doc='alpha')
        self.__model__.beta = Var(self.__model__.I,
                                  self.__model__.K,
                                  bounds=(0.0, None),
                                  doc='beta')
        self.__model__.delta = Var(self.__model__.I,
                                   self.__model__.J,
                                   bounds=(0.0, None),
                                   doc='delta')
        self.__model__.epsilon = Var(self.__model__.I, doc='residual')
        self.__model__.frontier = Var(self.__model__.I,
                                      bounds=(0.0, None),
                                      doc='estimated frontier')

        # Setup the objective function and constraints
        self.__model__.objective = Objective(rule=self._weakCNLSZG1__objective_rule(),
                                             sense=minimize,
                                             doc='objective function')
        self.__model__.regression_rule = Constraint(self.__model__.I,
                                                    rule=self.__regression_rule(),
                                                    doc='regression equation')
        if self.cet == CET_MULT:
            self.__model__.log_rule = Constraint(self.__model__.I,
                                                 rule=self._weakCNLSZG1__log_rule(),
                                                 doc='log-transformed regression equation')
        self.__model__.afriat_rule = Constraint(self.__model__.I,
                                                rule=self._weakCNLSZG1__afriat_rule(),
                                                doc='elementary Afriat approach')
        self.__model__.disposability_rule = Constraint(self.__model__.I,
                                                        rule=self._weakCNLSZG1__disposability_rule(),
                                                        doc='elementary weak disposibility')
        self.__model__.sweet_rule = Constraint(self.__model__.I,
                                               self.__model__.I,
                                               rule=self._weakCNLSZG1__sweet_rule(),
                                               doc='sweet spot approach')

        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0

    def __regression_rule(self):
        """Return the proper regression constraint"""
        if self.cet == CET_ADDI:
            if self.rts == RTS_VRS:
                def regression_rule(model, i):
                    return np.array((self.y.loc[i,])) \
                        == model.alpha[i] \
                        + sum(model.beta[i, k] * self.x.loc[i, self.xcol[k]] for k in model.K) \
                        + sum(model.delta[i, j] * self.b.loc[i, self.bcol[j]] for j in model.J) \
                        - sum(model.lamda[m] * self.z.loc[i, self.zcol[m]] for m in model.M) \
                        - model.epsilon[i]

                return regression_rule

            elif self.rts == RTS_CRS:
                def regression_rule(model, i):
                    return np.array((self.y.loc[i,])) == \
                        sum(model.beta[i, k] * self.x.loc[i, self.xcol[k]] for k in model.K) \
                        + sum(model.delta[i, j] * self.b.loc[i, self.bcol[j]] for j in model.J) \
                        - sum(model.lamda[m] * self.z.loc[i, self.zcol[m]] for m in model.M) \
                        - model.epsilon[i]

                return regression_rule

        elif self.cet == CET_MULT:
            def regression_rule(model, i):
                return log(np.array(self.y.loc[i,:]) ) == log(model.frontier[i] + 1) \
                        - sum(model.lamda[m] * self.z.loc[i,self.zcol[m]] for m in model.M) \
                        - model.epsilon[i]
            return regression_rule

        raise ValueError("Undefined model parameters.")


