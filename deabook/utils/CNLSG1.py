# import dependencies
from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, Constraint, log
from pyomo.core.expr.numvalue import NumericValue
import numpy as np
import pandas as pd
from ..constant import CET_ADDI, CET_MULT, FUN_PROD, FUN_COST, RTS_CRS, RTS_VRS1, OPT_DEFAULT, OPT_LOCAL
from .tools import optimize_model

from . import CNLSZG1


class CNLSG1(CNLSZG1.CNLSZG1):
    """initial Group-VC-added CNLS (CNLS+G) model
    """

    def __init__(self, y, x, cutactive, cet=CET_ADDI, fun=FUN_PROD, rts=RTS_VRS1):
        """CNLS+G model 1

        Args:
            y (float): output variable.
            x (float): input variables.
            cutactive (float): active concavity constraint.
            cet (String, optional): CET_ADDI (additive composite error term) or CET_MULT (multiplicative composite error term). Defaults to CET_ADDI.
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS.
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.x = x
        self.y = y
        self.cet = cet
        self.fun = fun
        self.rts = rts

        self.cutactive = cutactive

        # Initialize the CNLS model
        self.__model__ = ConcreteModel()

        # Initialize the sets
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))

        # Initialize the variables
        self.__model__.alpha = Var(self.__model__.I, doc='alpha')
        self.__model__.beta = Var(self.__model__.I,
                                  self.__model__.J,
                                  bounds=(0.0, None),
                                  doc='beta')
        self.__model__.epsilon = Var(self.__model__.I, doc='residual')
        self.__model__.frontier = Var(self.__model__.I,
                                      bounds=(0.0, None),
                                      doc='estimated frontier')

        # Setup the objective function and constraints
        self.__model__.objective = Objective(rule=self._CNLSZG1__objective_rule(),
                                             sense=minimize,
                                             doc='objective function')
        self.__model__.regression_rule = Constraint(self.__model__.I,
                                                    rule=self.__regression_rule(),
                                                    doc='regression equation')
        if self.cet == CET_MULT:
            self.__model__.log_rule = Constraint(self.__model__.I,
                                                 rule=self._CNLSZG1__log_rule(),
                                                 doc='log-transformed regression equation')
        self.__model__.afriat_rule = Constraint(self.__model__.I,
                                                rule=self._CNLSZG1__afriat_rule(),
                                                doc='elementary Afriat approach')
        self.__model__.sweet_rule = Constraint(self.__model__.I,
                                               self.__model__.I,
                                               rule=self._CNLSZG1__sweet_rule(),
                                               doc='sweet spot approach')

        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0




    def __regression_rule(self):
        """Return the proper regression constraint"""
        if self.cet == CET_ADDI:
            if self.rts == RTS_VRS1:
                def regression_rule(model, i):
                    return self.y[i] == model.alpha[i] \
                        + sum(model.beta[i, j] * self.x[i][j] for j in model.J) \
                        - model.epsilon[i]

                return regression_rule
            elif self.rts == RTS_CRS:
                def regression_rule(model, i):
                    return self.y[i] == sum(model.beta[i, j] * self.x[i][j] for j in model.J) \
                        - model.epsilon[i]

                return regression_rule
        elif self.cet == CET_MULT:

            def regression_rule(model, i):
                return log(
                    self.y[i]) == log(model.frontier[i] + 1) - model.epsilon[i]

            return regression_rule

        raise ValueError("Undefined model parameters.")

