# import dependencies
from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, Constraint, log
from pyomo.core.expr.numvalue import NumericValue
import numpy as np
import pandas as pd
from ..constant import CET_ADDI, CET_MULT, FUN_PROD, FUN_COST, RTS_CRS, RTS_VRS, OPT_DEFAULT, OPT_LOCAL
from .tools import  trans_list, to_2d_list
from . import weakCNLSxZG1,weakCNLSxZG2,weakCNLSxG1


class weakCNLSxG2(weakCNLSxZG2.weakCNLSxZG2,weakCNLSxG1.weakCNLSxG1,weakCNLSxZG1.weakCNLSxZG1):
    """weakCNLSx+G in iterative loop
    """

    def __init__(self, y, x, b, cutactive, active,activeweak, cet=CET_ADDI, fun=FUN_PROD, rts=RTS_VRS):
        """weakCNLSx+G model

        Args:
            y (float): output variable.
            x (float): input variables.
            b (float): undersiable variables.
            cutactive (float or list): active concavity constraint.
            active (float or ndarray): violated concavity constraint.
            activeweak (float or ndarray): violated concavity constraint for weak disposibility.
            cet (String, optional): CET_ADDI (additive composite error term) or CET_MULT (multiplicative composite error term). Defaults to CET_ADDI.
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS.
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.x = x
        self.y = y
        self.b = b
        self.cet = cet
        self.fun = fun
        self.rts = rts

        self.cutactive = cutactive
        self.active = to_2d_list(trans_list(active))
        self.activeweak = to_2d_list(trans_list(activeweak))

        # Initialize the CNLS model
        self.__model__ = ConcreteModel()

        # Initialize the sets
        self.__model__.I = Set(initialize=range(len(self.x)))
        self.__model__.J = Set(initialize=range(len(self.y[0])))
        self.__model__.L = Set(initialize=range(len(self.b[0])))

        # Initialize the variables
        self.__model__.alpha = Var(self.__model__.I, doc='alpha')
        self.__model__.delta = Var(self.__model__.I,
                                  self.__model__.J,
                                  bounds=(0.0, None),
                                  doc='delta')
        self.__model__.gamma = Var(self.__model__.I,
                                   self.__model__.L,
                                   bounds=(0.0, None),
                                   doc='gamma')
        self.__model__.epsilon = Var(self.__model__.I, doc='residual')
        self.__model__.frontier = Var(self.__model__.I,
                                      bounds=(0.0, None),
                                      doc='estimated frontier')

        # Setup the objective function and constraints
        self.__model__.objective = Objective(rule=self._weakCNLSxZG1__objective_rule(),
                                             sense=minimize,
                                             doc='objective function')
        self.__model__.regression_rule = Constraint(self.__model__.I,
                                                    rule=self._weakCNLSxG1__regression_rule(),
                                                    doc='regression equation')
        if self.cet == CET_MULT:
            self.__model__.log_rule = Constraint(self.__model__.I,
                                                 rule=self._weakCNLSxZG1__log_rule(),
                                                 doc='log-transformed regression equation')
        self.__model__.afriat_rule = Constraint(self.__model__.I,
                                                rule=self._weakCNLSxZG1__afriat_rule(),
                                                doc='elementary Afriat approach')
        self.__model__.disposability_rule = Constraint(self.__model__.I,
                                                        rule=self._weakCNLSxZG1__disposability_rule(),
                                                        doc='elementary weak disposibility')
        self.__model__.sweet_rule = Constraint(self.__model__.I,
                                               self.__model__.I,
                                               rule=self._weakCNLSxZG1__sweet_rule(),
                                               doc='sweet spot approach')
        self.__model__.sweet_rule2 = Constraint(self.__model__.I,
                                               self.__model__.I,
                                               rule=self._weakCNLSxZG2__sweet_rule2(),
                                               doc='sweet spot-2 approach')
        self.__model__.sweet_rule_weak2 = Constraint(self.__model__.I,
                                               self.__model__.I,
                                               rule=self._weakCNLSxZG2__sweet_rule_weak2(),
                                               doc='sweet spot-2 approach for weak dis')
        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0


    def __regression_rule(self):
        """Return the proper regression constraint"""
        if self.cet == CET_ADDI:
            if self.rts == RTS_VRS:

                def regression_rule(model, i):
                    return self.x[i][0] == - model.alpha[i] \
                        + sum(model.gamma[i, j] * self.y[i][j] for j in model.J) \
                        - sum(model.delta[i, l] * self.b[i][l] for l in model.L) \
                        + model.epsilon[i]

                return regression_rule
            elif self.rts == RTS_CRS:

                def regression_rule(model, i):
                    return self.x[i][0] == sum(model.gamma[i, j] * self.y[i][j] for j in model.J) \
                        - sum(model.delta[i, l] * self.b[i][l] for l in model.L) \
                        + model.epsilon[i]

                return regression_rule

        elif self.cet == CET_MULT:

            def regression_rule(model, i):
                return log(self.x[i][0]) == - log(model.frontier[i] + 1) \
                     + model.epsilon[i]

            return regression_rule

        raise ValueError("Undefined model parameters.")
