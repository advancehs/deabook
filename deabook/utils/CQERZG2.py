# import dependencies
from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, Constraint, log
from pyomo.core.expr.numvalue import NumericValue
import numpy as np
import pandas as pd
from ..constant import CET_ADDI, CET_MULT, FUN_PROD, FUN_COST, RTS_CRS, RTS_VRS1, OPT_DEFAULT, OPT_LOCAL
from .tools import optimize_model, trans_list, to_2d_list
from . import CQERZG1


class CQRZG2(CQERZG1.CQRZG1):
    """CQRZ+G in iterative loop
    """

    def __init__(self, y, x, z, tau, cutactive, active, cet=CET_ADDI, fun=FUN_PROD, rts=RTS_VRS1):
        """CQRZ+G model

        Args:
            y (float): output variable.
            x (float): input variables.
            z (float, optional): Contextual variable(s). Defaults to None.
            tau (float): quantile.
            cutactive (float): active concavity constraint.
            active (float): violated concavity constraint.
            cet (String, optional): CET_ADDI (additive composite error term) or CET_MULT (multiplicative composite error term). Defaults to CET_ADDI.
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS1.
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.x = x
        self.y = y
        self.z = z
        self.tau = tau
        self.cet = cet
        self.fun = fun
        self.rts = rts

        self.cutactive = cutactive
        self.active = to_2d_list(trans_list(active))

        # Initialize the CNLS model
        self.__model__ = ConcreteModel()

        # Initialize the sets
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))
        self.__model__.K = Set(initialize=range(len(self.z[0])))

        # Initialize the variables
        self.__model__.alpha = Var(self.__model__.I, doc='alpha')
        self.__model__.beta = Var(self.__model__.I,
                                  self.__model__.J,
                                  bounds=(0.0, None),
                                  doc='beta')
        self.__model__.lamda = Var(self.__model__.K, doc='zvalue')
        self.__model__.epsilon_plus = Var(
            self.__model__.I, bounds=(0.0, None), doc='positive error term')
        self.__model__.epsilon_minus = Var(
            self.__model__.I, bounds=(0.0, None), doc='negative error term')
        self.__model__.frontier = Var(self.__model__.I,
                                      bounds=(0.0, None),
                                      doc='estimated frontier')

        # Setup the objective function and constraints
        self.__model__.objective = Objective(rule=self._CQRZG1__objective_rule(),
                                             sense=minimize,
                                             doc='objective function')

        self.__model__.regression_rule = Constraint(self.__model__.I,
                                                    rule=self._CQRZG1__regression_rule(),
                                                    doc='regression equation')
        if self.cet == CET_MULT:
            self.__model__.log_rule = Constraint(self.__model__.I,
                                                 rule=self._CQRZG1__log_rule(),
                                                 doc='log-transformed regression equation')
        self.__model__.afriat_rule = Constraint(self.__model__.I,
                                                rule=self._CQRZG1__afriat_rule(),
                                                doc='elementary Afriat approach')
        self.__model__.sweet_rule = Constraint(self.__model__.I,
                                               self.__model__.I,
                                               rule=self._CQRZG1__sweet_rule(),
                                               doc='sweet spot approach')
        self.__model__.sweet_rule2 = Constraint(self.__model__.I,
                                                self.__model__.I,
                                                rule=self.__sweet_rule2(),
                                                doc='sweet spot-2 approach')

        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0


    def __sweet_rule2(self, ):
        """Return the proper sweet spot (step2) approach constraint"""
        if self.fun == FUN_PROD:
            __operator = NumericValue.__le__
        elif self.fun == FUN_COST:
            __operator = NumericValue.__ge__

        if self.cet == CET_ADDI:
            if self.rts == RTS_VRS1:

                def sweet_rule2(model, i, h):
                    if self.active[i][h]:
                        if i == h:
                            return Constraint.Skip
                        return __operator(model.alpha[i] + sum(model.beta[i, j] * self.x[i][j]
                                                               for j in model.J),
                                          model.alpha[h] + sum(model.beta[h, j] * self.x[i][j]
                                                               for j in model.J))
                    return Constraint.Skip

                return sweet_rule2
            elif self.rts == RTS_CRS:
                def sweet_rule2(model, i, h):
                    if self.active[i][h]:
                        if i == h:
                            return Constraint.Skip
                        return __operator(model.alpha[i] + sum(model.beta[i, j] * self.x[i][j]
                                                               for j in model.J),
                                          model.alpha[h] + sum(model.beta[h, j] * self.x[i][j]
                                                               for j in model.J))
                    return Constraint.Skip

                return sweet_rule2
        elif self.cet == CET_MULT:
            if self.rts == RTS_VRS1:

                def sweet_rule2(model, i, h):
                    if self.active[i][h]:
                        if i == h:
                            return Constraint.Skip
                        return __operator(model.alpha[i] + sum(model.beta[i, j] * self.x[i][j]
                                                               for j in model.J),
                                          model.alpha[h] + sum(model.beta[h, j] * self.x[i][j]
                                                               for j in model.J))
                    return Constraint.Skip

                return sweet_rule2
            elif self.rts == RTS_CRS:

                def sweet_rule2(model, i, h):
                    if self.active[i][h]:
                        if i == h:
                            return Constraint.Skip
                        return __operator(sum(model.beta[i, j] * self.x[i][j] for j in model.J),
                                          sum(model.beta[h, j] * self.x[i][j] for j in model.J))
                    return Constraint.Skip

                return sweet_rule2

        raise ValueError("Undefined model parameters.")


class CERZG2(CQRZG2):
    """CERZ+G in iterative loop
    """

    def __init__(self, y, x, z, tau, cutactive, active, cet=CET_ADDI, fun=FUN_PROD, rts=RTS_VRS1):
        """CERZ+G model

        Args:
            y (float): output variable.
            x (float): input variables.
            z (float, optional): Contextual variable(s). Defaults to None.
            tau (float): expectile.
            cutactive (float): active concavity constraint.
            active (float): violated concavity constraint.
            cet (String, optional): CET_ADDI (additive composite error term) or CET_MULT (multiplicative composite error term). Defaults to CET_ADDI.
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS1.
        """
        super().__init__(y, x, z, tau, cutactive, active, cet, fun, rts)
        self.__model__.objective.deactivate()
        self.__model__.squared_objective = Objective(
            rule=self.__squared_objective_rule(), sense=minimize, doc='squared objective rule')

    def __squared_objective_rule(self):
        def squared_objective_rule(model):
            return (1 - self.tau) * sum(model.epsilon_plus[i] ** 2 for i in model.I) \
                + self.tau * sum(model.epsilon_minus[i] ** 2 for i in model.I)
        return squared_objective_rule
