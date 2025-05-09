# import dependencies
from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, Constraint, log
from pyomo.core.expr.numvalue import NumericValue
import numpy as np
import pandas as pd
from ..constant import CET_ADDI, CET_MULT, FUN_PROD, FUN_COST, RTS_CRS, RTS_VRS1, OPT_DEFAULT, OPT_LOCAL
from .tools import optimize_model

class CQRZG1():
    """initial Group-VC-added CQRZ (CQRZ+G) model
    """

    def __init__(self, y, x, z, tau, cutactive, cet=CET_ADDI, fun=FUN_PROD, rts=RTS_VRS1):
        """CQRZ+G model

        Args:
            y (float): output variable.
            x (float): input variables.
            z (float, optional): Contextual variable(s). Defaults to None.
            tau (float): quantile.
            cutactive (float): active concavity constraint.
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
        self.__model__.lamda = Var(self.__model__.K, doc='Zvalue')
        self.__model__.epsilon_plus = Var(
            self.__model__.I, bounds=(0.0, None), doc='positive error term')
        self.__model__.epsilon_minus = Var(
            self.__model__.I, bounds=(0.0, None), doc='negative error term')
        self.__model__.frontier = Var(self.__model__.I,
                                      bounds=(0.0, None),
                                      doc='estimated frontier')

        # Setup the objective function and constraints
        self.__model__.objective = Objective(rule=self.__objective_rule(),
                                             sense=minimize,
                                             doc='objective function')

        self.__model__.regression_rule = Constraint(self.__model__.I,
                                                    rule=self.__regression_rule(),
                                                    doc='regression equation')
        if self.cet == CET_MULT:
            self.__model__.log_rule = Constraint(self.__model__.I,
                                                 rule=self.__log_rule(),
                                                 doc='log-transformed regression equation')
        self.__model__.afriat_rule = Constraint(self.__model__.I,
                                                rule=self.__afriat_rule(),
                                                doc='elementary Afriat approach')
        self.__model__.sweet_rule = Constraint(self.__model__.I,
                                               self.__model__.I,
                                               rule=self.__sweet_rule(),
                                               doc='sweet spot approach')

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
            self.__model__, email, self.cet, solver)

    def __objective_rule(self):
        """Return the proper objective function"""

        def objective_rule(model):
            return (1 - self.tau) * sum(model.epsilon_plus[i] for i in model.I) \
                + self.tau * sum(model.epsilon_minus[i] for i in model.I)
        return objective_rule


    def __regression_rule(self):
        """Return the proper regression constraint"""
        if self.cet == CET_ADDI:
            if self.rts == RTS_VRS1:

                def regression_rule(model, i):
                    return self.y[i] == model.alpha[i] \
                        + sum(model.beta[i, j] * self.x[i][j] for j in model.J) \
                        - sum(model.lamda[k] * self.z[i][k] for k in model.K) \
                        - model.epsilon_minus[i] + model.epsilon_plus[i]

                return regression_rule
            elif self.rts == RTS_CRS:

                def regression_rule(model, i):
                    return self.y[i] == sum(model.beta[i, j] * self.x[i][j] for j in model.J) \
                        + sum(model.lamda[k] * self.z[i][k] for k in model.K) \
                        - model.epsilon_minus[i] + model.epsilon_plus[i]
                return regression_rule

        elif self.cet == CET_MULT:

            def regression_rule(model, i):
                return log(self.y[i]) == log(model.frontier[i] + 1) \
                    - sum(model.lamda[k] * self.z[i][k] for k in model.K) \
                    - model.epsilon_minus[i] + model.epsilon_plus[i]
            return regression_rule

        raise ValueError("Undefined model parameters.")

    def __log_rule(self):
        """Return the proper log constraint"""
        if self.cet == CET_MULT:
            if self.rts == RTS_VRS1:
                def log_rule(model, i):
                    return model.frontier[i] == model.alpha[i] + sum(
                        model.beta[i, j] * self.x[i][j] for j in model.J) - 1
                return log_rule

            elif self.rts == RTS_CRS:
                def log_rule(model, i):
                    return model.frontier[i] == sum(
                        model.beta[i, j] * self.x[i][j] for j in model.J) - 1
                return log_rule

        raise ValueError("Undefined model parameters.")

    def __afriat_rule(self):
        """Return the proper elementary Afriat approach constraint"""
        if self.fun == FUN_PROD:
            __operator = NumericValue.__le__
        elif self.fun == FUN_COST:
            __operator = NumericValue.__ge__

        if self.cet == CET_ADDI:
            if self.rts == RTS_VRS1:
                def afriat_rule(model, i):
                    return __operator(
                        model.alpha[i] + sum(model.beta[i, j] * self.x[i][j]
                                             for j in model.J),
                        model.alpha[self.__model__.I.nextw(i)] +
                        sum(model.beta[self.__model__.I.nextw(i), j] * self.x[i][j]
                            for j in model.J))
                return afriat_rule

            elif self.rts == RTS_CRS:
                def afriat_rule(model, i):
                    return __operator(
                        sum(model.beta[i, j] * self.x[i][j]
                            for j in model.J),
                        sum(model.beta[self.__model__.I.nextw(i), j] * self.x[i][j]
                            for j in model.J))
                return afriat_rule

        elif self.cet == CET_MULT:
            if self.rts == RTS_VRS1:
                def afriat_rule(model, i):
                    return __operator(
                        model.alpha[i] + sum(model.beta[i, j] * self.x[i][j]
                                             for j in model.J),
                        model.alpha[self.__model__.I.nextw(i)] +
                        sum(model.beta[self.__model__.I.nextw(i), j] * self.x[i][j]
                            for j in model.J))
                return afriat_rule

            elif self.rts == RTS_CRS:
                def afriat_rule(model, i):
                    return __operator(
                        sum(model.beta[i, j] * self.x[i][j] for j in model.J),
                        sum(model.beta[self.__model__.I.nextw(i), j] * self.x[i][j] for j in model.J))
                return afriat_rule

        raise ValueError("Undefined model parameters.")

    def __sweet_rule(self, ):
        """Return the proper sweet spot approach constraint"""
        if self.fun == FUN_PROD:
            __operator = NumericValue.__le__
        elif self.fun == FUN_COST:
            __operator = NumericValue.__ge__

        if self.cet == CET_ADDI:
            if self.rts == RTS_VRS1:
                def sweet_rule(model, i, h):
                    if self.cutactive[i][h]:
                        if i == h:
                            return Constraint.Skip
                        return __operator(model.alpha[i] + sum(model.beta[i, j] * self.x[i][j]
                                                               for j in model.J),
                                          model.alpha[h] + sum(model.beta[h, j] * self.x[i][j]
                                                               for j in model.J))
                    return Constraint.Skip
                return sweet_rule

            elif self.rts == RTS_CRS:
                def sweet_rule(model, i, h):
                    if self.cutactive[i][h]:
                        if i == h:
                            return Constraint.Skip
                        return __operator(sum(model.beta[i, j] * self.x[i][j]
                                              for j in model.J),
                                          sum(model.beta[h, j] * self.x[i][j]
                                              for j in model.J))
                    return Constraint.Skip
                return sweet_rule

        elif self.cet == CET_MULT:
            if self.rts == RTS_VRS1:
                def sweet_rule(model, i, h):
                    if self.cutactive[i][h]:
                        if i == h:
                            return Constraint.Skip
                        return __operator(model.alpha[i] + sum(model.beta[i, j] * self.x[i][j]
                                                               for j in model.J),
                                          model.alpha[h] + sum(model.beta[h, j] * self.x[i][j]
                                                               for j in model.J))
                    return Constraint.Skip
                return sweet_rule

            elif self.rts == RTS_CRS:
                def sweet_rule(model, i, h):
                    if self.cutactive[i][h]:
                        if i == h:
                            return Constraint.Skip
                        return __operator(sum(model.beta[i, j] * self.x[i][j] for j in model.J),
                                          sum(model.beta[h, j] * self.x[i][j] for j in model.J))
                    return Constraint.Skip
                return sweet_rule

        raise ValueError("Undefined model parameters.")

    def get_alpha(self):
        """Return alpha value by array"""
        if self.optimization_status == 0:
            self.optimize()
        alpha = list(self.__model__.alpha[:].value)
        return np.asarray(alpha)

    def get_beta(self):
        """Return beta value by array"""
        if self.optimization_status == 0:
            self.optimize()
        beta = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.beta),
                                                          list(self.__model__.beta[:, :].value))])
        beta = pd.DataFrame(beta, columns=['Name', 'Key', 'Value'])
        beta = beta.pivot(index='Name', columns='Key', values='Value')
        return beta.to_numpy()


class CERZG1(CQRZG1):
    """initial Group-VC-added CERZ (CERZ+G) model
    """

    def __init__(self, y, x, z, tau, cutactive, cet=CET_ADDI, fun=FUN_PROD, rts=RTS_VRS1):
        """CERZ+G model

        Args:
            y (float): output variable.
            x (float): input variables.
            z (float, optional): Contextual variable(s). Defaults to None.
            tau (float): expectile.
            cutactive (float): active concavity constraint.
            cet (String, optional): CET_ADDI (additive composite error term) or CET_MULT (multiplicative composite error term). Defaults to CET_ADDI.
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS1.
        """
        super().__init__(y, x, z, tau, cutactive, cet, fun, rts)
        self.__model__.objective.deactivate()
        self.__model__.squared_objective = Objective(
            rule=self.__squared_objective_rule(), sense=minimize, doc='squared objective rule')

    def __squared_objective_rule(self):
        def squared_objective_rule(model):
            return (1 - self.tau) * sum(model.epsilon_plus[i] ** 2 for i in model.I) \
                + self.tau * sum(model.epsilon_minus[i] ** 2 for i in model.I)

        return squared_objective_rule
