# import dependencies
from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, Constraint, log
from pyomo.core.expr.numvalue import NumericValue
import numpy as np
import pandas as pd

from .constant import CET_ADDI, CET_MULT, FUN_PROD, FUN_COST, RTS_CRS, RTS_VRS1, OPT_LOCAL, OPT_DEFAULT
from .utils import tools, interpolation


class CQR:
    """Convex quantile regression (CQR)
    """

    def __init__(self, y, x, tau, z=None, cet=CET_ADDI, fun=FUN_PROD, rts=RTS_VRS1):
        """CQR model

        Args:
            y (float): output variable.
            x (float): input variables.
            tau (float): quantile.
            z (float, optional): Contextual variable(s). Defaults to None.
            cet (String, optional): CET_ADDI (additive composite error term) or CET_MULT (multiplicative composite error term). Defaults to CET_ADDI.
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS1.
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.y, self.x, self.z = tools.assert_valid_basic_data(y, x, z)
        self.tau = tau
        self.cet = cet
        self.fun = fun
        self.rts = rts

        # Initialize the CQR model
        self.__model__ = ConcreteModel()

        if type(self.z) != type(None):
            # Initialize the set of z
            self.__model__.K = Set(initialize=range(len(self.z[0])))

            # Initialize the variables for z variable
            self.__model__.lamda = Var(self.__model__.K, doc='z coefficient')

        # Initialize the sets
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))

        # Initialize the variables
        self.__model__.alpha = Var(self.__model__.I, doc='alpha')
        self.__model__.beta = Var(self.__model__.I,
                                  self.__model__.J,
                                  bounds=(0.0, None),
                                  doc='beta')
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
                                                self.__model__.I,
                                                rule=self.__afriat_rule(),
                                                doc='afriat inequality')

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
        self.problem_status, self.optimization_status = tools.optimize_model(
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
                if type(self.z) != type(None):
                    def regression_rule(model, i):
                        return self.y[i] == model.alpha[i] \
                            + sum(model.beta[i, j] * self.x[i][j] for j in model.J) \
                            - sum(model.lamda[k] * self.z[i][k]
                                  for k in model.K) + model.epsilon_plus[i] - model.epsilon_minus[i]
                    return regression_rule

                def regression_rule(model, i):
                    return self.y[i] == model.alpha[i] \
                        + sum(model.beta[i, j] * self.x[i][j] for j in model.J) \
                        + model.epsilon_plus[i] - model.epsilon_minus[i]
                return regression_rule

            elif self.rts == RTS_CRS:
                if type(self.z) != type(None):
                    def regression_rule(model, i):
                        return self.y[i] == sum(model.beta[i, j] * self.x[i][j] for j in model.J) \
                            - sum(model.lamda[k] * self.z[i][k] for k in model.K) \
                            + model.epsilon_plus[i] - model.epsilon_minus[i]
                    return regression_rule

                def regression_rule(model, i):
                    return self.y[i] == sum(model.beta[i, j] * self.x[i][j] for j in model.J) \
                        + model.epsilon_plus[i] - model.epsilon_minus[i]
                return regression_rule

        elif self.cet == CET_MULT:
            if type(self.z) != type(None):
                def regression_rule(model, i):
                    return log(self.y[i]) == log(model.frontier[i] + 1) \
                        - sum(model.lamda[k] * self.z[i][k] for k in model.K) \
                        + model.epsilon_plus[i] - model.epsilon_minus[i]
                return regression_rule

            def regression_rule(model, i):
                return log(self.y[i]) == log(model.frontier[i] + 1) \
                        + model.epsilon_plus[i] - model.epsilon_minus[i]
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
        """Return the proper afriat inequality constraint"""
        if self.fun == FUN_PROD:
            __operator = NumericValue.__le__
        elif self.fun == FUN_COST:
            __operator = NumericValue.__ge__

        if self.cet == CET_ADDI:
            if self.rts == RTS_VRS1:

                def afriat_rule(model, i, h):
                    if i == h:
                        return Constraint.Skip
                    return __operator(
                        model.alpha[i] + sum(model.beta[i, j] * self.x[i][j]
                                             for j in model.J),
                        model.alpha[h] + sum(model.beta[h, j] * self.x[i][j]
                                             for j in model.J))

                return afriat_rule
            elif self.rts == RTS_CRS:
                def afriat_rule(model, i, h):
                    if i == h:
                        return Constraint.Skip
                    return __operator(
                        sum(model.beta[i, j] * self.x[i][j]
                            for j in model.J),
                        sum(model.beta[h, j] * self.x[i][j]
                            for j in model.J))

                return afriat_rule
        elif self.cet == CET_MULT:
            if self.rts == RTS_VRS1:

                def afriat_rule(model, i, h):
                    if i == h:
                        return Constraint.Skip
                    return __operator(
                        model.alpha[i] + sum(model.beta[i, j] * self.x[i][j]
                                             for j in model.J),
                        model.alpha[h] + sum(model.beta[h, j] * self.x[i][j]
                                             for j in model.J))

                return afriat_rule
            elif self.rts == RTS_CRS:

                def afriat_rule(model, i, h):
                    if i == h:
                        return Constraint.Skip
                    return __operator(
                        sum(model.beta[i, j] * self.x[i][j] for j in model.J),
                        sum(model.beta[h, j] * self.x[i][j] for j in model.J))

                return afriat_rule

        raise ValueError("Undefined model parameters.")

    def display_status(self):
        """Display the status of problem"""
        print(self.optimization_status)

    def display_alpha(self):
        """Display alpha value"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_various_return_to_scale(self.rts)
        self.__model__.alpha.display()

    def display_beta(self):
        """Display beta value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.beta.display()

    def display_lamda(self):
        """Display lamda value"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_contextual_variable(self.z)
        self.__model__.lamda.display()

    # def display_residual(self):
    #     """Display residual value"""
    #     tools.assert_optimized(self.optimization_status)
    #     self.__model__.epsilon_plus.display()

    def display_positive_residual(self):
        """Display positive residual value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.epsilon_plus.display()

    def display_negative_residual(self):
        """Display negative residual value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.epsilon_minus.display()

    def get_status(self):
        """Return status"""
        return self.optimization_status

    def get_alpha(self):
        """Return alpha value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_various_return_to_scale(self.rts)
        alpha = pd.Series(self.__model__.alpha.extract_values(),name='alpha')
        return alpha

    def get_beta(self):
        """Return beta value by array"""
        tools.assert_optimized(self.optimization_status)
        beta = pd.Series(self.__model__.beta.extract_values(),index=self.__model__.beta.extract_values().keys())
        # if the series is multi-indexed we need to unstack it...
        if type(beta.index[0]) == tuple:  # it is multi-indexed
            beta = beta.unstack(level=1)
        else:
            beta = pd.DataFrame(beta)  # force transition from Series -> df
        # multi-index the columns
        beta.columns = map(lamda x: "beta"+str(x) ,beta.columns)
        return beta

    def get_lamda(self):
        """Return beta value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_contextual_variable(self.z)
        lamda = pd.DataFrame(self.__model__.lamda.extract_values(),index=self.z.index)
        lamda.columns = map(lamda x: "lamda"+str(x) ,lamda.columns)
        return lamda

    def get_residual(self):
        """Return residual value by array"""
        tools.assert_optimized(self.optimization_status)
        residual = pd.Series(self.get_positive_residual()- self.get_negative_residual(),name='epsilon')

        return residual

    def get_positive_residual(self):
        """Return positive residual value by array"""
        tools.assert_optimized(self.optimization_status)
        residual_plus = pd.Series(self.__model__.epsilon_plus.extract_values(),name='epsilon_plus')
        return residual_plus

    def get_negative_residual(self):
        """Return negative residual value by array"""
        tools.assert_optimized(self.optimization_status)
        epsilon_minus = pd.Series(self.__model__.epsilon_minus.extract_values(),name='epsilon_minus')
        return epsilon_minus

    # def get_frontier(self):
    #     """Return estimated frontier value by array"""
    #     tools.assert_optimized(self.optimization_status)
    #     if self.cet == CET_MULT and type(self.z) == type(None):
    #         frontier = np.asarray(list(self.__model__.frontier[:].value)) + 1
    #     elif self.cet == CET_MULT and type(self.z) != type(None):
    #         frontier = list(np.divide(np.exp(
    #              self.get_residual() + self.get_lamda() * np.asarray(self.z)[:, 0]),  self.b) - 1)
    #     elif self.cet == CET_ADDI:
    #         frontier = np.asarray(self.y) + self.get_residual()
    #     return np.asarray(frontier)

    # def get_predict(self, x_test):
    #     """Return the estimated function in testing sample"""
    #     tools.assert_optimized(self.optimization_status)
    #     return interpolation.interpolation(self.get_alpha(), self.get_beta(), x_test, fun=self.fun)

    # def get_delta(self):
    #     """Return delta value by array"""
    #     tools.assert_optimized(self.optimization_status)
    #     tools.assert_undesirable_output(self.b)
    #     delta = pd.Series(self.__model__.delta.extract_values(),index=self.__model__.delta.extract_values().keys())
    #     # if the series is multi-indexed we need to unstack it...
    #     if type(delta.index[0]) == tuple:  # it is multi-indexed
    #         delta = delta.unstack(level=1)
    #     else:
    #         delta = pd.DataFrame(delta)  # force transition from Series -> df
    #     # multi-index the columns
    #     delta.columns = map(lamda x: "beta"+str(x) ,delta.columns)
    #     return delta

class CER(CQR):
    """Convex expectile regression (CER)
    """

    def __init__(self, y, x, tau, z=None, cet=CET_ADDI, fun=FUN_PROD, rts=RTS_VRS1):
        """CER model

        Args:
            y (float): output variable.
            x (float): input variables.
            tau (float): expectile.
            z (float, optional): Contextual variable(s). Defaults to None.
            cet (String, optional): CET_ADDI (additive composite error term) or CET_MULT (multiplicative composite error term). Defaults to CET_ADDI.
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS1.
        """
        super().__init__(y, x, tau, z, cet, fun, rts)
        self.__model__.objective.deactivate()
        self.__model__.squared_objective = Objective(
            rule=self.__squared_objective_rule(), sense=minimize, doc='squared objective rule')

    def __squared_objective_rule(self):
        def squared_objective_rule(model):
            return (1 - self.tau) * sum(model.epsilon_plus[i] ** 2 for i in model.I) \
                + self.tau* sum(model.epsilon_minus[i] ** 2 for i in model.I)

        return squared_objective_rule
