# import dependencies
from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, Constraint
from pyomo.core.expr.numvalue import NumericValue
from .constant import FUN_PROD, FUN_COST, RTS_VRS1,RTS_CRS, CET_ADDI,OPT_DEFAULT, OPT_LOCAL
from . import  CQER
from .utils import tools
import numpy as np
import pandas as pd

class CQRDDF(CQER.CQR):
    """Convex quantile regression with directional distance function
    """

    def __init__(self, y, x, tau, b=None, z=None, gy=[1], gx=[1], gb=None, fun=FUN_PROD,rts=RTS_VRS1):
        """CQR DDF

        Args:
            y (float): output variable.
            x (float): input variables.
            tau (float): quantile.
            b (float), optional): undesirable output variables. Defaults to None.
            z (float, optional): Contextual variable(s). Defaults to None.
            gy (list, optional): output directional vector. Defaults to [1].
            gx (list, optional): input directional vector. Defaults to [1].
            gb (list, optional): undesirable output directional vector. Defaults to None.
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS1.
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.y, self.x, self.b, self.z,self.gy, self.gx, self.gb = \
            tools.assert_valid_direciontal_data_with_z(y,x,b,z,gy,gx,gb)
        self.tau = tau

        self.fun = fun
        self.rts = rts

        self.__model__ = ConcreteModel()

        # Initialize the sets
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))
        self.__model__.K = Set(initialize=range(len(self.y[0])))

        # Initialize the variables
        self.__model__.alpha = Var(self.__model__.I, doc='alpha')
        self.__model__.beta = Var(
            self.__model__.I, self.__model__.J, bounds=(0.0, None), doc='beta')
        self.__model__.gamma = Var(
            self.__model__.I, self.__model__.K, bounds=(0.0, None), doc='gamma')

        self.__model__.epsilon_plus = Var(
            self.__model__.I, bounds=(0.0, None), doc='positive error term')
        self.__model__.epsilon_minus = Var(
            self.__model__.I, bounds=(0.0, None), doc='negative error term')

        if type(self.b) != type(None):
            self.__model__.L = Set(initialize=range(len(self.b[0])))
            self.__model__.delta = Var(
                self.__model__.I, self.__model__.L, bounds=(0.0, None), doc='delta')

        self.__model__.objective = Objective(rule=self._CQR__objective_rule(),
                                             sense=minimize,
                                             doc='objective function')


        self.__model__.regression_rule = Constraint(self.__model__.I,
                                                    rule=self.__regression_rule(),
                                                    doc='regression equation')

        self.__model__.translation_rule = Constraint(self.__model__.I,
                                                     rule=self.__translation_property(),
                                                     doc='translation property')

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
            self.__model__, email, CET_ADDI, solver)

    def __regression_rule(self):
        """Return the proper regression constraint"""
        if self.rts == RTS_VRS1:
            if type(self.b) == type(None):
                if type(self.z) != type(None):
                    def regression_rule(model, i):
                        return sum(model.gamma[i, k] * self.y[i][k] for k in model.K) \
                            == model.alpha[i] \
                            + sum(model.beta[i, j] * self.x[i][j] for j in model.J) \
                            - sum(model.lamda[m] * self.z[i][m] for m in model.M) \
                            - model.epsilon_minus[i] + model.epsilon_plus[i]
                    return regression_rule

                elif type(self.z) == type(None):
                    def regression_rule(model, i):
                        return sum(model.gamma[i, k] * self.y[i][k] for k in model.K) \
                            == model.alpha[i] \
                            + sum(model.beta[i, j] * self.x[i][j] for j in model.J) \
                            - model.epsilon_minus[i] + model.epsilon_plus[i]
                    return regression_rule

            elif type(self.b) != type(None):
                if type(self.z) != type(None):
                    def regression_rule(model, i):
                        return sum(model.gamma[i, k] * self.y[i][k] for k in model.K) \
                            == model.alpha[i] \
                            + sum(model.beta[i, j] * self.x[i][j] for j in model.J) \
                            + sum(model.delta[i, l] * self.b[i][l] for l in model.L) \
                            - sum(model.lamda[m] * self.z[i][m] for m in model.M) \
                            - model.epsilon_minus[i] + model.epsilon_plus[i]
                    return regression_rule

                elif type(self.z) == type(None):
                    def regression_rule(model, i):
                        return sum(model.gamma[i, k] * self.y[i][k] for k in model.K) \
                            == model.alpha[i] \
                            + sum(model.beta[i, j] * self.x[i][j] for j in model.J) \
                            + sum(model.delta[i, l] * self.b[i][l] for l in model.L) \
                            - model.epsilon_minus[i] + model.epsilon_plus[i]
                    return regression_rule

        elif self.rts == RTS_CRS:
            if type(self.b) == type(None):
                if type(self.z) != type(None):
                    def regression_rule(model, i):
                        return sum(model.gamma[i, k] * self.y[i][k] for k in model.K) \
                            == sum(model.beta[i, j] * self.x[i][j] for j in model.J) \
                            - sum(model.lamda[m] * self.z[i][m] for m in model.M) \
                            - model.epsilon_minus[i] + model.epsilon_plus[i]
                    return regression_rule

                elif type(self.z) == type(None):
                    def regression_rule(model, i):
                        return sum(model.gamma[i, k] * self.y[i][k] for k in model.K) \
                            == sum(model.beta[i, j] * self.x[i][j] for j in model.J) \
                            - model.epsilon_minus[i] + model.epsilon_plus[i]
                    return regression_rule

            elif type(self.b) != type(None):
                if type(self.z) != type(None):
                    def regression_rule(model, i):
                        return sum(model.gamma[i, k] * self.y[i][k] for k in model.K) \
                            == sum(model.beta[i, j] * self.x[i][j] for j in model.J) \
                            + sum(model.delta[i, l] * self.b[i][l] for l in model.L) \
                            - sum(model.lamda[m] * self.z[i][m] for m in model.M) \
                            - model.epsilon_minus[i] + model.epsilon_plus[i]
                    return regression_rule

                elif type(self.z) == type(None):
                    def regression_rule(model, i):
                        return sum(model.gamma[i, k] * self.y[i][k] for k in model.K) \
                            == sum(model.beta[i, j] * self.x[i][j] for j in model.J) \
                            + sum(model.delta[i, l] * self.b[i][l] for l in model.L) \
                            - model.epsilon_minus[i] + model.epsilon_plus[i]
                    return regression_rule

        raise ValueError("Undefined model parameters.")

    def __translation_property(self):
        """Return the proper translation property"""
        if type(self.b) == type(None):
            def translation_rule(model, i):
                return sum(model.beta[i, j] * self.gx[j] for j in model.J) \
                    + sum(model.gamma[i, k] * self.gy[k] for k in model.K) == 1

            return translation_rule

        elif type(self.b) != type(None):
            def translation_rule(model, i):
                return sum(model.beta[i, j] * self.gx[j] for j in model.J) \
                    + sum(model.gamma[i, k] * self.gy[k] for k in model.K) \
                    + sum(model.delta[i, l] * self.gb[l] for l in model.L) == 1

            return translation_rule

    def __afriat_rule(self):
        """Return the proper afriat inequality constraint"""
        if self.fun == FUN_PROD:
            __operator = NumericValue.__le__
        elif self.fun == FUN_COST:
            __operator = NumericValue.__ge__

        if type(self.b) == type(None):
            def afriat_rule(model, i, h):
                if i == h:
                    return Constraint.Skip
                return __operator(model.alpha[i]
                                  + sum(model.beta[i, j] * self.x[i][j]
                                        for j in model.J)
                                  - sum(model.gamma[i, k] * self.y[i][k]
                                        for k in model.K),
                                  model.alpha[h]
                                  + sum(model.beta[h, j] * self.x[i][j]
                                        for j in model.J)
                                  - sum(model.gamma[h, k] * self.y[i][k] for k in model.K))
            return afriat_rule

        def afriat_rule(model, i, h):
            if i == h:
                return Constraint.Skip
            return __operator(model.alpha[i]
                              + sum(model.beta[i, j] * self.x[i][j]
                                    for j in model.J)
                              + sum(model.delta[i, l] * self.b[i][l]
                                    for l in model.L)
                              - sum(model.gamma[i, k] * self.y[i][k]
                                    for k in model.K),
                              model.alpha[h]
                              + sum(model.beta[h, j] * self.x[i][j]
                                    for j in model.J)
                              + sum(model.delta[h, l] * self.b[i][l]
                                    for l in model.L)
                              - sum(model.gamma[h, k] * self.y[i][k] for k in model.K))
        return afriat_rule

    def display_gamma(self):
        """Display beta value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.gamma.display()

    def get_frontier(self):
        """Return estimated frontier value by array"""
        raise ValueError("DDF has no frontier.")

    def get_gamma(self):
        """Return gamma value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_desirable_output(self.y)
        gamma = pd.Series(self.__model__.gamma.extract_values(),index=self.__model__.gamma.extract_values().keys())
        # if the series is multi-indexed we need to unstack it...
        if type(gamma.index[0]) == tuple:  # it is multi-indexed
            gamma = gamma.unstack(level=1)
        else:
            gamma = pd.DataFrame(gamma)  # force transition from Series -> df
        # multi-index the columns
        gamma.columns = map(lambda x: "gamma"+str(x) ,gamma.columns)
        return gamma

    def display_delta(self):
        """Display beta value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.delta.display()

    def get_delta(self):
        """Return delta value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_undesirable_output(self.b)
        delta = pd.Series(self.__model__.delta.extract_values(),index=self.__model__.delta.extract_values().keys())
        # if the series is multi-indexed we need to unstack it...
        if type(delta.index[0]) == tuple:  # it is multi-indexed
            delta = delta.unstack(level=1)
        else:
            delta = pd.DataFrame(delta)  # force transition from Series -> df
        # multi-index the columns
        delta.columns = map(lambda x: "delta"+str(x) ,delta.columns)
        return delta



class CERDDF(CQRDDF):
    """Convex expectile regression with DDF formulation
    """

    def __init__(self, y, x, tau, b=None, z=None, gy=[1], gx=[1], gb=None, fun=FUN_PROD,rts=RTS_VRS1):
        """CER DDF

        y (float): output variable.
        x (float): input variables.
        tau (float): quantile.
        b (float), optional): undesirable output variables. Defaults to None.
        z (float, optional): Contextual variable(s). Defaults to None.
        gy (list, optional): output directional vector. Defaults to [1].
        gx (list, optional): input directional vector. Defaults to [1].
        gb (list, optional): undesirable output directional vector. Defaults to None.
        fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
        rts (String, optional): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS1.
        """

        super().__init__(y, x,tau, b,z, gy, gx, gb, fun, rts)
        self.__model__.objective.deactivate()
        self.__model__.squared_objective = Objective(
            rule=self.__squared_objective_rule(), sense=minimize, doc='squared objective rule')

    def __squared_objective_rule(self):
        def squared_objective_rule(model):
            return (1 - self.tau) * sum(model.epsilon_plus[i] ** 2 for i in model.I) \
                +  self.tau* sum(model.epsilon_minus[i] ** 2 for i in model.I)
        return squared_objective_rule
