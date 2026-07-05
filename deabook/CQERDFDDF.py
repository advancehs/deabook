# -*- coding: utf-8 -*-
"""Combined CQER/CER and CQERDDF/CERDDF models.

This module combines the base convex quantile/expectile regression models,
their directional distance function variants, and the corresponding G
computational variants from ``CQERG.py`` and ``CQERDDFG.py``. The original
source modules are left unchanged.
"""

# -*- coding: utf-8 -*-
# import dependencies
import numpy as np
import pandas as pd
from .constant import CET_ADDI, CET_MULT, FUN_PROD, FUN_COST, RTS_CRS, RTS_VRS1, OPT_LOCAL, OPT_DEFAULT
import time

from .constant import CET_ADDI, FUN_PROD, FUN_COST, OPT_DEFAULT, RTS_CRS, RTS_VRS,OPT_LOCAL


# ---------------------------------------------------------------------------
# CQER / CER and DDF base models
# ---------------------------------------------------------------------------
"""Combined CQER and CQERDDF models.

This module merges the convex quantile/expectile regression models from
``CQER.py`` and their directional distance function variants from
``CQERDDF.py`` into one file. The original source modules are left unchanged.
"""

# import dependencies
from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, Constraint, log
from pyomo.core.expr.numvalue import NumericValue
import numpy as np
import pandas as pd

from .constant import CET_ADDI, CET_MULT, FUN_PROD, FUN_COST, RTS_CRS, RTS_VRS1, OPT_LOCAL, OPT_DEFAULT
from .utils import tools

from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, Constraint
from .constant import FUN_PROD, FUN_COST, RTS_VRS1,RTS_CRS, CET_ADDI,OPT_DEFAULT, OPT_LOCAL


# ---------------------------------------------------------------------------
# CQER / CER models
# ---------------------------------------------------------------------------
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
            return self.tau * sum(model.epsilon_plus[i] for i in model.I) \
                + (1 - self.tau) * sum(model.epsilon_minus[i] for i in model.I)

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
    #     """Dispaly residual value"""
    #     tools.assert_optimized(self.optimization_status)
    #     self.__model__.epsilon_plus.display()

    def display_positive_residual(self):
        """Dispaly positive residual value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.epsilon_plus.display()

    def display_negative_residual(self):
        """Dispaly negative residual value"""
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
        beta.columns = map(lambda x: "beta"+str(x) ,beta.columns)
        return beta

    def get_lamda(self):
        """Return beta value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_contextual_variable(self.z)
        lamda = pd.DataFrame(self.__model__.lamda.extract_values(),index=self.z.index)
        lamda.columns = map(lambda x: "lamda"+str(x) ,lamda.columns)
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
    #     delta.columns = map(lambda x: "beta"+str(x) ,delta.columns)
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
            return self.tau * sum(model.epsilon_plus[i] ** 2 for i in model.I) \
                + (1 - self.tau) * sum(model.epsilon_minus[i] ** 2 for i in model.I)

        return squared_objective_rule

# -----------------------------------------------------------------------------
# Constraint-generation counterparts
# -----------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# CQERDDF / CERDDF models
# ---------------------------------------------------------------------------
class CQRDDF(CQR):
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
        if type(self.z) != type(None):
            self.__model__.M = Set(initialize=range(len(self.z[0])))
            self.__model__.lamda = Var(self.__model__.M, doc='z coefficient')

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
        """Return the Afriat convexity inequality."""
        if self.fun == FUN_PROD:
            __operator = NumericValue.__le__
        elif self.fun == FUN_COST:
            __operator = NumericValue.__ge__
        else:
            raise ValueError("Undefined model parameters.")

        if self.rts == RTS_VRS1:
            if type(self.b) == type(None):
                def afriat_rule(model, i, h):
                    if i == h:
                        return Constraint.Skip
                    return __operator(
                        model.alpha[i]
                        + sum(model.beta[i, j] * self.x[i][j] for j in model.J)
                        - sum(model.gamma[i, k] * self.y[i][k] for k in model.K),
                        model.alpha[h]
                        + sum(model.beta[h, j] * self.x[i][j] for j in model.J)
                        - sum(model.gamma[h, k] * self.y[i][k] for k in model.K),
                    )
                return afriat_rule

            def afriat_rule(model, i, h):
                if i == h:
                    return Constraint.Skip
                return __operator(
                    model.alpha[i]
                    + sum(model.beta[i, j] * self.x[i][j] for j in model.J)
                    + sum(model.delta[i, l] * self.b[i][l] for l in model.L)
                    - sum(model.gamma[i, k] * self.y[i][k] for k in model.K),
                    model.alpha[h]
                    + sum(model.beta[h, j] * self.x[i][j] for j in model.J)
                    + sum(model.delta[h, l] * self.b[i][l] for l in model.L)
                    - sum(model.gamma[h, k] * self.y[i][k] for k in model.K),
                )
            return afriat_rule

        elif self.rts == RTS_CRS:
            if type(self.b) == type(None):
                def afriat_rule(model, i, h):
                    if i == h:
                        return Constraint.Skip
                    return __operator(
                        sum(model.beta[i, j] * self.x[i][j] for j in model.J)
                        - sum(model.gamma[i, k] * self.y[i][k] for k in model.K),
                        sum(model.beta[h, j] * self.x[i][j] for j in model.J)
                        - sum(model.gamma[h, k] * self.y[i][k] for k in model.K),
                    )
                return afriat_rule

            def afriat_rule(model, i, h):
                if i == h:
                    return Constraint.Skip
                return __operator(
                    sum(model.beta[i, j] * self.x[i][j] for j in model.J)
                    + sum(model.delta[i, l] * self.b[i][l] for l in model.L)
                    - sum(model.gamma[i, k] * self.y[i][k] for k in model.K),
                    sum(model.beta[h, j] * self.x[i][j] for j in model.J)
                    + sum(model.delta[h, l] * self.b[i][l] for l in model.L)
                    - sum(model.gamma[h, k] * self.y[i][k] for k in model.K),
                )
            return afriat_rule

        raise ValueError("Undefined model parameters.")

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

    def get_positive_residual(self):
        """Return positive residual component by array."""
        tools.assert_optimized(self.optimization_status)
        return np.asarray(list(getattr(self, "__model__").epsilon_plus[:].value))

    def get_negative_residual(self):
        """Return negative residual component by array."""
        tools.assert_optimized(self.optimization_status)
        return np.asarray(list(getattr(self, "__model__").epsilon_minus[:].value))

    def get_residual(self):
        """Return quantile residual epsilon = epsilon_plus - epsilon_minus."""
        return self.get_positive_residual() - self.get_negative_residual()




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
            return self.tau * sum(model.epsilon_plus[i] ** 2 for i in model.I) \
                +  (1 - self.tau) * sum(model.epsilon_minus[i] ** 2 for i in model.I)
        return squared_objective_rule

# -----------------------------------------------------------------------------
# Constraint-generation counterparts
# -----------------------------------------------------------------------------


from .utils import CQERG1, CQERG2, CQERZG1, CQERZG2, sweet, interpolation

# ---------------------------------------------------------------------------
# G variants from CQERG.py
# ---------------------------------------------------------------------------
class CQRG:
    """Convex quantile regression (CQR) with Genetic algorithm
    """

    def __init__(self, y, x, tau, z=None, cet=CET_ADDI, fun=FUN_PROD, rts=RTS_VRS1):
        """CQRG model

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
        self.cutactive = sweet.sweet(x)
        self.y, self.x, self.z = tools.assert_valid_basic_data(y, x, z)
        self.tau = tau
        self.cet = cet
        self.fun = fun
        self.rts = rts

        # active (added) violated concavity constraint by iterative procedure
        self.active = np.zeros((len(x), len(x)))
        # violated concavity constraint
        self.active2 = np.zeros((len(x), len(x)))

        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Optimize the function by requested method"""
        # TODO(error/warning handling): Check problem status after optimization
        self.t0 = time.time()
        if type(self.z) != type(None):
            model1 = CQERZG1.CQRZG1(
                self.y, self.x, self.z, self.tau, self.cutactive, self.cet, self.fun, self.rts)
        else:
            model1 = CQERG1.CQRG1(
                self.y, self.x, self.tau, self.cutactive, self.cet, self.fun, self.rts)
        model1.optimize(email, solver)
        self.alpha = model1.get_alpha()
        self.beta = model1.get_beta()
        self.__model__ = model1.__model__

        self.count = 0
        while self.__convergence_test(self.alpha, self.beta) > 0.0001:
            if type(self.z) != type(None):
                model2 = CQERZG2.CQRZG2(
                    self.y, self.x, self.z, self.tau, self.cutactive,self.active, self.cet, self.fun, self.rts)
            else:
                model2 = CQERG2.CQRG2(
                    self.y, self.x, self.tau, self.cutactive, self.active, self.cet, self.fun, self.rts)
            model2.optimize(email, solver)
            self.alpha = model2.get_alpha()
            self.beta = model2.get_beta()
            # TODO: Replace print with log system
            # print("Genetic Algorithm Convergence : %8f" %
            #       (self.__convergence_test(self.alpha, self.beta)))
            self.__model__ = model2.__model__
            self.count += 1
        self.optimization_status = 1
        self.tt = time.time() - self.t0

    def __convergence_test(self, alpha, beta):
        alpha = np.asarray(alpha, dtype=float).reshape(-1)
        beta = np.asarray(beta, dtype=float)
        x = np.asarray(self.x)
        activetmp1 = 0.0

        # go into the loop
        for i in range(len(x)):
            activetmp = 0.0
            # go into the sub-loop and find the violated concavity constraints
            for j in range(len(x)):
                if self.cet == CET_ADDI:
                    if self.rts == RTS_VRS1:
                        if self.fun == FUN_PROD:
                            self.active2[i, j] = alpha[i] + np.sum(beta[i, :] * x[i, :]) - \
                                alpha[j] - np.sum(beta[j, :] * x[i, :])
                        elif self.fun == FUN_COST:
                            self.active2[i, j] = - alpha[i] - np.sum(beta[i, :] * x[i, :]) + \
                                alpha[j] + np.sum(beta[j, :] * x[i, :])
                    if self.rts == RTS_CRS:
                        if self.fun == FUN_PROD:
                            self.active2[i, j] = np.sum(beta[i, :] * x[i, :]) \
                                - np.sum(beta[j, :] * x[i, :])
                        elif self.fun == FUN_COST:
                            self.active2[i, j] = - np.sum(beta[i, :] * x[i, :]) \
                                + np.sum(beta[j, :] * x[i, :])
                if self.cet == CET_MULT:
                    if self.rts == RTS_VRS1:
                        if self.fun == FUN_PROD:
                            self.active2[i, j] = alpha[i] + np.sum(beta[i, :] * x[i, :]) - \
                                alpha[j] - np.sum(beta[j, :] * x[i, :])
                        elif self.fun == FUN_COST:
                            self.active2[i, j] = - alpha[i] - np.sum(beta[i, :] * x[i, :]) + \
                                alpha[j] + np.sum(beta[j, :] * x[i, :])
                    if self.rts == RTS_CRS:
                        if self.fun == FUN_PROD:
                            self.active2[i, j] = np.sum(beta[i, :] * x[i, :]) - \
                                np.sum(beta[j, :] * x[i, :])
                        elif self.fun == FUN_COST:
                            self.active2[i, j] = - np.sum(beta[i, :] * x[i, :]) + \
                                np.sum(beta[j, :] * x[i, :])
                if self.active2[i, j] > activetmp:
                    activetmp = self.active2[i, j]
            # find the maximal violated constraint in sub-loop and added into the active matrix
            for j in range(len(x)):
                if self.active2[i, j] >= activetmp and activetmp > 0:
                    self.active[i, j] = 1
            if activetmp > activetmp1:
                activetmp1 = activetmp
        return activetmp

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
    #     """Dispaly residual value"""
    #     tools.assert_optimized(self.optimization_status)
    #     self.__model__.epsilon.display()

    def display_positive_residual(self):
        """Dispaly positive residual value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.epsilon_plus.display()

    def display_negative_residual(self):
        """Dispaly negative residual value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.epsilon_minus.display()

    def get_status(self):
        """Return status"""
        return self.optimization_status

    def get_alpha(self):
        """Return alpha value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_various_return_to_scale(self.rts)
        alpha = list(self.__model__.alpha[:].value)
        return np.asarray(alpha)

    def get_beta(self):
        """Return beta value by array"""
        tools.assert_optimized(self.optimization_status)
        beta = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.beta),
                                                          list(self.__model__.beta[:, :].value))])
        beta = pd.DataFrame(beta, columns=['Name', 'Key', 'Value'])
        beta = beta.pivot(index='Name', columns='Key', values='Value')
        return beta.to_numpy()

    def get_residual(self):
        """Return residual value by array"""
        tools.assert_optimized(self.optimization_status)
        residual = np.asarray(list(self.__model__.epsilon_minus[:].value))\
                   - np.asarray(list(self.__model__.epsilon_plus[:].value))
        return residual

    def get_positive_residual(self):
        """Return positive residual value by array"""
        tools.assert_optimized(self.optimization_status)
        residual_plus = list(self.__model__.epsilon_plus[:].value)
        return np.asarray(residual_plus)

    def get_negative_residual(self):
        """Return negative residual value by array"""
        tools.assert_optimized(self.optimization_status)
        residual_minus = list(self.__model__.epsilon_minus[:].value)
        return np.asarray(residual_minus)

    def get_lamda(self):
        """Return beta value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_contextual_variable(self.z)
        lamda = list(self.__model__.lamda[:].value)
        return np.asarray(lamda)

    def get_frontier(self):
        """Return estimated frontier value by array"""
        tools.assert_optimized(self.optimization_status)
        if self.cet == CET_MULT and type(self.z) == type(None):
            frontier = np.asarray(list(self.__model__.frontier[:].value)) + 1
        elif self.cet == CET_MULT and type(self.z) != type(None):
            frontier = list(np.multiply(self.y, np.exp(
                self.get_residual() + self.get_lamda() * np.asarray(self.z)[:, 0])) - 1)
        elif self.cet == CET_ADDI:
            frontier = np.asarray(self.y) + self.get_residual()
        return np.asarray(frontier)

    # def get_totalconstr(self):
    #     """Return the number of total constraints"""
    #     tools.assert_optimized(self.optimization_status)
    #     activeconstr = np.sum(self.active) - np.trace(self.active)
    #     cutactiveconstr = np.sum(self.cutactive) - np.trace(self.cutactive)
    #     totalconstr = activeconstr + cutactiveconstr + 2 * len(self.active) + 1
    #     return totalconstr
    #
    # def get_runningtime(self):
    #     """Return the running time"""
    #     tools.assert_optimized(self.optimization_status)
    #     return self.tt
    #
    # def get_blocks(self):
    #     """Return the number of blocks"""
    #     tools.assert_optimized(self.optimization_status)
    #     return self.count
    #
    # def get_predict(self, x_test):
    #     """Return the estimated function in testing sample"""
    #     tools.assert_optimized(self.optimization_status)
    #     return interpolation.interpolation(self.get_alpha(), self.get_beta(), x_test, fun=self.fun)
    #

class CERG:
    """Convex expectile regression (CER) with Genetic algorithm
    """

    def __init__(self, y, x, tau, z=None, cet=CET_ADDI, fun=FUN_PROD, rts=RTS_VRS1):
        """CERG model

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
        self.cutactive = sweet.sweet(x)
        self.y, self.x, self.z = tools.assert_valid_basic_data(y, x, z)
        self.tau = tau
        self.cet = cet
        self.fun = fun
        self.rts = rts

        # active (added) violated concavity constraint by iterative procedure
        self.active = np.zeros((len(x), len(x)))
        # violated concavity constraint
        self.active2 = np.zeros((len(x), len(x)))

        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Optimize the function by requested method"""
        # TODO(error/warning handling): Check problem status after optimization
        self.t0 = time.time()
        if type(self.z) != type(None):
            model1 = CQERZG1.CERZG1(
                self.y, self.x, self.z, self.tau, self.cutactive, self.cet, self.fun, self.rts)
        else:
            model1 = CQERG1.CERG1(
                self.y, self.x, self.tau, self.cutactive, self.cet, self.fun, self.rts)
        model1.optimize(email, solver)
        self.alpha = model1.get_alpha()
        self.beta = model1.get_beta()
        self.__model__ = model1.__model__

        self.count = 0
        while self.__convergence_test(self.alpha, self.beta) > 0.0001:
            if type(self.z) != type(None):
                model2 = CQERZG2.CERZG2(
                    self.y, self.x, self.z, self.tau, self.cutactive, self.active, self.cet, self.fun, self.rts)
            else:
                model2 = CQERG2.CERG2(
                    self.y, self.x, self.tau, self.cutactive, self.active, self.cet, self.fun, self.rts)
            model2.optimize(email, solver)
            self.alpha = model2.get_alpha()
            self.beta = model2.get_beta()
            # TODO: Replace print with log system
            # print("Genetic Algorithm Convergence : %8f" %
            #       (self.__convergence_test(self.alpha, self.beta)))
            self.__model__ = model2.__model__
            self.count += 1
        self.optimization_status = 1
        self.tt = time.time() - self.t0

    def __convergence_test(self, alpha, beta):
        alpha = np.asarray(alpha, dtype=float).reshape(-1)
        beta = np.asarray(beta, dtype=float)
        x = np.asarray(self.x)
        activetmp1 = 0.0

        # go into the loop
        for i in range(len(x)):
            activetmp = 0.0
            # go into the sub-loop and find the violated concavity constraints
            for j in range(len(x)):
                if self.cet == CET_ADDI:
                    if self.rts == RTS_VRS1:
                        if self.fun == FUN_PROD:
                            self.active2[i, j] = alpha[i] + np.sum(beta[i, :] * x[i, :]) - \
                                alpha[j] - np.sum(beta[j, :] * x[i, :])
                        elif self.fun == FUN_COST:
                            self.active2[i, j] = - alpha[i] - np.sum(beta[i, :] * x[i, :]) + \
                                alpha[j] + np.sum(beta[j, :] * x[i, :])
                    if self.rts == RTS_CRS:
                        if self.fun == FUN_PROD:
                            self.active2[i, j] = np.sum(beta[i, :] * x[i, :]) \
                                - np.sum(beta[j, :] * x[i, :])
                        elif self.fun == FUN_COST:
                            self.active2[i, j] = - np.sum(beta[i, :] * x[i, :]) \
                                + np.sum(beta[j, :] * x[i, :])
                if self.cet == CET_MULT:
                    if self.rts == RTS_VRS1:
                        if self.fun == FUN_PROD:
                            self.active2[i, j] = alpha[i] + np.sum(beta[i, :] * x[i, :]) - \
                                alpha[j] - np.sum(beta[j, :] * x[i, :])
                        elif self.fun == FUN_COST:
                            self.active2[i, j] = - alpha[i] - np.sum(beta[i, :] * x[i, :]) + \
                                alpha[j] + np.sum(beta[j, :] * x[i, :])
                    if self.rts == RTS_CRS:
                        if self.fun == FUN_PROD:
                            self.active2[i, j] = np.sum(beta[i, :] * x[i, :]) - \
                                np.sum(beta[j, :] * x[i, :])
                        elif self.fun == FUN_COST:
                            self.active2[i, j] = - np.sum(beta[i, :] * x[i, :]) + \
                                np.sum(beta[j, :] * x[i, :])
                if self.active2[i, j] > activetmp:
                    activetmp = self.active2[i, j]
            # find the maximal violated constraint in sub-loop and added into the active matrix
            for j in range(len(x)):
                if self.active2[i, j] >= activetmp and activetmp > 0:
                    self.active[i, j] = 1
            if activetmp > activetmp1:
                activetmp1 = activetmp
        return activetmp

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
    #     """Dispaly residual value"""
    #     tools.assert_optimized(self.optimization_status)
    #     self.__model__.epsilon.display()

    def display_positive_residual(self):
        """Dispaly positive residual value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.epsilon_plus.display()

    def display_negative_residual(self):
        """Dispaly negative residual value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.epsilon_minus.display()

    def get_status(self):
        """Return status"""
        return self.optimization_status

    def get_alpha(self):
        """Return alpha value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_various_return_to_scale(self.rts)
        alpha = list(self.__model__.alpha[:].value)
        return np.asarray(alpha)

    def get_beta(self):
        """Return beta value by array"""
        tools.assert_optimized(self.optimization_status)
        beta = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.beta),
                                                          list(self.__model__.beta[:, :].value))])
        beta = pd.DataFrame(beta, columns=['Name', 'Key', 'Value'])
        beta = beta.pivot(index='Name', columns='Key', values='Value')
        return beta.to_numpy()

    def get_residual(self):
        """Return residual value by array"""
        tools.assert_optimized(self.optimization_status)
        residual = np.asarray(list(self.__model__.epsilon_minus[:].value))\
                   - np.asarray(list(self.__model__.epsilon_plus[:].value))
        return residual

    def get_positive_residual(self):
        """Return positive residual value by array"""
        tools.assert_optimized(self.optimization_status)
        residual_plus = list(self.__model__.epsilon_plus[:].value)
        return np.asarray(residual_plus)

    def get_negative_residual(self):
        """Return negative residual value by array"""
        tools.assert_optimized(self.optimization_status)
        residual_minus = list(self.__model__.epsilon_minus[:].value)
        return np.asarray(residual_minus)

    def get_lamda(self):
        """Return beta value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_contextual_variable(self.z)
        lamda = list(self.__model__.lamda[:].value)
        return np.asarray(lamda)

    def get_frontier(self):
        """Return estimated frontier value by array"""
        tools.assert_optimized(self.optimization_status)
        if self.cet == CET_MULT and type(self.z) == type(None):
            frontier = np.asarray(list(self.__model__.frontier[:].value)) + 1
        elif self.cet == CET_MULT and type(self.z) != type(None):
            frontier = list(np.multiply(self.y, np.exp(
                self.get_residual() + self.get_lamda() * np.asarray(self.z)[:, 0])) - 1)
        elif self.cet == CET_ADDI:
            frontier = np.asarray(self.y) + self.get_residual()
        return np.asarray(frontier)

    # def get_totalconstr(self):
    #     """Return the number of total constraints"""
    #     tools.assert_optimized(self.optimization_status)
    #     activeconstr = np.sum(self.active) - np.trace(self.active)
    #     cutactiveconstr = np.sum(self.cutactive) - np.trace(self.cutactive)
    #     totalconstr = activeconstr + cutactiveconstr + 2 * len(self.active) + 1
    #     return totalconstr
    #
    # def get_runningtime(self):
    #     """Return the running time"""
    #     tools.assert_optimized(self.optimization_status)
    #     return self.tt
    #
    # def get_blocks(self):
    #     """Return the number of blocks"""
    #     tools.assert_optimized(self.optimization_status)
    #     return self.count
    #
    # def get_predict(self, x_test):
    #     """Return the estimated function in testing sample"""
    #     tools.assert_optimized(self.optimization_status)
    #     return interpolation.interpolation(self.get_alpha(), self.get_beta(), x_test, fun=self.fun)
    #


from .utils import CQERDDFG1, CQERDDFG2, CQERDDFZG1, CQERDDFZG2

# ---------------------------------------------------------------------------
# G-DDF variants from CQERDDFG.py
# ---------------------------------------------------------------------------
class CQRDDFG:
    """Convex Nonparametric Least Square with weak disposability (weakCNLSDDF) and Genetic algorithm
    """
    def __init__(self, y, x, tau, b=None, z=None, gy=[1], gx=[1], gb=[1], fun=FUN_PROD, rts=RTS_VRS):
        """weakCNLSG model

        Args:
            y (ndarray): output variable.
            x (ndarray): input variables.
            tau (float): quantile.
            b (ndarray, optional): undersiable variables.
            z (ndarray, optional): Contextual variable(s). Defaults to None.
            gy (list, optional): output directional vector. Defaults to [1].
            gx (list, optional): input directional vector. Defaults to [1].
            gb (list, optional): undesirable output directional vector. Defaults to [1].
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS.
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.cutactive = sweet.sweet(x if type(b) == type(None) else np.column_stack((x,b)))
        self.y, self.x, self.b, self.z,self.gy, self.gx, self.gb = \
            tools.assert_valid_direciontal_data_with_z(y,x,b,z,gy,gx,gb)

        self.tau = tau

        self.fun = fun
        self.rts = rts

        # active (added) violated concavity constraint by iterative procedure
        self.active = np.zeros((len(x), len(x)))
        # violated concavity constraint
        self.active2 = np.zeros((len(x), len(x)))

        # active (added) violated concavity constraint for weak disposbility constrains by iterative procedure
        self.activeweak = np.zeros((len(x), len(x)))
        # violated concavity constraint for weak disposbility constrains
        self.activeweak2 = np.zeros((len(x), len(x)))

        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0



    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Optimize the function by requested method"""
        # TODO(error/warning handling): Check problem status after optimization
        self.t0 = time.time()
        if type(self.z) != type(None):
            model1 = CQERDDFZG1.CQRDDFZG1(
                self.y, self.x, self.b, self.z, self.tau, self.cutactive, self.gy, self.gx, self.gb, self.fun, self.rts)
        else:
            model1 = CQERDDFG1.CQRDDFG1(
                self.y, self.x, self.b, self.tau, self.cutactive, self.gy, self.gx, self.gb, self.fun, self.rts)
        model1.optimize(email, solver)
        self.alpha = model1.get_alpha()
        self.beta = model1.get_beta()
        self.gamma = model1.get_gamma()
        self.delta = model1.get_delta() if type(self.b) != type(None) else np.zeros((len(self.x), 0))
        self.__model__ = model1.__model__

        self.count = 0
        while self.__convergence_test(self.alpha, self.beta, self.gamma, self.delta) > 0.0001:
            if type(self.z) != type(None):
                model2 = CQERDDFZG2.CQRDDFZG2(
                    self.y, self.x, self.b, self.z, self.tau, self.cutactive, self.active,
                    self.gy, self.gx, self.gb,  self.fun, self.rts)
            else:
                model2 = CQERDDFG2.CQRDDFG2(
                    self.y, self.x, self.b, self.tau, self.cutactive, self.active,
                    self.gy, self.gx, self.gb, self.fun, self.rts)
            model2.optimize(email, solver)
            self.alpha = model2.get_alpha()
            self.beta = model2.get_beta()
            self.gamma = model2.get_gamma()
            self.delta = model2.get_delta() if type(self.b) != type(None) else np.zeros((len(self.x), 0))
            # TODO: Replace print with log system
            # print("Genetic Algorithm Convergence : %8f" %
            #       (self.__convergence_test(self.alpha, self.beta)))
            self.__model__ = model2.__model__
            self.count += 1
        self.optimization_status = 1
        self.tt = time.time() - self.t0

    def __convergence_test(self, alpha, beta, gamma, delta):
        alpha = np.asarray(alpha, dtype=float).reshape(-1)
        beta = np.asarray(beta, dtype=float)
        gamma = np.asarray(gamma, dtype=float)
        delta = np.asarray(delta, dtype=float)
        x = np.asarray(self.x)
        y = np.asarray(self.y)
        b = np.zeros((len(x), 0)) if type(self.b) == type(None) else np.asarray(self.b)
        activetmp1 = 0.0
        if self.rts == RTS_VRS and self.fun == FUN_PROD:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = alpha[i] + np.sum(beta[i, :] * x[i, :]) \
                        + np.sum(delta[i, :] * b[i, :]) - np.sum(gamma[i, :] * y[i, :]) \
                                         - alpha[j] - np.sum(beta[j, :] * x[i, :]) \
                        - np.sum(delta[j, :] * b[i, :] + np.sum(gamma[j, :] * y[i, :]))

                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp


        elif self.rts == RTS_VRS and self.fun == FUN_COST:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = -alpha[i] - np.sum(beta[i, :] * x[i, :]) \
                        - np.sum(delta[i, :] * b[i, :]) + np.sum(gamma[i, :] * y[i, :]) \
                                         + alpha[j] + np.sum(beta[j, :] * x[i, :]) \
                        + np.sum(delta[j, :] * b[i, :] - np.sum(gamma[j, :] * y[i, :]))

                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp


        elif self.rts == RTS_CRS and self.fun == FUN_PROD:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = np.sum(beta[i, :] * x[i, :]) \
                        + np.sum(delta[i, :] * b[i, :]) - np.sum(gamma[i, :] * y[i, :]) \
                                          - np.sum(beta[j, :] * x[i, :]) \
                        - np.sum(delta[j, :] * b[i, :] + np.sum(gamma[j, :] * y[i, :]))

                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_CRS and self.fun == FUN_COST:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = - np.sum(beta[i, :] * x[i, :]) \
                        - np.sum(delta[i, :] * b[i, :]) + np.sum(gamma[i, :] * y[i, :]) \
                                          + np.sum(beta[j, :] * x[i, :]) \
                        + np.sum(delta[j, :] * b[i, :] - np.sum(gamma[j, :] * y[i, :]))
                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp
        return activetmp


    def display_status(self):
        """Display the status of problem"""
        tools.assert_optimized(self.optimization_status)
        print(self.display_status)

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

    def display_residual(self):
        """Dispaly residual value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.epsilon.display()

    def display_delta(self):
        """Display delta value"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_undesirable_output(self.b)
        self.__model__.delta.display()

    def display_gamma(self):
        """Display gamma value"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_desirable_output(self.y)
        self.__model__.gamma.display()

    def get_status(self):
        """Return status"""
        return self.optimization_status

    def get_alpha(self):
        """Return alpha value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_various_return_to_scale(self.rts)
        alpha = list(self.__model__.alpha[:].value)
        return np.asarray(alpha)

    def get_beta(self):
        """Return beta value by array"""
        tools.assert_optimized(self.optimization_status)
        beta = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.beta),
                                                          list(self.__model__.beta[:, :].value))])
        beta = pd.DataFrame(beta, columns=['Name', 'Key', 'Value'])
        beta = beta.pivot(index='Name', columns='Key', values='Value')
        return beta.to_numpy()

    def get_residual(self):
        """Return residual value by array"""
        tools.assert_optimized(self.optimization_status)
        residual = np.asarray(list(self.__model__.epsilon_minus[:].value))\
                   - np.asarray(list(self.__model__.epsilon_plus[:].value))
        return residual

    def get_lamda(self):
        """Return beta value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_contextual_variable(self.z)
        lamda = list(self.__model__.lamda[:].value)
        return np.asarray(lamda)

    def get_frontier(self):
        """Return estimated frontier value by array"""
        raise ValueError("DDF hsa no frontier.")

    def get_delta(self):
        """Return delta value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_undesirable_output(self.b)
        delta = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.delta),
                                                           list(self.__model__.delta[:, :].value))])
        delta = pd.DataFrame(delta, columns=['Name', 'Key', 'Value'])
        delta = delta.pivot(index='Name', columns='Key', values='Value')
        return delta.to_numpy()

    def get_gamma(self):
        """Return gamma value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_desirable_output(self.y)
        gamma = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.gamma),
                                                           list(self.__model__.gamma[:, :].value))])
        gamma = pd.DataFrame(gamma, columns=['Name', 'Key', 'Value'])
        gamma = gamma.pivot(index='Name', columns='Key', values='Value')
        return gamma.to_numpy()



class CERDDFG(CQRDDFG):
    """Convex Nonparametric Least Square with weak disposability (weakCNLSDDF) and Genetic algorithm
    """
    def __init__(self, y, x, tau, b=None, z=None, gy=[1], gx=[1], gb=[1], fun=FUN_PROD, rts=RTS_VRS):
        """weakCNLSG model

        Args:
            y (ndarray): output variable.
            x (ndarray): input variables.
            tau (float): quantile.
            b (ndarray, optional): undersiable variables.
            z (ndarray, optional): Contextual variable(s). Defaults to None.
            gy (list, optional): output directional vector. Defaults to [1].
            gx (list, optional): input directional vector. Defaults to [1].
            gb (list, optional): undesirable output directional vector. Defaults to [1].
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS.
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.cutactive = sweet.sweet(x if type(b) == type(None) else np.column_stack((x,b)))
        self.y, self.x, self.b, self.z,self.gy, self.gx, self.gb = \
            tools.assert_valid_direciontal_data_with_z(y,x,b,z,gy,gx,gb)

        self.tau = tau

        self.fun = fun
        self.rts = rts

        # active (added) violated concavity constraint by iterative procedure
        self.active = np.zeros((len(x), len(x)))
        # violated concavity constraint
        self.active2 = np.zeros((len(x), len(x)))

        # active (added) violated concavity constraint for weak disposbility constrains by iterative procedure
        self.activeweak = np.zeros((len(x), len(x)))
        # violated concavity constraint for weak disposbility constrains
        self.activeweak2 = np.zeros((len(x), len(x)))

        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0



    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Optimize the function by requested method"""
        # TODO(error/warning handling): Check problem status after optimization
        self.t0 = time.time()
        if type(self.z) != type(None):
            model1 = CQERDDFZG1.CERDDFZG1(
                self.y, self.x, self.b, self.z, self.tau, self.cutactive, self.gy, self.gx, self.gb, self.fun, self.rts)
        else:
            model1 = CQERDDFG1.CERDDFG1(
                self.y, self.x, self.b, self.tau, self.cutactive, self.gy, self.gx, self.gb, self.fun, self.rts)
        model1.optimize(email, solver)
        self.alpha = model1.get_alpha()
        self.beta = model1.get_beta()
        self.gamma = model1.get_gamma()
        self.delta = model1.get_delta() if type(self.b) != type(None) else np.zeros((len(self.x), 0))
        self.__model__ = model1.__model__

        self.count = 0
        while self.__convergence_test(self.alpha, self.beta, self.gamma, self.delta) > 0.0001:
            if type(self.z) != type(None):
                model2 = CQERDDFZG2.CERDDFZG2(
                    self.y, self.x, self.b, self.z, self.tau, self.cutactive, self.active,
                    self.gy, self.gx, self.gb,  self.fun, self.rts)
            else:
                model2 = CQERDDFG2.CERDDFG2(
                    self.y, self.x, self.b, self.tau, self.cutactive, self.active,
                    self.gy, self.gx, self.gb, self.fun, self.rts)
            model2.optimize(email, solver)
            self.alpha = model2.get_alpha()
            self.beta = model2.get_beta()
            self.gamma = model2.get_gamma()
            self.delta = model2.get_delta() if type(self.b) != type(None) else np.zeros((len(self.x), 0))
            # TODO: Replace print with log system
            # print("Genetic Algorithm Convergence : %8f" %
            #       (self.__convergence_test(self.alpha, self.beta)))
            self.__model__ = model2.__model__
            self.count += 1
        self.optimization_status = 1
        self.tt = time.time() - self.t0

    def __convergence_test(self, alpha, beta, gamma, delta):
        alpha = np.asarray(alpha, dtype=float).reshape(-1)
        beta = np.asarray(beta, dtype=float)
        gamma = np.asarray(gamma, dtype=float)
        delta = np.asarray(delta, dtype=float)
        x = np.asarray(self.x)
        y = np.asarray(self.y)
        b = np.zeros((len(x), 0)) if type(self.b) == type(None) else np.asarray(self.b)
        activetmp1 = 0.0
        if self.rts == RTS_VRS and self.fun == FUN_PROD:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = alpha[i] + np.sum(beta[i, :] * x[i, :]) \
                        + np.sum(delta[i, :] * b[i, :]) - np.sum(gamma[i, :] * y[i, :]) \
                                         - alpha[j] - np.sum(beta[j, :] * x[i, :]) \
                        - np.sum(delta[j, :] * b[i, :] + np.sum(gamma[j, :] * y[i, :]))

                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp


        elif self.rts == RTS_VRS and self.fun == FUN_COST:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = -alpha[i] - np.sum(beta[i, :] * x[i, :]) \
                        - np.sum(delta[i, :] * b[i, :]) + np.sum(gamma[i, :] * y[i, :]) \
                                         + alpha[j] + np.sum(beta[j, :] * x[i, :]) \
                        + np.sum(delta[j, :] * b[i, :] - np.sum(gamma[j, :] * y[i, :]))

                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp


        elif self.rts == RTS_CRS and self.fun == FUN_PROD:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = np.sum(beta[i, :] * x[i, :]) \
                        + np.sum(delta[i, :] * b[i, :]) - np.sum(gamma[i, :] * y[i, :]) \
                                          - np.sum(beta[j, :] * x[i, :]) \
                        - np.sum(delta[j, :] * b[i, :] + np.sum(gamma[j, :] * y[i, :]))

                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_CRS and self.fun == FUN_COST:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = - np.sum(beta[i, :] * x[i, :]) \
                        - np.sum(delta[i, :] * b[i, :]) + np.sum(gamma[i, :] * y[i, :]) \
                                          + np.sum(beta[j, :] * x[i, :]) \
                        + np.sum(delta[j, :] * b[i, :] - np.sum(gamma[j, :] * y[i, :]))
                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp
        return activetmp

# ---------------------------------------------------------------------------
# Deterministic efficiency helpers for CQER/CER models
# ---------------------------------------------------------------------------
def _cqer_as_2d(value):
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 0:
        return arr.reshape(1, 1)
    if arr.ndim == 1:
        return arr.reshape(-1, 1)
    return arr


def _cqer_as_1d(value):
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 2 and arr.shape[1] == 1:
        return arr[:, 0]
    return arr.reshape(-1)


def _cqer_to_numpy(value):
    if hasattr(value, "to_numpy"):
        return np.asarray(value.to_numpy(), dtype=float)
    return np.asarray(value, dtype=float)


def _cqer_alpha(self):
    n = len(_cqer_as_2d(self.y))
    try:
        alpha = _cqer_to_numpy(self.get_alpha()).reshape(-1)
        if alpha.size == n:
            return alpha
    except Exception:
        pass
    return np.zeros(n, dtype=float)


def _cqer_lamda_effect(self):
    if getattr(self, "z", None) is None:
        return None
    try:
        z = _cqer_as_2d(self.z)
        lamda = _cqer_to_numpy(self.get_lamda()).reshape(-1)
        return z @ lamda
    except Exception:
        return None


def _cqer_linear_frontier(self):
    tools.assert_optimized(self.optimization_status)
    x = _cqer_as_2d(self.x)
    beta = _cqer_as_2d(_cqer_to_numpy(self.get_beta()))
    value = _cqer_alpha(self) + np.sum(beta * x, axis=1)
    if getattr(self, "b", None) is not None and hasattr(self, "get_delta"):
        try:
            b = _cqer_as_2d(self.b)
            delta = _cqer_as_2d(_cqer_to_numpy(self.get_delta()))
            if delta.shape[1] == b.shape[1]:
                value = value + np.sum(delta * b, axis=1)
        except Exception:
            pass
    z_effect = _cqer_lamda_effect(self)
    if z_effect is not None:
        if getattr(self, "cet", CET_ADDI) == CET_MULT:
            value = value * np.exp(-z_effect)
        else:
            value = value - z_effect
    return np.asarray(value, dtype=float)


def _cqer_get_frontier(self):
    """Return the fitted deterministic CQR/CER frontier at observed DMUs."""
    return _cqer_linear_frontier(self)


def _cqer_ratio(numerator, denominator):
    numerator = np.asarray(numerator, dtype=float)
    denominator = np.asarray(denominator, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.divide(
            numerator,
            denominator,
            out=np.full_like(numerator, np.nan, dtype=float),
            where=np.isfinite(numerator) & np.isfinite(denominator) & (denominator > 1e-12),
        )
    return np.clip(ratio, 0.0, 1.0)


def _cqer_get_technical_efficiency(self):
    """Return deterministic front-to-observation efficiency for CQR/CER.

    This is not a StoNED/JLMS efficiency: no RED_MOM, RED_QLE, or RED_KDE
    residual decomposition is applied.  For production frontiers it reports
    y / fitted_frontier; for cost frontiers it reports fitted_frontier / y.
    """
    y = _cqer_as_2d(self.y)[:, 0]
    frontier = _cqer_as_1d(self.get_frontier())
    if self.fun == FUN_PROD:
        te = _cqer_ratio(y, frontier)
    elif self.fun == FUN_COST:
        te = _cqer_ratio(frontier, y)
    else:
        raise ValueError("Undefined model parameters.")
    self.TE = te
    return te


def _cqer_get_technical_efficiency_components(self):
    """Return deterministic CQR/CER frontier efficiency components."""
    y = _cqer_as_2d(self.y)[:, 0]
    frontier = _cqer_as_1d(self.get_frontier())
    te = self.get_technical_efficiency()
    return pd.DataFrame({
        "y": y,
        "frontier": frontier,
        "residual": _cqer_as_1d(self.get_residual()),
        "TE": te,
    })


def _cqer_ddf_raw_distance(self):
    tools.assert_optimized(self.optimization_status)
    y = _cqer_as_2d(self.y)
    gamma = _cqer_as_2d(_cqer_to_numpy(self.get_gamma()))
    return _cqer_linear_frontier(self) - np.sum(gamma * y, axis=1)


def _cqer_ddf_get_directional_distance(self):
    """Return nonnegative deterministic DDF distance to the fitted CQR/CER frontier."""
    raw = _cqer_ddf_raw_distance(self)
    distance = np.where(np.isfinite(raw), np.maximum(raw, 0.0), np.nan)
    self.directional_distance = distance
    return distance


def _cqer_ddf_components(self):
    distance = _cqer_as_1d(self.get_directional_distance())
    x = _cqer_as_2d(self.x)
    y = _cqer_as_2d(self.y)
    b = None if getattr(self, "b", None) is None else _cqer_as_2d(self.b)
    gx = np.asarray(getattr(self, "gx", []), dtype=float).reshape(-1)
    gy = np.asarray(getattr(self, "gy", []), dtype=float).reshape(-1)
    gb = np.asarray(getattr(self, "gb", [] if b is None else np.ones(b.shape[1])), dtype=float).reshape(-1)
    columns = {
        "raw_directional_distance": _cqer_ddf_raw_distance(self),
        "directional_distance": distance,
    }
    scores = []
    for j, g in enumerate(gx[:x.shape[1]]):
        if g > 0:
            projected = x[:, j] - distance * g
            te = _cqer_ratio(projected, x[:, j])
            columns[f"inefficiency_x_{j}"] = distance * g
            columns[f"TE_x_{j}"] = te
            scores.append(te)
    for k, g in enumerate(gy[:y.shape[1]]):
        if g > 0:
            projected = y[:, k] + distance * g
            te = _cqer_ratio(y[:, k], projected)
            columns[f"inefficiency_y_{k}"] = distance * g
            columns[f"TE_y_{k}"] = te
            scores.append(te)
    if b is not None:
        for l, g in enumerate(gb[:b.shape[1]]):
            if g > 0:
                projected = b[:, l] - distance * g
                te = _cqer_ratio(projected, b[:, l])
                columns[f"inefficiency_b_{l}"] = distance * g
                columns[f"TE_b_{l}"] = te
                scores.append(te)
    if scores:
        columns["TE"] = np.nanmin(np.column_stack(scores), axis=1)
    else:
        columns["TE"] = np.full_like(distance, np.nan, dtype=float)
    return pd.DataFrame(columns)


def _cqer_ddf_get_technical_efficiency_components(self):
    """Return deterministic DDF efficiency components for CQRDDF/CERDDF."""
    comp = _cqer_ddf_components(self)
    self.TE = comp["TE"].to_numpy(dtype=float)
    return comp


def _cqer_ddf_get_technical_efficiency(self):
    """Return binding deterministic DDF efficiency for CQRDDF/CERDDF."""
    return self.get_technical_efficiency_components()["TE"].to_numpy(dtype=float)


for _cls in [CQR, CER, CQRG, CERG]:
    _cls.get_frontier = _cqer_get_frontier
    _cls.get_technical_efficiency = _cqer_get_technical_efficiency
    _cls.get_technical_efficiency_components = _cqer_get_technical_efficiency_components

for _cls in [CQRDDF, CERDDF, CQRDDFG, CERDDFG]:
    _cls.get_directional_distance = _cqer_ddf_get_directional_distance
    _cls.get_technical_efficiency = _cqer_ddf_get_technical_efficiency
    _cls.get_technical_efficiency_components = _cqer_ddf_get_technical_efficiency_components

