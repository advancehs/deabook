# import dependencies
from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, Constraint, log
from pyomo.core.expr.numvalue import NumericValue
import numpy as np
import pandas as pd

from .constant import CET_ADDI, CET_MULT, FUN_PROD, FUN_COST, OPT_DEFAULT, RTS_CRS, RTS_VRS1, OPT_LOCAL
from .utils import tools
from .core import CNLSResultMixin

class CNLSSD(CNLSResultMixin):
    """Convex Nonparametric Least Square for Shephard distance function(CNLSSD)
    """

    def __init__(self, data, sent, z=None, gy=[1], gx=[0], cet=CET_MULT, fun=FUN_PROD, rts=RTS_VRS1):
        """CNLS model

        Args:
            # y (float): output variable.
            # x (float): input variables.
            data
            sent
            z (float, optional): Contextual variable(s). Defaults to None.
            gy (list, optional): output distance vector. Defaults to [1].
            gx (list, optional): input distance vector. Defaults to [0].
            cet (String, optional): CET_ADDI (additive composite error term) or CET_MULT (multiplicative composite error term). Defaults to CET_ADDI.
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS.
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.y, self.x, self.z, self.gy, self.gx, self.basexy= tools.assert_CNLSSD(data, sent, z, gy, gx) ## 注意，这里的数据已经除以x1，或者y1了。

        self.cet = cet
        self.fun = fun
        self.rts = rts

        # Initialize the CNLS model
        self.__model__ = ConcreteModel()

        if type(self.z) != type(None):
            # Initialize the set of z
            self.__model__.M = Set(initialize=range(len(self.z[0])))

            # Initialize the variables for z variable
            self.__model__.omega = Var(self.__model__.M, doc='z coefficient')

        # Initialize the sets
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))
        self.__model__.K = Set(initialize=range(len(self.y[0])))

        # Initialize the variables
        self.__model__.delta = Var(self.__model__.I,self.__model__.J,bounds=(0.0, None),doc='delta')
        self.__model__.gamma = Var(self.__model__.I,self.__model__.K,bounds=(0.0, None),doc='gamma')
        self.__model__.epsilon = Var(self.__model__.I, doc='residual')
        self.__model__.frontier = Var(self.__model__.I,
                                      bounds=( 0.0, None),
                                      doc='estimated frontier')
        if self.rts == RTS_VRS1:
            self.__model__.alpha = Var(self.__model__.I, doc='alpha')

        # Setup the objective function and constraints
        self.__model__.objective = Objective(rule=self.__objective_rule(),sense=minimize,
                                             doc='objective function')
        self.__model__.regression_rule = Constraint(self.__model__.I,rule=self.__regression_rule(),
                                                    doc='regression equation')
        if self.cet == CET_MULT:
            self.__model__.log_rule = Constraint(self.__model__.I, rule=self.__log_rule(),
                                                 doc='log-transformed regression equation')

        self.__model__.afriat_rule = Constraint(self.__model__.I,self.__model__.I, rule=self.__afriat_rule(),
                                                doc='afriat inequality')

        self.__model__.translation_rule = Constraint(self.__model__.I,
                                                     rule=self.__translation_property(),
                                                     doc='translation property')
        # self.__model__.log_rule.pprint()

        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Optimize the function by requested method

        Args:
            email (string): The email address for remote optimization. It will optimize locally if OPT_LOCAL is given.
            solver (string): The solver chosen for optimization. It will optimize with default solver if OPT_DEFAULT is given.
        """
        self._optimize_cnls_model(email=email, solver=solver, cet=self.cet)

    def __objective_rule(self):
        """Return the proper objective function"""

        def objective_rule(model):
            return sum(model.epsilon[i] ** 2 for i in model.I)

        return objective_rule

    def __regression_rule(self):
        """Return the proper regression constraint"""
        if self.cet == CET_MULT:
            if sum(self.gx) >= 1 and sum(self.gy) ==0:
                if type(self.z) != type(None):
                    def regression_rule(model, i):
                        return  - log(self.basexy[i]) == log(model.frontier[i] +1) \
                            + sum(model.omega[m] * self.z[i][m] for m in model.M) - model.epsilon[i]
                    return regression_rule

                def regression_rule(model, i):
                    return  - log(self.basexy[i]) == log(model.frontier[i]+1 ) - model.epsilon[i]
                return regression_rule
            elif sum(self.gy) >= 1 and sum(self.gx) ==0:
                if type(self.z) != type(None):
                    def regression_rule(model, i):
                        return  - log(self.basexy[i]) == log(model.frontier[i] +1) \
                            + sum(model.omega[m] * self.z[i][m] for m in model.M) + model.epsilon[i]
                    return regression_rule

                def regression_rule(model, i):
                    return  - log(self.basexy[i]) == log(model.frontier[i]+1 ) - model.epsilon[i]
                return regression_rule
        raise ValueError("Undefined model parameters.")

    def __log_rule(self):
        """Return the proper log constraint"""
        if self.cet == CET_MULT:
            if self.rts == RTS_VRS1:
                def log_rule(model, i):
                    return model.frontier[i] == model.alpha[i] +\
                        sum(model.delta[i, j] * self.x[i][j] for j in model.J) -\
                        sum(model.gamma[i, k] * self.y[i][k] for k in model.K)
                return log_rule

            elif self.rts == RTS_CRS:
                def log_rule(model, i):
                    return model.frontier[i] == \
                        sum(model.delta[i, j] * self.x[i][j] for j in model.J) -\
                        sum(model.gamma[i, k] * self.y[i][k] for k in model.K)
                return log_rule

        raise ValueError("Undefined model parameters.")

    def __afriat_rule(self):
        """Return the proper afriat inequality constraint"""
        if self.fun == FUN_PROD:
            __operator = NumericValue.__le__
        elif self.fun == FUN_COST:
            __operator = NumericValue.__ge__


        if self.cet == CET_MULT:
            if self.rts == RTS_VRS1:

                def afriat_rule(model, i, h):
                    if i == h:
                        return Constraint.Skip
                    return __operator(
                        (model.frontier[i]+1)  ,
                        model.alpha[h] + sum(model.delta[h, j] * self.x[i][j]for j in model.J) - \
                                                        sum(model.gamma[h, k] * self.y[i][k]for k in model.K)      )

                return afriat_rule
            elif self.rts == RTS_CRS:

                def afriat_rule(model, i, h):
                    if i == h:
                        return Constraint.Skip
                    return __operator(
                        (model.frontier[i]+1),
                        sum(model.delta[h, j] * self.x[i][j]for j in model.J) - \
                                                        sum(model.gamma[h, k] * self.y[i][k]for k in model.K)      )

                return afriat_rule

        raise ValueError("Undefined model parameters.")

    def __translation_property(self):
        """Return the proper translation property

        Following Chambers, Chung, Färe (1996) and Kuosmanen & Johnson (2017,
        EJOR 262(2):792-801), the translation property normalizes the
        coefficients against the direction vector only — not against data
        values.  See also the pyStoNED reference implementation
        (ds2010/pyStoNED, CNLSDDF.__translation_property).
        """
        def translation_rule(model, i):
            return sum(model.delta[i, j] * self.gx[j]  for j in model.J) \
                + sum(model.gamma[i, k] * self.gy[k]  for k in model.K) == 1

        return translation_rule











    def get_residual(self):
        """Return residual value by array"""
        tools.assert_optimized(self.optimization_status)
        residual = list(  self._single_model().epsilon[:].value)
        if sum(self.gx) >= 1 and sum(self.gy) ==0:
            residual = [+1*v for v in residual]
        elif sum(self.gy) >= 1 and sum(self.gx) ==0:
            residual = [1*v for v in residual]
        # print("aaasss#,",list(self._single_model().epsilon[:].value))
        return np.asarray(residual)


    def get_frontier(self):
        """Return estimated frontier value by array"""
        tools.assert_optimized(self.optimization_status)
        if self.cet == CET_MULT :
            frontier_values = [v.value for v in self._single_model().frontier[:]]
            frontier = [v + 1 for v in frontier_values]
            return frontier
        raise ValueError("CET_ADDI can not be employed")




    # def get_predict(self, x_test):
    #     """Return the estimated function in testing sample"""
    #     tools.assert_optimized(self.optimization_status)
    #     return interpolation.interpolation(self.get_alpha(), self.get_delta(), x_test, fun=self.fun)



class CNLSSDG(CNLSSD):
    """Sparse-G variant of :class:`CNLSSD`.

    This class keeps the original CNLSSD Shephard-distance model, residual
    decomposition, translation property and solver path unchanged, but replaces
    the full pairwise Afriat constraints with the Sweet-spot sparse constraint
    set used by the package's G variants.
    """

    def __init__(self, data, sent, z=None, gy=[1], gx=[0],
                 cet=CET_MULT, fun=FUN_PROD, rts=RTS_VRS1):
        from .utils import sweet
        from .utils._g_common import cg_pairs, replace_constraint

        super().__init__(data, sent, z=z, gy=gy, gx=gx,
                         cet=cet, fun=fun, rts=rts)

        self.cutactive = sweet.sweet(self.x)
        self.active = np.zeros((len(self.x), len(self.x)))
        self.active2 = np.zeros((len(self.x), len(self.x)))

        replace_constraint(
            self,
            'afriat_rule',
            cg_pairs(self.cutactive),
            self._CNLSSD__afriat_rule(),
            'active afriat inequality'
        )

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Optimize the sparse CNLSSD-G model."""
        self._optimize_cnls_model(email=email, solver=solver, cet=self.cet)


class CNLSDDF(CNLSSD):
    """Convex Nonparametric Least Square with directional distance function
    """

    def __init__(self, data, sent, z=None, gy=[1], gx=[1], fun=FUN_PROD, rts=RTS_VRS1):
        """CNLS DDF model

        Args:
            # y (float): output variable.
            # x (float): input variables.
            data
            sent
            z (float, optional): Contextual variable(s). Defaults to None.
            gy (list, optional): output directional vector. Defaults to [1].
            gx (list, optional): input directional vector. Defaults to [1].
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS.
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.y, self.x, self.z, self.gy, self.gx, self.basexy = tools.assert_CNLSDDF(data, sent, z, gy, gx)
        self.fun = fun
        self.rts = rts

        # self.data = data
        # self.sent = sent
        # self.z = z
        # self.gy = gy
        # self.gx = gx

        self.__model__ = ConcreteModel()

        # Initialize the sets
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))
        self.__model__.K = Set(initialize=range(len(self.y[0])))

        # Initialize the variables
        if self.rts == RTS_VRS1:
            self.__model__.alpha = Var(self.__model__.I, doc='alpha')
        self.__model__.delta = Var(
            self.__model__.I, self.__model__.J, bounds=(0.0, None), doc='delta')
        self.__model__.epsilon = Var(self.__model__.I, doc='residuals')
        self.__model__.gamma = Var(
            self.__model__.I, self.__model__.K, bounds=(0.0, None), doc='gamma')

        if type(self.z) != type(None):
            # Initialize the set of z
            self.__model__.M = Set(initialize=range(len(self.z[0])))
            # Initialize the variables for z variable
            self.__model__.omega = Var(self.__model__.M, doc='z coefficient')

        # Setup the objective function and constraints
        self.__model__.objective = Objective(rule=self._CNLSSD__objective_rule(),
                                             sense=minimize,
                                             doc='objective function')
        self.__model__.regression_rule = Constraint(self.__model__.I,
                                                    rule=self.__regression_rule(),
                                                    doc='regression equation')
        self.__model__.afriat_rule = Constraint(self.__model__.I,
                                                self.__model__.I,
                                                rule=self.__afriat_rule(),
                                                doc='afriat inequality')

        self.__model__.translation_rule = Constraint(self.__model__.I,
                                                     rule=self.__translation_property(),
                                                     doc='translation property')
        # self.__model__.translation_rule.pprint()


        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Optimize the function by requested method

        Args:
            email (string): The email address for remote optimization. It will optimize locally if OPT_LOCAL is given.
            solver (string): The solver chosen for optimization. It will optimize with default solver if OPT_DEFAULT is given.
        """
        self._optimize_cnls_model(email=email, solver=solver, cet=CET_ADDI)

    def __regression_rule(self):
        """Return the proper regression constraint"""
        if self.rts == RTS_VRS1:
            if type(self.z) != type(None):
                def regression_rule(model, i):
                    return model.epsilon[i] + self.basexy[i] \
                        == model.alpha[i] \
                        + sum(model.delta[i, j] * self.x[i][j] for j in model.J) \
                        - sum(model.gamma[i, k] * self.y[i][k] for k in model.K) \
                        - sum(model.omega[m] * self.z[i][m] for m in model.M)

                return regression_rule
            elif type(self.z) == type(None):
                def regression_rule(model, i):
                    return model.epsilon[i] + self.basexy[i] \
                        == model.alpha[i] \
                        + sum(model.delta[i, j] * self.x[i][j] for j in model.J) \
                        - sum(model.gamma[i, k] * self.y[i][k] for k in model.K)

                return regression_rule



        elif self.rts == RTS_CRS:
            if type(self.z) != type(None):
                def regression_rule(model, i):
                    return self.basexy[i] \
                        == sum(model.delta[i, j] * self.x[i][j] for j in model.J) \
                        - sum(model.gamma[i, k] * self.y[i][k] for k in model.K) \
                        - sum(model.omega[m] * self.z[i][m] for m in model.M)-model.epsilon[i]

                return regression_rule
            elif type(self.z) == type(None):
                def regression_rule(model, i):
                    return self.basexy[i] \
                        == sum(model.delta[i, j] * self.x[i][j] for j in model.J) \
                        - sum(model.gamma[i, k] * self.y[i][k] for k in model.K)-model.epsilon[i]

                return regression_rule

        raise ValueError("Undefined model parameters.")

    def __afriat_rule(self):
        """Return the proper afriat inequality constraint"""
        if self.fun == FUN_PROD:
            __operator = NumericValue.__le__
        elif self.fun == FUN_COST:
            __operator = NumericValue.__ge__

        if self.rts == RTS_VRS1:
            def afriat_rule(model, i, h):
                if i == h:
                    return Constraint.Skip
                return __operator(model.alpha[i]
                                  + sum(model.delta[i, j] * self.x[i][j] for j in model.J) \
                                  - sum(model.gamma[i, k] * self.y[i][k] for k in model.K),
                                  model.alpha[h]
                                  + sum(model.delta[h, j] * self.x[i][j] for j in model.J) \
                                  - sum(model.gamma[h, k] * self.y[i][k] for k in model.K))

            return afriat_rule



        elif self.rts == RTS_CRS:
            def afriat_rule(model, i, h):
                if i == h:
                    return Constraint.Skip
                return __operator(sum(model.delta[i, j] * self.x[i][j] for j in model.J) \
                                  - sum(model.gamma[i, k] * self.y[i][k] for k in model.K),
                                  sum(model.delta[h, j] * self.x[i][j] for j in model.J) \
                                  - sum(model.gamma[h, k] * self.y[i][k] for k in model.K))

            return afriat_rule

        raise ValueError("Undefined model parameters.")

    def __translation_property(self):
        """Return the proper translation property

        Standard form: coefficients × direction vector, independent of data
        (Chambers, Chung, Färe 1996; Kuosmanen & Johnson 2017, EJOR 262(2):792).
        Matches pyStoNED reference (ds2010/pyStoNED, CNLSDDF).
        """

        def translation_rule(model, i):
            return sum(model.delta[i, j] * self.gx[j]  for j in model.J) \
                + sum(model.gamma[i, k]  * self.gy[k]   for k in model.K) == 1

        return translation_rule







    def get_frontier(self):
        """Return estimated frontier value by array"""
        tools.assert_optimized(self.optimization_status)

        frontier_values = [v.value for v in self.__model__.epsilon[:]]
        frontier = [  i+v  for i,v in zip(self.basexy, frontier_values)]
        return frontier

    def get_residual(self):
        """Return residual value by array.

        Uses the same DDF sign convention as CNLSDDFweak: if an output
        direction gy is active, the regression residual is already on the
        production composite-error scale; otherwise input directions use the
        opposite sign.
        """
        tools.assert_optimized(self.optimization_status)
        residual = list(self._single_model().epsilon[:].value)
        if sum(self.gy) >= 1:
            residual = [+1 * v for v in residual]
        elif sum(self.gx) >= 1:
            residual = [-1 * v for v in residual]
        return np.asarray(residual)

# -----------------------------------------------------------------------------
# Constraint-generation counterpart
# -----------------------------------------------------------------------------

class CNLSDDFG(CNLSDDF):
    """CNLSDDF with delayed / sparse Afriat constraints.

    The class keeps the public DataFrame + ``sent`` interface of ``CNLSDDF``.
    It initializes with the standard CNLSDDF equations and replaces the full
    pairwise Afriat block by the Sweet-spot sparse constraint set. This gives the
    same computational-method interface as the ``G`` family without changing the
    original full-constraint ``CNLSDDF`` class.
    """

    def __init__(self, data, sent, z=None, gy=[1], gx=[1], fun=FUN_PROD, rts=RTS_VRS1):
        from .utils import sweet
        from .utils._g_common import cg_pairs, replace_constraint
        super().__init__(data, sent, z=z, gy=gy, gx=gx, fun=fun, rts=rts)
        self.cutactive = sweet.sweet(self.x)
        self.active = np.zeros((len(self.x), len(self.x)))
        self.active2 = np.zeros((len(self.x), len(self.x)))
        replace_constraint(self, 'afriat_rule', cg_pairs(self.cutactive), self._CNLSDDF__afriat_rule(), 'active afriat inequality')

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        self._optimize_cnls_model(email=email, solver=solver, cet=CET_ADDI)

