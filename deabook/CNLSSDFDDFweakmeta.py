# import dependencies
from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, Constraint, log
from pyomo.core.expr.numvalue import NumericValue
from .constant import CET_ADDI,CET_MULT, FUN_PROD,FUN_COST, RTS_CRS, RTS_VRS1, RED_MOM,RED_QLE,RED_KDE

import pandas as pd
import numpy as np

from . import StoNED
from .CNLSSDFDDFweak import CNLSSDweak,CNLSDDFweak
from .constant import CET_ADDI, CET_MULT, FUN_COST, FUN_PROD, RTS_VRS1, RTS_VRS2, RTS_CRS, OPT_DEFAULT, OPT_LOCAL
from .utils import tools
class CNLSSDweakmeta(CNLSSDweak):
    """Convex Nonparametric Least Square for Shephard distance function(CNLSSD)
    """

    def __init__(self, data, sent, z=None, gy=[1], gx=[0], gb=[0], cet=CET_MULT, fun=FUN_PROD, rts=RTS_VRS1):
        """CNLS model

        Args:
            # y (float): output variable.
            # x (float): input variables.
            # b (float): undesirable output variables.
            data
            sent
            z (float, optional): Contextual variable(s). Defaults to None.
            gy (list, optional): output distance vector. Defaults to [1].
            gx (list, optional): input distance vector. Defaults to [0].
            gb (list, optional): undesirable output directional vector. Defaults to [0].
            cet (String, optional): CET_ADDI (additive composite error term) or CET_MULT (multiplicative composite error term). Defaults to CET_ADDI.
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS.
        """
        # TODO(error/warning handling): Check the configuration of the model exist

        model = CNLSSDweak(data, sent=sent, z=z, gy=gy, gx=gx, gb=gb, cet=cet, fun=fun, rts=rts)
        model.optimize(solver="ipopt")
        rd = StoNED.StoNED(model)
        self.gce = rd.get_technical_efficiency(RED_QLE)
        self.y, self.x, self.b, self.z, self.gy, self.gx, self.gb, self.basexyb \
            = tools.assert_CNLSSDFweak(data, sent, z, gy, gx, gb) ## 注意，这里的数据已经除以x1，或者y1了。

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
        self.__model__.L = Set(initialize=range(len(self.b[0])))

        # Initialize the variables
        self.__model__.delta = Var(self.__model__.I,self.__model__.J,bounds=(0.0, None),doc='delta')
        self.__model__.gamma = Var(self.__model__.I,self.__model__.K,bounds=(0.0, None),doc='gamma')
        self.__model__.kappa = Var(self.__model__.I,self.__model__.K,bounds=(0.0, None),doc='kappa')

        self.__model__.epsilon = Var(self.__model__.I, doc='residual')
        self.__model__.frontier = Var(self.__model__.I,
                                      bounds=(0.0, None),
                                      doc='estimated frontier')
        if self.rts == RTS_VRS1:
            self.__model__.alpha = Var(self.__model__.I, bounds=(0.0, None),doc='alpha')
        elif self.rts == RTS_VRS2:
            self.__model__.alpha = Var(self.__model__.I, bounds=(None, None),doc='alpha')
        # Setup the objective function and constraints
        self.__model__.objective = Objective(rule=self.__objective_rule(),sense=minimize,
                                             doc='objective function')
        self.__model__.regression_rule = Constraint(self.__model__.I,rule=self.__regression_rule(),
                                                    doc='regression equation')
        if self.cet != CET_MULT:
            raise ValueError("cet must be CET_MULT.")

        self.__model__.log_rule = Constraint(self.__model__.I, rule=self.__log_rule(),
                                                 doc='log-transformed regression equation')

        self.__model__.afriat_rule = Constraint(self.__model__.I,self.__model__.I, rule=self.__afriat_rule(),
                                                doc='afriat inequality')

        if self.rts == RTS_VRS2:
            self.__model__.red_factor_rule = Constraint(self.__model__.I, self.__model__.I,
                                                    rule=self.__reduction_factor_rule(),
                                                    doc='reduction factor inequality')

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
            if type(self.z) != type(None):
                def regression_rule(model, i):
                    return model.epsilon[i] - log(self.gce[i])- log(self.basexyb[i]) \
                        == log(model.frontier[i] + 1) \
                        + sum(model.omega[m] * self.z[i][m] for m in model.M)
                return regression_rule

            def regression_rule(model, i):
                # print(log(self.basexy[i]),"################")
                return model.epsilon[i] - log(self.gce[i])- log(self.basexyb[i]) == log(model.frontier[i] + 1)
            return regression_rule

        raise ValueError("Undefined model parameters.")

    def __log_rule(self):
        """Return the proper log constraint"""
        if self.cet == CET_MULT:
            if self.rts == RTS_VRS1:
                def log_rule(model, i):
                    return (model.frontier[i] == - 1 + \
                        sum(model.delta[i, j] * self.x[i][j] for j in model.J) -\
                        sum(model.gamma[i, k] * self.y[i][k] for k in model.K) +\
                        sum(model.kappa[i, l] * self.b[i][l] for l in model.L))

                return log_rule

            elif self.rts == RTS_VRS2:
                def log_rule(model, i):
                    return (model.frontier[i] == - 1 + model.alpha[i] + \
                        sum(model.delta[i, j] * self.x[i][j] for j in model.J) -\
                        sum(model.gamma[i, k] * self.y[i][k] for k in model.K) +\
                        sum(model.kappa[i, l] * self.b[i][l] for l in model.L) )

                return log_rule

            elif self.rts == RTS_CRS:
                def log_rule(model, i):
                    return model.frontier[i] == - 1 +\
                        sum(model.delta[i, j] * self.x[i][j] for j in model.J) -\
                        sum(model.gamma[i, k] * self.y[i][k] for k in model.K) +\
                        sum(model.kappa[i, l] * self.b[i][l] for l in model.L)

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
                        - model.alpha[h] + sum(model.delta[h, j] * self.x[i][j]for j in model.J) - \
                                        sum(model.gamma[h, k] * self.y[i][k]for k in model.K) + \
                                        sum(model.kappa[h, l] * self.b[i][l]for l in model.L) )

                return afriat_rule

            elif self.rts == RTS_VRS2:

                def afriat_rule(model, i, h):
                    if i == h:
                        return Constraint.Skip
                    return __operator(
                        (model.frontier[i]+1)  ,
                        + model.alpha[h] + sum(model.delta[h, j] * self.x[i][j]for j in model.J) - \
                                        sum(model.gamma[h, k] * self.y[i][k]for k in model.K) + \
                                        sum(model.kappa[h, l] * self.b[i][l]for l in model.L) )

                return afriat_rule

            elif self.rts == RTS_CRS:

                def afriat_rule(model, i, h):
                    if i == h:
                        return Constraint.Skip
                    return __operator(
                        (model.frontier[i]+1),
                        sum(model.delta[h, j] * self.x[i][j]for j in model.J) - \
                                        sum(model.gamma[h, k] * self.y[i][k]for k in model.K) + \
                                        sum(model.kappa[h, l] * self.b[i][l]for l in model.L) )
                return afriat_rule

        raise ValueError("Undefined model parameters.")

    def __reduction_factor_rule(self):
        """Return the proper reduction_factor inequality constraint"""
        def reduction_factor_rule(model, i, h):
            return model.alpha[h] + sum(model.delta[h, j] * self.x[i][j] for j in model.J) >= 0

        return reduction_factor_rule
    def __translation_property(self):
        """Return the proper translation property"""
        def translation_rule(model, i):
            return sum(model.delta[i, j] * self.gx[j] for j in model.J) \
                + sum(model.gamma[i, k] * self.gy[k] for k in model.K) \
                + sum(model.kappa[i, l] * self.gb[l] for l in model.L) == 1

        return translation_rule












    # def get_residual(self):
    #     """Return residual value by array"""
    #     tools.assert_optimized(self.optimization_status)
    #     residual = list(self._single_model().epsilon[:].value)
    #     return np.asarray(residual)

    def get_residual(self):
        """Return residual value by array"""
        tools.assert_optimized(self.optimization_status)
        residual = list(  self._single_model().epsilon[:].value)
        if sum(self.gx) >= 1 and sum(self.gy) ==0 and sum(self.gb) == 0:
            residual = [+1*v for v in residual]
        elif sum(self.gb) >= 1 and sum(self.gx) ==0 and sum(self.gy) == 0:
            residual = [+1*v for v in residual]
        elif sum(self.gy) >= 1 and sum(self.gx) ==0 and sum(self.gb) == 0:
            residual = [-1*v for v in residual]
        # print("aaasss#,",list(self._single_model().epsilon[:].value))
        return np.asarray(residual)










class CNLSDDFweakmeta(CNLSDDFweak):
    """Convex Nonparametric Least Square with directional distance function
    """

    def __init__(self, data, sent, z=None, gy=[1], gx=[1], gb=[1], fun=FUN_PROD, rts=RTS_VRS1, kind = 'col_name',
        deltap = None, gammap = None, kappap = None):

        """CNLS DDF meta-frontier model.

        First, group frontiers are estimated within each value of ``kind`` and
        their conditional inefficiency estimates are recovered by StoNED. Those
        group inefficiencies are then removed along the DDF direction before the
        common meta-frontier is fitted. This implements the matafrontier idea
        that the meta gap is estimated from the adjusted group-frontier point,
        not by adding group inefficiency as an extra residual term.
        """
        if kind not in data.columns:
            raise ValueError(f"kind column '{kind}' not found in data.")

        nobs = len(data)
        self.gddf = np.zeros(nobs, dtype=float)
        self.gddf_er = [None] * nobs
        self.gresidual = np.zeros(nobs, dtype=float)
        self.group_models = {}

        for kind0 in data[kind].drop_duplicates():
            mask = data[kind] == kind0
            pos = np.flatnonzero(mask.to_numpy())
            data0 = data.loc[mask].copy()
            model0 = CNLSDDFweak(data0, sent=sent, z=z, gy=gy, gx=gx, gb=gb, fun=fun, rts=rts,
                                 deltap=deltap, gammap=gammap, kappap=kappap)
            model0.optimize(solver="mosek")
            rd0 = StoNED.StoNED(model0)
            gddf0 = np.asarray(rd0.get_technical_inefficiency(RED_QLE), dtype=float).reshape(-1)
            gddf_er0 = rd0.get_technical_efficiency(RED_QLE)
            gresidual0 = np.asarray(model0.get_residual(), dtype=float).reshape(-1)

            if len(gddf0) != len(pos):
                raise ValueError("Group inefficiency length does not match group sample size.")
            self.gddf[pos] = gddf0
            self.gresidual[pos] = gresidual0
            if hasattr(gddf_er0, 'iloc'):
                values = gddf_er0.reset_index(drop=True).to_dict('records')
            else:
                values = np.asarray(gddf_er0).reshape(-1).tolist()
            for k, v in zip(pos, values):
                self.gddf_er[k] = v
            self.group_models[kind0] = model0

        self.gddf = self.gddf.tolist()
        self.gresidual = self.gresidual.tolist()

        self.y, self.x, self.b, self.z,self.gy, self.gx, self.gb, self.basexyb, self.basexyb_old \
            = tools.assert_CNLSDDFweakmeta(data,sent,z,self.gddf,gy,gx,gb)
        self.fun = fun
        self.rts = rts

        self.__model__ = ConcreteModel()

        # Initialize the sets
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))
        self.__model__.K = Set(initialize=range(len(self.y[0])))
        self.__model__.L = Set(initialize=range(len(self.b[0])))

        # Initialize the variables
        if self.rts == RTS_VRS1:
            self.__model__.alpha = Var(self.__model__.I, bounds=(0.0, None), doc='alpha')
        elif self.rts == RTS_VRS2:
            self.__model__.alpha = Var(self.__model__.I, bounds=(None, None), doc='alpha')

        if type(deltap) != type(None):
            self.__model__.delta = Var(
                self.__model__.I, self.__model__.J, bounds=(deltap[0], deltap[1]), doc='delta')
        else:
            self.__model__.delta = Var(
                self.__model__.I, self.__model__.J, bounds=(0.0, None), doc='delta')

        if type(gammap) != type(None):
            self.__model__.gamma = Var(
                self.__model__.I, self.__model__.K, bounds=(gammap[0], gammap[1]), doc='gamma')
        else:
            self.__model__.gamma = Var(
                self.__model__.I, self.__model__.K, bounds=(0.0, None), doc='gamma')

        if type(kappap) != type(None):
            self.__model__.kappa = Var(
                self.__model__.I, self.__model__.L, bounds=(kappap[0], kappap[1]), doc='kappa')
        else:
            self.__model__.kappa = Var(
                self.__model__.I, self.__model__.L, bounds=(0.0, None), doc='kappa')

        self.__model__.epsilon = Var(self.__model__.I, bounds=(None, None), doc='residuals')

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
        if self.rts == RTS_VRS2:
            self.__model__.red_factor_rule = Constraint(self.__model__.I,
                                                    self.__model__.I,
                                                    rule=self.__reduction_factor_rule(),
                                                    doc='reduction factor inequality')

        self.__model__.translation_rule = Constraint(self.__model__.I,
                                                     rule=self.__translation_property(),
                                                     doc='translation property')

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
        """Return the meta-frontier regression constraint.

        The group inefficiency has already been removed from x/y/b by
        assert_CNLSDDFweakmeta(...). Therefore the second-stage residual is the
        meta-frontier composite error / technology gap term and should be fitted
        against the adjusted base value, exactly as in the ordinary CNLSDDFweak
        regression.
        """
        if self.rts == RTS_VRS1:
            if type(self.z) != type(None):
                def regression_rule(model, i):
                    return model.epsilon[i] + self.basexyb[i] == \
                            sum(model.delta[i, j] * self.x[i][j] for j in model.J)- \
                            sum(model.gamma[i, k] * self.y[i][k] for k in model.K)+ \
                            sum(model.kappa[i, l] * self.b[i][l] for l in model.L)- \
                            sum(model.omega[m] * self.z[i][m] for m in model.M)
                return regression_rule

            def regression_rule(model, i):
                return model.epsilon[i] + self.basexyb[i] == \
                        sum(model.delta[i, j] * self.x[i][j] for j in model.J)- \
                        sum(model.gamma[i, k] * self.y[i][k] for k in model.K)+ \
                        sum(model.kappa[i, l] * self.b[i][l] for l in model.L)
            return regression_rule

        if self.rts == RTS_VRS2:
            if type(self.z) != type(None):
                def regression_rule(model, i):
                    return model.epsilon[i] + self.basexyb[i] == \
                            sum(model.delta[i, j] * self.x[i][j] for j in model.J)- \
                            sum(model.gamma[i, k] * self.y[i][k] for k in model.K)+ \
                            sum(model.kappa[i, l] * self.b[i][l] for l in model.L)- \
                            sum(model.omega[m] * self.z[i][m] for m in model.M)+model.alpha[i]
                return regression_rule

            def regression_rule(model, i):
                return model.epsilon[i] + self.basexyb[i] == \
                        sum(model.delta[i, j] * self.x[i][j] for j in model.J)- \
                        sum(model.gamma[i, k] * self.y[i][k] for k in model.K)+ \
                        sum(model.kappa[i, l] * self.b[i][l] for l in model.L)+model.alpha[i]
            return regression_rule

        elif self.rts == RTS_CRS:
            if type(self.z) != type(None):
                def regression_rule(model, i):
                    return model.epsilon[i] + self.basexyb[i] == \
                        sum(model.delta[i, j] * self.x[i][j] for j in model.J) - \
                        sum(model.gamma[i, k] * self.y[i][k] for k in model.K) + \
                        sum(model.kappa[i, l] * self.b[i][l] for l in model.L) - \
                        sum(model.omega[m] * self.z[i][m] for m in model.M)
                return regression_rule

            def regression_rule(model, i):
                return model.epsilon[i] + self.basexyb[i] == \
                    sum(model.delta[i, j] * self.x[i][j] for j in model.J) - \
                    sum(model.gamma[i, k] * self.y[i][k] for k in model.K) + \
                    sum(model.kappa[i, l] * self.b[i][l] for l in model.L)
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
                return __operator(-model.alpha[i] \
                                  + sum(model.delta[i, j] * self.x[i][j] for j in model.J) \
                                  - sum(model.gamma[i, k] * self.y[i][k] for k in model.K)  \
                                  + sum(model.kappa[i, l] * self.b[i][l] for l in model.L)
                                  ,
                                  -model.alpha[h]
                                  + sum(model.delta[h, j] * self.x[i][j] for j in model.J) \
                                  - sum(model.gamma[h, k] * self.y[i][k] for k in model.K)  \
                                  + sum(model.kappa[h, l] * self.b[i][l] for l in model.L)
                                  )

            return afriat_rule

        elif self.rts == RTS_VRS2:
            def afriat_rule(model, i, h):
                if i == h:
                    return Constraint.Skip
                return __operator(model.alpha[i] \
                                  + sum(model.delta[i, j] * self.x[i][j] for j in model.J) \
                                  - sum(model.gamma[i, k] * self.y[i][k] for k in model.K)  \
                                  + sum(model.kappa[i, l] * self.b[i][l] for l in model.L)
                                  ,
                                  model.alpha[h]
                                  + sum(model.delta[h, j] * self.x[i][j] for j in model.J) \
                                  - sum(model.gamma[h, k] * self.y[i][k] for k in model.K)  \
                                  + sum(model.kappa[h, l] * self.b[i][l] for l in model.L)
                                  )

            return afriat_rule

        elif self.rts == RTS_CRS:
            def afriat_rule(model, i, h):
                if i == h:
                    return Constraint.Skip
                return __operator(0+ sum(model.delta[i, j] * self.x[i][j] for j in model.J) \
                                  - sum(model.gamma[i, k] * self.y[i][k] for k in model.K)  \
                                  + sum(model.kappa[i, l] * self.b[i][l] for l in model.L), \
                                  0+ sum(model.delta[h, j] * self.x[i][j] for j in model.J) \
                                  - sum(model.gamma[h, k] * self.y[i][k] for k in model.K)  \
                                  + sum(model.kappa[h, l] * self.b[i][l] for l in model.L)     )

            return afriat_rule

        raise ValueError("Undefined model parameters.")

    def __reduction_factor_rule(self):
        """Return the proper reduction_factor inequality constraint"""
        def reduction_factor_rule(model, i, h):
            return model.alpha[h] + sum(model.delta[h, j] *self.x[i][j]  for j in model.J) >= 0

        return reduction_factor_rule


    def __translation_property(self):
        """Return the proper translation property"""

        def translation_rule(model, i):
            return sum(model.delta[i, j] * self.gx[j] for j in model.J) \
                + sum(model.gamma[i, k] * self.gy[k] for k in model.K) \
                + sum(model.kappa[i, l] * self.gb[l] for l in model.L) == 1

        return translation_rule



    # def __translation_property2(self):
    #     """Return the proper translation property"""

    #     def translation_rule2(model, i):
    #         return -self.basexyb[i] >= model.epsilon[i] + self.gddf[i]

    #     return translation_rule2





    def get_residual(self):
        """Return residual value by array"""
        tools.assert_optimized(self.optimization_status)
        residual = list(  self._single_model().epsilon[:].value)
        if sum(self.gy) >= 1:
            residual = [+1*v for v in residual]
        elif sum(self.gx) >= 1 or sum(self.gb) >= 1:
            residual = [-1*v for v in residual]
        return np.asarray(residual)
