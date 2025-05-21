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
        print("gce",self.gce)
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
        # TODO(error/warning handling): Check problem status after optimization
        self.problem_status, self.optimization_status = tools.optimize_model(
            self.__model__, email, self.cet, solver)

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

    def display_status(self):
        """Display the status of problem"""
        tools.assert_optimized(self.optimization_status)
        print(self.display_status)

    def display_alpha(self):
        """Display alpha value"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_various_return_to_scale(self.rts)
        self.__model__.alpha.display()

    def display_delta(self):
        """Display delta value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.delta.display()

    def display_gamma(self):
        """Display gamma value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.gamma.display()

    def display_kappa(self):
        """Display kappa value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.kappa.display()
    def display_omega(self):
        """Display omega value"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_contextual_variable(self.z)
        self.__model__.omega.display()

    def display_residual(self):
        """Dispaly residual value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.epsilon.display()

    def get_status(self):
        """Return status"""
        return self.optimization_status

    def get_alpha(self):
        """Return alpha value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_various_return_to_scale(self.rts)
        alpha = list(self.__model__.alpha[:].value)
        return np.asarray(alpha)

    def get_delta(self):
        """Return delta value by array"""
        tools.assert_optimized(self.optimization_status)
        delta = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.delta),
                                                          list(self.__model__.delta[:, :].value))])
        delta = pd.DataFrame(delta, columns=['Name', 'Key', 'Value'])
        delta = delta.pivot(index='Name', columns='Key', values='Value')
        return delta.to_numpy()

    def get_gamma(self):
        """Return gamma value by array"""
        tools.assert_optimized(self.optimization_status)
        gamma = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.gamma),
                                                          list(self.__model__.gamma[:, :].value))])
        gamma = pd.DataFrame(gamma, columns=['Name', 'Key', 'Value'])
        gamma = gamma.pivot(index='Name', columns='Key', values='Value')
        return gamma.to_numpy()

    def get_kappa(self):
        """Return kappa value by array"""
        tools.assert_optimized(self.optimization_status)
        kappa = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.kappa),
                                                          list(self.__model__.kappa[:, :].value))])
        kappa = pd.DataFrame(kappa, columns=['Name', 'Key', 'Value'])
        kappa = kappa.pivot(index='Name', columns='Key', values='Value')
        return kappa.to_numpy()

    # def get_residual(self):
    #     """Return residual value by array"""
    #     tools.assert_optimized(self.optimization_status)
    #     residual = list(self.__model__.epsilon[:].value)
    #     return np.asarray(residual)

    def get_residual(self):
        """Return residual value by array"""
        tools.assert_optimized(self.optimization_status)
        residual = list(  self.__model__.epsilon[:].value)
        if sum(self.gx) >= 1 and sum(self.gy) ==0 and sum(self.gb) == 0:
            residual = [+1*v for v in residual]
        elif sum(self.gb) >= 1 and sum(self.gx) ==0 and sum(self.gy) == 0:
            residual = [+1*v for v in residual]
        elif sum(self.gy) >= 1 and sum(self.gx) ==0 and sum(self.gb) == 0:
            residual = [-1*v for v in residual]
        # print("aaasss#,",list(self.__model__.epsilon[:].value))
        return np.asarray(residual)


    def get_omega(self):
        """Return omega value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_contextual_variable(self.z)
        omega = list(self.__model__.omega[:].value)
        return np.asarray(omega)

    def get_adjusted_residual(self):
        """Return the shifted residuals(epsilon) tern by CCNLS"""
        tools.assert_optimized(self.optimization_status)
        return self.get_residual() - np.amax(self.get_residual())

    def get_adjusted_alpha(self):
        """Return the shifted constatnt(alpha) term by CCNLS"""
        tools.assert_optimized(self.optimization_status)
        return self.get_alpha() + np.amax(self.get_residual())






class CNLSDDFweakmeta(CNLSDDFweak):
    """Convex Nonparametric Least Square with directional distance function
    """

    def __init__(self, data, sent, z=None, gy=[1], gx=[1], gb=[1], fun=FUN_PROD, rts=RTS_VRS1, kind = 'col_name',
        deltap = None, gammap = None, kappap = None):

        """CNLS DDF model

        Args:
            # y (float): output variable.
            # x (float): input variables.
            # b (float): undesirable output variables.
            data
            sent
            z (float, optional): Contextual variable(s). Defaults to None.
            gy (list, optional): output directional vector. Defaults to [1].
            gx (list, optional): input directional vector. Defaults to [1].
            gb (list, optional): undesirable output directional vector. Defaults to [1].
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS.
            kind(string): 'col_name' or 'constatnt'. Defaults to 'col_name'.
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        kind_lt = data[kind].unique()
        print("kind_lt",kind_lt)
        self.gddf = []
        self.gddf_er = []
        self.gresidual = []

        for kind0 in kind_lt:
            data0 = data[data[kind] == kind0]
            model0 = CNLSDDFweak(data0, sent=sent, z=z, gy=gy, gx=gx, gb=gb, fun=fun, rts=rts,
                                 deltap=deltap, gammap=deltap, kappap=deltap  )
            model0.optimize(solver="mosek")
            rd0 = StoNED.StoNED(model0)
            gddf0 = rd0.get_technical_inefficiency(RED_QLE)
            print("gddf0", len(gddf0))
            self.gddf = np.hstack((self.gddf, gddf0)).tolist()

            gddf_er0 = rd0.get_technical_efficiency_ratio(RED_QLE)
            self.gddf_er = np.hstack((self.gddf_er, gddf_er0)).tolist()

            gresidual0 = model0.get_residual()
            self.gresidual = np.hstack((self.gresidual, gresidual0)).tolist()

        print("gddf_er", self.gddf_er)

        # 将列表转换为NumPy数组并调整形状为一列
        print("gddf", self.gddf)
        print("gddf", len(self.gddf))

        self.y, self.x, self.b, self.z,self.gy, self.gx, self.gb, self.basexyb, self.basexyb_old \
            = tools.assert_CNLSDDFweakmeta(data,sent,z,self.gddf,gy,gx,gb)
        self.fun = fun
        self.rts = rts
        print("basexyb1", (self.basexyb))

        # self.basexyb = [bb+gf for bb,gf in zip(self.basexyb , self.gddf)]
        print("basexyb2", len(self.basexyb))

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


        self.__model__.epsilon = Var(self.__model__.I, bounds=(0.0, None), doc='residuals')


        if type(self.z) != type(None):
            # Initialize the set of z
            self.__model__.M = Set(initialize=range(len(self.z[0])))
            # Initialize the variables for z variable
            self.__model__.omega = Var(self.__model__.M, doc='z coefficient')

        # Setup the objective function and constraints
        self.__model__.objective = Objective(rule=self._CNLSSD__objective_rule(),
                                             sense=minimize,
                                             doc='objective function')
        # self.__model__.objective.pprint()

        self.__model__.regression_rule = Constraint(self.__model__.I,
                                                    rule=self.__regression_rule(),
                                                    doc='regression equation')
        # self.__model__.regression_rule.pprint()
        self.__model__.afriat_rule = Constraint(self.__model__.I,
                                                self.__model__.I,
                                                rule=self.__afriat_rule(),
                                                doc='afriat inequality')
        # self.__model__.afriat_rule.pprint()
        if self.rts == RTS_VRS2:
            self.__model__.red_factor_rule = Constraint(self.__model__.I,
                                                    self.__model__.I,
                                                    rule=self.__reduction_factor_rule(),
                                                    doc='reduction factor inequality')

        self.__model__.translation_rule = Constraint(self.__model__.I,
                                                     rule=self.__translation_property(),
                                                     doc='translation property')
        # if sum(self.gy) < 1:
        #     if (sum(self.gx) >= 1) or (sum(self.gb) >= 1):
        #         # print("aaaaaaaaaaaaaaa")
        #         self.__model__.translation_rule2 = Constraint(self.__model__.I,
        #                                                  rule=self.__translation_property2(),
        #                                                  doc='translation property2')
        # self.__model__.translation_rule2.pprint()

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

            if type(self.z) != type(None):
                def regression_rule(model, i):
                    return model.epsilon[i]   +self.gddf[i] == \
                            sum(model.delta[i, j] * self.x[i][j] for j in model.J)- \
                            sum(model.gamma[i, k] * self.y[i][k] for k in model.K)+ \
                            sum(model.kappa[i, l] * self.b[i][l] for l in model.L)- \
                            sum(model.omega[m] * self.z[i][m] for m in model.M)

                return regression_rule

            elif type(self.z) == type(None):
                def regression_rule(model, i):
                    return model.epsilon[i]   +self.gddf[i] == \
                            sum(model.delta[i, j] * self.x[i][j] for j in model.J)- \
                            sum(model.gamma[i, k] * self.y[i][k] for k in model.K)+ \
                            sum(model.kappa[i, l] * self.b[i][l] for l in model.L)

                return regression_rule

        if self.rts == RTS_VRS2:

            if type(self.z) != type(None):
                def regression_rule(model, i):
                    return model.epsilon[i]  +self.gddf[i]   == \
                            sum(model.delta[i, j] * self.x[i][j] for j in model.J)- \
                            sum(model.gamma[i, k] * self.y[i][k] for k in model.K)+ \
                            sum(model.kappa[i, l] * self.b[i][l] for l in model.L)- \
                            sum(model.omega[m] * self.z[i][m] for m in model.M)+model.alpha[i]

                return regression_rule

            elif type(self.z) == type(None):
                def regression_rule(model, i):
                    return model.epsilon[i]  +self.gddf[i]   == \
                            sum(model.delta[i, j] * self.x[i][j] for j in model.J)- \
                            sum(model.gamma[i, k] * self.y[i][k] for k in model.K)+ \
                            sum(model.kappa[i, l] * self.b[i][l] for l in model.L)+model.alpha[i]

                return regression_rule

        elif self.rts == RTS_CRS:
            if type(self.z) != type(None):
                def regression_rule(model, i):
                    return model.epsilon[i]  +self.gddf[i]  == \
                        sum(model.delta[i, j] * self.x[i][j] for j in model.J) - \
                        sum(model.gamma[i, k] * self.y[i][k] for k in model.K) + \
                        sum(model.kappa[i, l] * self.b[i][l] for l in model.L) - \
                        sum(model.omega[m] * self.z[i][m] for m in model.M)

                return regression_rule

            elif type(self.z) == type(None):
                def regression_rule(model, i):
                    return model.epsilon[i]   +self.gddf[i]  == \
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

    def display_gamma(self):
        """Display gamma value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.gamma.display()

    def display_delta(self):
        """Display delta value"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_undesirable_output(self.b)
        self.__model__.delta.display()

    def get_gamma(self):
        """Return gamma value by array"""
        tools.assert_optimized(self.optimization_status)
        gamma = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.gamma),
                                                           list(self.__model__.gamma[:, :].value))])
        gamma = pd.DataFrame(gamma, columns=['Name', 'Key', 'Value'])
        gamma = gamma.pivot(index='Name', columns='Key', values='Value')
        return gamma.to_numpy()

    def get_delta(self):
        """Return delta value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_undesirable_output(self.b)
        delta = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.delta),
                                                           list(self.__model__.delta[:, :].value))])
        delta = pd.DataFrame(delta, columns=['Name', 'Key', 'Value'])
        delta = delta.pivot(index='Name', columns='Key', values='Value')
        return delta.to_numpy()

    def get_residual(self):
        """Return residual value by array"""
        tools.assert_optimized(self.optimization_status)
        residual = list(  self.__model__.epsilon[:].value)
        if  sum(self.gx) >= 1 and sum(self.gy) ==0 and sum(self.gb) == 0:
            residual = [-1*v for v in residual]
        elif sum(self.gb) >= 1 and sum(self.gx) ==0 and sum(self.gy) == 0:  
            residual = [-1*v for v in residual]
        elif sum(self.gy) >= 1 and sum(self.gx) ==0 and sum(self.gb) == 0:
            residual = [+1*v for v in residual]
        return np.asarray(residual)

