# import dependencies
from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, maximize, Constraint
import numpy as np
import pandas as pd
from .constant import CET_ADDI, ORIENT_IO, ORIENT_OO, ORIENT_BO, RTS_VRS1, RTS_VRS2, RTS_CRS, OPT_DEFAULT, OPT_LOCAL
from .utils import tools
from .DEA import DEA


class DEAweak(DEA):
    """weak dispsbnility of Data Envelopment Analysis (DEA)
    """

    def __init__(self, data, sent, gy, gx,gb , rts, baseindex=None, refindex=None):
        """DEA: Envelopment problem

        Args:
            data
            sent
            gy (list, optional): output distance vector. Defaults to [1].
            gx (list, optional): input distance vector. Defaults to [0].
            gb (list, optional): undesirable output directional vector. Defaults to [1].
            rts (String): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.y, self.x, self.b,  self.gy, self.gx, self.gb, self.yref, self.xref, self.bref = \
            tools.assert_DEAweak(data, sent, gy, gx, gb, baseindex, refindex)
        self.rts = rts

        # Initialize DEA model
        self.__model__ = ConcreteModel()
        self.__model__.R = Set(initialize=range(len(self.yref)))

        # Initialize sets
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))
        self.__model__.K = Set(initialize=range(len(self.y[0])))
        self.__model__.L = Set(initialize=range(len(self.b[0])))

        # Initialize variable
        self.__model__.rho = Var(self.__model__.I, doc='efficiency')
        if self.rts == RTS_VRS1:
            self.__model__.theta = Var(self.__model__.I, bounds=(0.0, 1.0),doc='emission reduction factor')
        elif self.rts == RTS_VRS2:
            self.__model__.mu = Var(self.__model__.I, self.__model__.R, bounds=(0.0, None),doc='emission reduction factor2')

        self.__model__.lamda = Var(self.__model__.I, self.__model__.R, bounds=(0.0, None), doc='intensity variables')

        # Setup the objective function and constraints
        if sum(self.gy) >= 1:
            self.__model__.objective = Objective(
                rule=self.__objective_rule(), sense=maximize, doc='objective function')
        else:
            self.__model__.objective = Objective(
                rule=self.__objective_rule(), sense=minimize, doc='objective function')

        self.__model__.input = Constraint(
            self.__model__.I, self.__model__.J, rule=self.__input_rule(), doc='input constraint')
        self.__model__.output = Constraint(
            self.__model__.I, self.__model__.K, rule=self.__output_rule(), doc='output constraint')
        self.__model__.unoutput = Constraint(
            self.__model__.I, self.__model__.L, rule=self.__unoutput_rule(), doc='undesirable output constraint')
        if self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
            self.__model__.vrs = Constraint(self.__model__.I, rule=self.__vrs_rule(), doc='variable return to scale rule')

        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0

    def __objective_rule(self):
        """Return the proper objective function"""
        def objective_rule(model):
            return sum(model.rho[i] for i in model.I)
        return objective_rule

    def __input_rule(self):
        """Return the proper input constraint"""
        if sum(self.gx)>=1:
            if self.rts == RTS_VRS2:
                def input_rule(model, o, j):
                    if self.gx[j]==1:
                        return sum((model.lamda[o, r]+model.mu[o,r])*self.xref[r][j] for r in model.R)<=model.rho[o]*self.x[o][j]
                    else:
                        return sum((model.lamda[o, r]+model.mu[o,r])*self.xref[r][j] for r in model.R)<=self.x[o][j]

                return input_rule
            else:
                def input_rule(model, o, j):
                    if self.gx[j]==1:
                        return sum(model.lamda[o, r]*self.xref[r][j] for r in model.R)<=model.rho[o]*self.x[o][j]
                    else:
                        return sum(model.lamda[o, r]*self.xref[r][j] for r in model.R)<=self.x[o][j]
                return input_rule
        else:
            if self.rts == RTS_VRS2:
                def input_rule(model, o, j):
                    return sum((model.lamda[o, r]+model.mu[o,r])*self.xref[r][j] for r in model.R)<=self.x[o][j]
                return input_rule
            else:
                def input_rule(model, o, j):
                    return sum(model.lamda[o, r] * self.xref[r][j] for r in model.R) <= self.x[o][j]

                return input_rule



    def __output_rule(self):
        """Return the proper output constraint"""
        if sum(self.gy)>=1:

            def output_rule(model, o, k):
                if self.gy[k] == 1:
                    return sum(model.lamda[o, r]*self.yref[r][k] for r in model.R) >= model.rho[o]*self.y[o][k]
                else:
                    return sum(model.lamda[o, r] * self.yref[r][k] for r in model.R) >= self.y[o][k]

            return output_rule

        else:
            def output_rule(model, o, k):
                return sum(model.lamda[o, r] * self.yref[r][k] for r in model.R) >= self.y[o][k]
            return output_rule

    def __unoutput_rule(self):
        """Return the proper undesirable output constraint"""
        if sum(self.gb)>=1:
            def unoutput_rule(model, o, l):
                if self.gb[l] == 1:
                    return sum(model.lamda[o, r]*self.bref[r][l] for r in model.R) == model.rho[o]*self.b[o][l]
                else:
                    return sum(model.lamda[o, r] * self.bref[r][l] for r in model.R) == self.b[o][l]

            return unoutput_rule
        else:
            def unoutput_rule(model, o, l):
                return sum(model.lamda[o, r] * self.bref[r][l] for r in model.R) == self.b[o][l]
            return unoutput_rule

    def __vrs_rule(self):
        if self.rts==RTS_VRS1:
            def vrs_rule(model, o):
                return sum(model.lamda[o, r] for r in model.R) == model.theta[o]
            return vrs_rule

        elif self.rts==RTS_VRS2:
            def vrs_rule(model, o):
                return sum((model.lamda[o, r]+model.mu[o, r] )for r in model.R) == 1

            return vrs_rule

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Optimize the function by requested method

        Args:
            email (string): The email address for remote optimization. It will optimize locally if OPT_LOCAL is given.
            solver (string): The solver chosen for optimization. It will optimize with default solver if OPT_DEFAULT is given.
        """
        # TODO(error/warning handling): Check problem status after optimization
        self.problem_status, self.optimization_status = tools.optimize_model(
            self.__model__, email, CET_ADDI, solver)

    def display_status(self):
        """Display the status of problem"""
        print(self.optimization_status)

    def display_theta(self):
        """Display theta value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.theta.display()

    def display_rho(self):
        """Display rho value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.rho.display()

    def display_lamda(self):
        """Display lamda value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.lamda.display()

    def get_status(self):
        """Return status"""
        return self.optimization_status

    def get_theta(self):
        """Return theta value by array"""
        tools.assert_optimized(self.optimization_status)
        theta = list(self.__model__.theta[:].value)
        return np.asarray(theta)

    def get_rho(self):
        """Return rho value by array"""
        tools.assert_optimized(self.optimization_status)
        rho = list(self.__model__.rho[:].value)
        return np.asarray(rho)

    def get_lamda(self):
        """Return lamda value by array"""
        tools.assert_optimized(self.optimization_status)
        lamda = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.lamda),
                                                          list(self.__model__.lamda[:, :].value))])
        lamda = pd.DataFrame(lamda, columns=['Name', 'Key', 'Value'])
        lamda = lamda.pivot(index='Name', columns='Key', values='Value')
        return lamda.to_numpy()


class DDFweak(DEAweak):
    def __init__(self,  data, sent, gy=[1], gx=[1], gb=[1], rts=RTS_VRS1, baseindex=None, refindex=None):
        """DDFweak: Directional distance function with undesirable output

        Args:
            data
            sent
            gy (list, optional): output directional vector. Defaults to [1].
            gx (list, optional): input directional vector. Defaults to [1].
            gb (list, optional): undesirable output directional vector. Defaults to [1].
            rts (String): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """

        self.y, self.x, self.b,  self.gy, self.gx, self.gb, self.yref, self.xref, self.bref = \
            tools.assert_DDFweak(data,sent, gy, gx, gb,baseindex,refindex)
        self.rts = rts

        # Initialize DEA model
        self.__model__ = ConcreteModel()

        # Initialize sets
        self.__model__.R = Set(initialize=range(len(self.yref)))
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))
        self.__model__.K = Set(initialize=range(len(self.y[0])))
        self.__model__.L = Set(initialize=range(len(self.b[0])))

        # Initialize variable

        self.__model__.rho = Var(self.__model__.I, doc='directional distance')
        if self.rts == RTS_VRS1:
            self.__model__.theta = Var(self.__model__.I, bounds=(0.0, 1.0),doc='emission reduction factor')
        elif self.rts == RTS_VRS2:
            self.__model__.mu = Var(self.__model__.I, self.__model__.R, bounds=(0.0, None),
                                       doc='emission reduction factor2')

        self.__model__.lamda = Var(self.__model__.I, self.__model__.R, bounds=(0.0, None), doc='intensity variables')

        # Setup the objective function and constraints
        self.__model__.objective = Objective(
            rule=self._DEAweak__objective_rule(), sense=maximize, doc='objective function')
        self.__model__.input = Constraint(
            self.__model__.I, self.__model__.J, rule=self.__input_rule(), doc='input constraint')
        self.__model__.output = Constraint(
            self.__model__.I, self.__model__.K, rule=self.__output_rule(), doc='output constraint')
        self.__model__.undesirable_output = Constraint(
            self.__model__.I, self.__model__.L, rule=self.__undesirable_output_rule(), doc='undesirable output constraint')

        if self.rts == RTS_VRS1 or self.rts == RTS_VRS2:
            self.__model__.vrs = Constraint(self.__model__.I, rule=self.__vrs_rule(), doc='various return to scale rule')

        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0

    def __input_rule(self):
        """Return the proper input constraint"""
        if self.rts == RTS_VRS2:

            def input_rule(model, o, j):
                return sum((model.lamda[o, r]+model.mu[o, r]) * self.xref[r][j] for r in model.R) \
                            <= self.x[o][j] - model.rho[o]*self.gx[j]*self.x[o][j]
            return input_rule
        else:

            def input_rule(model, o, j):
                return sum(model.lamda[o, r] * self.xref[r][j] for r in model.R) \
                            <= self.x[o][j] - model.rho[o]*self.gx[j]*self.x[o][j]
            return input_rule

    def __output_rule(self):
        """Return the proper output constraint"""
        def output_rule(model, o, k):
            return sum(model.lamda[o, r] * self.yref[r][k] for r in model.R) >= self.y[o][k] + model.rho[o]*self.gy[k]*self.y[o][k]
        return output_rule

    def __undesirable_output_rule(self):
        """Return the proper undesirable output constraint"""
        def undesirable_output_rule(model, o, l):
            return sum(model.lamda[o, r] * self.bref[r][l] for r in model.R) == self.b[o][l] - model.rho[o]*self.gb[l]*self.b[o][l]
        return undesirable_output_rule

    def __vrs_rule(self):
        """Return the VRS constraint"""
        if self.rts == RTS_VRS1:

            def vrs_rule(model, o):
                return sum(model.lamda[o, r] for r in model.R) == model.theta[o]
            return vrs_rule
        elif self.rts == RTS_VRS2:

            def vrs_rule(model, o):
                return sum((model.lamda[o, r]+model.mu[o, r] )for r in model.R) == 1
            return vrs_rule


class DEAweakDUAL(DEAweak):

    def __init__(self, data, sent, gy, gx,gb , rts, baseindex=None, refindex=None):
        """DEA: Envelopment problem

        Args:
            data
            sent
            gy (list, optional): output distance vector. Defaults to [1].
            gx (list, optional): input distance vector. Defaults to [0].
            gb (list, optional): undesirable output directional vector. Defaults to [1].
            rts (String): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.y, self.x, self.b,  self.gy, self.gx, self.gb, self.yref, self.xref, self.bref = \
            tools.assert_DEAweak(data, sent, gy, gx, gb, baseindex, refindex)
        self.rts = rts


        # Initialize DEA model
        self.__model__ = ConcreteModel()

        # Initialize sets
        self.__model__.R = Set(initialize=range(len(self.yref)))
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))
        self.__model__.K = Set(initialize=range(len(self.y[0])))
        self.__model__.L = Set(initialize=range(len(self.b[0])))

        # Initialize variable
        self.__model__.delta = Var(self.__model__.I, self.__model__.J, bounds=(
            0.0, None), doc='multiplier x')
        self.__model__.gamma = Var(self.__model__.I, self.__model__.K, bounds=(
            0.0, None), doc='multiplier y')
        self.__model__.kappa = Var(self.__model__.I, self.__model__.L,  doc='multiplier b')

        if self.rts == RTS_VRS1:
            self.__model__.alpha = Var(self.__model__.I,               bounds=(
            0.0, None), doc='variable return to scale')
        elif self.rts == RTS_VRS2:
            self.__model__.alpha = Var(self.__model__.I,               bounds=(
            None, None), doc='variable return to scale')
        # Setup the objective function and constraints
        self.__model__.objective = Objective(
            rule=self.__objective_rule(), sense=minimize, doc='objective function')

        self.__model__.first = Constraint(
            self.__model__.I, self.__model__.R, rule=self.__first_rule(), doc='technology constraint')
        self.__model__.second = Constraint(
            self.__model__.I,                   rule=self.__second_rule(), doc='normalization constraint')
        if self.rts == RTS_VRS2:
            self.__model__.third = Constraint(
                self.__model__.I, rule=self.__third_rule(), doc='weak disposability constraint')
        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0

    def __objective_rule(self):
        """Return the proper objective function"""
        def objective_rule(model):
            return sum(sum(model.delta[o, j] * self.x[o][j]* (1-self.gx[j]) for o in model.I) for j in model.J) - \
                sum(sum(model.gamma[o, k] * self.y[o][k]* (1-self.gb[k]) for o in model.I) for k in model.K) +\
                sum(sum(model.kappa[o,l] * self.b[o][l]* (1-self.gb[l]) for o in model.I) for l in model.L)

        return objective_rule


    def __first_rule(self):
        """Return the proper technology constraint"""
        if self.rts == RTS_VRS1:
            def first_rule(model, o, r):
                return sum(model.delta[o, j] * self.xref[r][j] for j in model.J) -\
                    sum(model.gamma[o, k] * self.yref[r][k] for k in model.K) + \
                    sum(model.kappa[o, l] * self.bref[r][l] for l in model.L) - model.alpha[o] >= 0
            return first_rule

        elif self.rts == RTS_VRS2:
            def first_rule(model, o, r):
                return sum(model.delta[o, j] * self.xref[r][j] for j in model.J) -\
                    sum(model.gamma[o, k] * self.yref[r][k] for k in model.K) + \
                    sum(model.kappa[o, l] * self.bref[r][l] for l in model.L) + model.alpha[o] >= 0
            return first_rule

        elif self.rts == RTS_CRS:
            def first_rule(model, o, r):
                return sum(model.delta[o, j] * self.xref[r][j] for j in model.J) -\
                    sum(model.gamma[o, k] * self.yref[r][k] for k in model.K) + \
                    sum(model.kappa[o, l] * self.bref[r][l] for l in model.L) >= 0
            return first_rule

    def __second_rule(self):
        """Return the proper normalization constraint"""
        if sum(self.gx) >= 1:
            def second_rule(model, o):
                return sum(model.delta[o, j] * self.x[o][j]* self.gx[j] for j in model.J) == 1
            return second_rule
        elif sum(self.gy) >= 1:

            def second_rule(model, o):
                return sum(model.gamma[o, k] * self.y[o][k]* self.gy[k] for k in model.K) == 1
            return second_rule
        elif sum(self.gb) >= 1:

            def second_rule(model, o):
                return sum(model.kappa[o, l] * self.b[o][l]* self.gb[l] for l in model.L) == 1
            return second_rule

    def __third_rule(self):
        """Return the proper weak disposability constraint"""

        def third_rule(model, o, r):
            return sum(model.delta[o, j] * self.xref[r][j] for j in model.J) + model.alpha[o] >= 0

        return third_rule

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

    def display_alpha(self):
         """Display omega value"""
         tools.assert_optimized(self.optimization_status)
         tools.assert_various_return_to_scale_alpha(self.rts)
         self.__model__.alpha.display()

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

    def get_alpha(self):
        """Return omega value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_various_return_to_scale_alpha(self.rts)
        alpha = list(self.__model__.alpha[:].value)
        return np.asarray(alpha)

    def get_efficiency(self):
        """Return efficiency value by array"""
        tools.assert_optimized(self.optimization_status)
        if sum(self.gx) >= 1:
            if self.rts == RTS_CRS:
                return (np.sum(self.get_delta()*self.y, axis=1)).reshape(len(self.y), 1)
            else:
                return (np.sum(self.get_delta()*self.y, axis=1)).reshape(len(self.y), 1) + self.get_alpha().reshape(len(self.y), 1)
        if sum(self.gy) >= 1:
            if self.rts == RTS_CRS:
                return (np.sum(self.get_gamma()*self.x, axis=1)).reshape(len(self.x), 1)
            else:
                return (np.sum(self.get_gamma()*self.x, axis=1)).reshape(len(self.x), 1) + self.get_alpha().reshape(len(self.x), 1)
        if sum(self.gb) >= 1:
            if self.rts == RTS_CRS:
                return (np.sum(self.get_delta()*self.y, axis=1)).reshape(len(self.y), 1)
            else:
                return (np.sum(self.get_delta()*self.y, axis=1)).reshape(len(self.y), 1) + self.get_alpha().reshape(len(self.y), 1)




class DDFweakDUAL(DEAweakDUAL):
    def __init__(self, y, x, b, gy=[1], gx=[1], gb=[1], rts=RTS_VRS1, yref=None, xref=None, bref=None):
        """DDFweak: Directional distance function with undesirable output

        Args:
            data
            sent
            gy (list, optional): output directional vector. Defaults to [1].
            gx (list, optional): input directional vector. Defaults to [1].
            gb (list, optional): undesirable output directional vector. Defaults to [1].
            rts (String): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """

        self.y, self.x, self.b,  self.gy, self.gx, self.gb, self.yref, self.xref, self.bref = \
            tools.assert_DDFweak(data,sent, gy, gx, gb,baseindex,refindex)
        self.rts = rts


        # Initialize DEA model
        self.__model__ = ConcreteModel()
        self.__model__.R = Set(initialize=range(len(self.yref)))

        # Initialize sets
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))
        self.__model__.K = Set(initialize=range(len(self.y[0])))
        self.__model__.L = Set(initialize=range(len(self.b[0])))

        # Initialize variable
        self.__model__.delta = Var(self.__model__.I, self.__model__.J, bounds=(
            0.0, None), doc='multiplier x')
        self.__model__.gamma = Var(self.__model__.I, self.__model__.K, bounds=(
            0.0, None), doc='multiplier y')
        self.__model__.kappa = Var(self.__model__.I, self.__model__.L,  doc='multiplier b')

        if self.rts == RTS_VRS1:
            self.__model__.alpha = Var(self.__model__.I,               bounds=(
            0.0, None), doc='variable return to scale')
        elif self.rts == RTS_VRS2:
            self.__model__.alpha = Var(self.__model__.I,               bounds=(
            None, None), doc='variable return to scale')

        # Setup the objective function and constraints
        self.__model__.objective = Objective(
            rule=self.__objective_rule(), sense=minimize, doc='objective function')

        self.__model__.first = Constraint(
            self.__model__.I, self.__model__.R, rule=self.__first_rule(), doc='technology constraint')
        self.__model__.second = Constraint(
            self.__model__.I, rule=self.__second_rule(), doc='normalization constraint')
        if self.rts == RTS_VRS2:
            self.__model__.third = Constraint(
                self.__model__.I, rule=self.__third_rule(), doc='weak disposability constraint')
        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0

    def __objective_rule(self):
        """Return the proper objective function"""
        if self.rts == RTS_VRS1:

            def objective_rule(model):
                return sum(sum(model.delta[o, j] * self.x[o][j] for o in model.I) for j in model.J) -\
                    sum(sum(model.gamma[o, k] * self.y[o][k] for o in model.I) for k in model.K) + \
                    sum(sum(model.kappa[o, l] * self.b[o][l] for o in model.I) for l in model.L)
            return objective_rule

        elif self.rts == RTS_VRS2:

            def objective_rule(model):
                return sum(sum(model.delta[o, j] * self.x[o][j] for o in model.I) for j in model.J) -\
                    sum(sum(model.gamma[o, k] * self.y[o][k] for o in model.I) for k in model.K) + \
                    sum(sum(model.kappa[o, l] * self.b[o][l] for o in model.I) for l in model.L) + \
                    sum(model.alpha[o] for o in model.I)
            return objective_rule
    def __first_rule(self):
        """Return the proper technology constraint"""
        if self.rts == RTS_VRS1:
            def first_rule(model, o, r):
                return sum(model.delta[o, j] * self.xref[r][j] for j in model.J) -\
                    sum(model.gamma[o, k] * self.yref[r][k] for k in model.K) + \
                    sum(model.kappa[o, l] * self.bref[r][l] for l in model.L) - model.alpha[o] >= 0
            return first_rule
        elif self.rts == RTS_VRS2:
            def first_rule(model, o, r):
                return sum(model.delta[o, j] * self.xref[r][j] for j in model.J) -\
                    sum(model.gamma[o, k] * self.yref[r][k] for k in model.K) + \
                    sum(model.kappa[o, l] * self.bref[r][l] for l in model.L) + model.alpha[o] >= 0
            return first_rule
        elif self.rts == RTS_CRS:
            def first_rule(model, o, r):
                return sum(model.delta[o, j] * self.xref[r][j] for j in model.J) -\
                    sum(model.gamma[o, k] * self.yref[r][k] for k in model.K) + \
                    sum(model.kappa[o, l] * self.bref[r][l] for l in model.L) >= 0
            return first_rule


    def __second_rule(self):
        """Return the proper normalization constraint"""
        def second_rule(model, o):
            return sum(model.delta[o, j] *self.gx[j]* self.x[o][j] for j in model.J) +\
                sum(model.gamma[o, k] *self.gy[k]* self.y[o][k] for k in model.K) +\
                sum(model.kappa[o, l] *self.gb[l]* self.b[o][l] for l in model.L) == 1
        return second_rule

    def __third_rule(self):
        """Return the proper weak disposability constraint"""

        def third_rule(model, o, r):
            return sum(model.delta[o, j] * self.xref[r][j] for j in model.J) + model.alpha[o] >= 0

        return third_rule

    def display_delta(self):
        """Display delta value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.delta.display()

    def display_gamma(self):
        """Display gamma value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.gamma.display()

    def display_alpha(self):
         """Display omega value"""
         tools.assert_optimized(self.optimization_status)
         tools.assert_various_return_to_scale_alpha(self.rts)
         self.__model__.alpha.display()

    def get_delta(self):
        """Return delta value by array"""
        tools.assert_optimized(self.optimization_status)
        delta = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.delta),
                                                          list(self.__model__.delta[:, :].value))])
        delta = pd.DataFrame(delta, columns=['Name', 'Key', 'Value'])
        delta = delta.pivot(index='Name', columns='Key', values='Value')
        return delta.to_numpy()

    def get_gamma(self):
        """Return nu value by array"""
        tools.assert_optimized(self.optimization_status)
        gamma = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.gamma),
                                                          list(self.__model__.gamma[:, :].value))])
        gamma = pd.DataFrame(gamma, columns=['Name', 'Key', 'Value'])
        gamma = gamma.pivot(index='Name', columns='Key', values='Value')
        return gamma.to_numpy()

    def get_alpha(self):
        """Return omega value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_various_return_to_scale_alpha(self.rts)
        alpha = list(self.__model__.alpha[:].value)
        return np.asarray(alpha)


