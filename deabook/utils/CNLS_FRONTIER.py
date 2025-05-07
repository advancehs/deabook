# import dependencies
from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, maximize, Constraint
import numpy as np
import pandas as pd
from ..constant import CET_ADDI, ORIENT_IO, ORIENT_OO, RTS_VRS1, RTS_CRS, OPT_DEFAULT, OPT_LOCAL
from . import tools


class DEA:
    """Data Envelopment Analysis (DEA)
    """

    def __init__(self, data, sent, gy=[1], gx=[0], rts=RTS_VRS1, baseindex=None, refindex=None):
        """DEA: Envelopment problem

        Args:
            data
            sent
            gy (list, optional): output distance vector. Defaults to [1].
            gx (list, optional): input distance vector. Defaults to [0].
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.y, self.x, self.yref, self.xref, self.gy, self.gx,self.data_index \
            = tools.assert_DEA(data, sent, gy, gx,baseindex,refindex)
        self.rts = rts

        self.modeldict = {}
        for i in self.data_index:
            self.ii = i
            # Initialize DEA model
            self.__model__ = ConcreteModel()

            # Initialize sets
            self.__model__.R = Set(initialize=range(len(self.yref)))
            self.__model__.J = Set(initialize=range(len(self.x[0])))
            self.__model__.K = Set(initialize=range(len(self.y[0])))

            # Initialize variable
            self.__model__.rho = Var( doc='efficiency')
            self.__model__.lamda = Var(self.__model__.R, bounds=(0.0, None), doc='intensity variables')

            # Setup the objective function and constraints
            if sum(self.gx) >= 1:
                self.__model__.objective = Objective(
                    rule=self.__objective_rule(), sense=minimize, doc='objective function')
            elif sum(self.gy) >=1:
                self.__model__.objective = Objective(
                    rule=self.__objective_rule(), sense=maximize, doc='objective function')
            else:
                raise ValueError("gx and gy must either be 1 or 0.")
            # self.__model__.objective.pprint()

            self.__model__.input = Constraint(
                 self.__model__.J, rule=self.__input_rule(), doc='input constraint')
            # self.__model__.input.pprint()
            self.__model__.output = Constraint(
                 self.__model__.K, rule=self.__output_rule(), doc='output constraint')
            # self.__model__.output.pprint()

            if self.rts == RTS_VRS1:
                self.__model__.vrs = Constraint(
                    rule=self.__vrs_rule(), doc='variable return to scale rule')

            # Optimize model
            self.optimization_status = 0
            self.problem_status = 0
            self.modeldict[i] = self.__model__

    def __objective_rule(self):
        """Return the proper objective function"""
        def objective_rule(model):
            return model.rho
        return objective_rule

    def __input_rule(self):
        """Return the proper input constraint"""
        if sum(self.gx) >= 1:
            def input_rule(model, j):
                if self.gx[j] == 1:

                    return (sum(model.lamda[ r]*self.xref[r][j] for r in model.R) <= \
                            model.rho*self.x[self.ii][j]   )
                else:
                    return (sum(model.lamda[ r]*self.xref[r][j] for r in model.R) <= \
                            self.x[self.ii][j] )
            return input_rule
        elif sum(self.gy) >=1:
            def input_rule(model,  j):
                return sum(model.lamda[ r] * self.xref[r][j] for r in model.R) <= self.x[self.ii][j]
            return input_rule

    def __output_rule(self):
        """Return the proper output constraint"""
        if sum(self.gx) >= 1:
            def output_rule(model,  k):
                return sum(model.lamda[ r] * self.yref[r][k] for r in model.R) >= self.y[self.ii][k]
            return output_rule
        elif sum(self.gy) >= 1:
            def output_rule(model, k):
                if sum(self.gy) >= 1:

                    return (sum(model.lamda[ r]*self.yref[r][k] for r in model.R) >= \
                            model.rho*self.y[self.ii][k]   )
                else:
                    return sum(model.lamda[r] * self.yref[r][k] for r in model.R) >= self.y[self.ii][k]

            return output_rule

    def __vrs_rule(self):
        def vrs_rule(model):
            return sum(model.lamda[ r] for r in model.R) == 1
        return vrs_rule

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Optimize the function by requested method

        Args:
            email (string): The email address for remote optimization. It will optimize locally if OPT_LOCAL is given.
            solver (string): The solver chosen for optimization. It will optimize with default solver if OPT_DEFAULT is given.
        """
        # TODO(error/warning handling): Check problem status after optimization
        self.optimization_status = pd.DataFrame()
        for ind, problem in self.modeldict.items():
            _,self.optimization_status.loc[ind,'optimization_status'] = tools.optimize_model(problem, email, CET_ADDI, solver)

        # self.problem_status, self.optimization_status = tools.optimize_model(
        #     self.__model__, email, CET_ADDI, solver)

    def display_status(self):
        """Display the status of problem"""
        for ind, problem in self.modeldict.items():

            print(self.optimization_status.loc[ind,'optimization_status'])

    def display_rho(self):
        """Display rho value"""
        for ind, problem in self.modeldict.items():
            tools.assert_optimized(self.optimization_status.loc[ind,'optimization_status'])
            problem.rho.display()

    def display_lamda(self):
        """Display lamda value"""
        for ind, problem in self.modeldict.items():
            tools.assert_optimized(self.optimization_status.loc[ind,'optimization_status'])
            problem.lamda.display()


    def get_status(self):
        """Return status"""

        return self.optimization_status

    def get_rho(self):
        """Return rho value by array"""
        rho_lt = []
        for ind, problem in self.modeldict.items():
            tools.assert_optimized(self.optimization_status.loc[ind,'optimization_status'])
            rho = list(problem.rho[:].value)
            rho_lt.append(rho)
        print("rho_lt", len(rho_lt[0]))
        return np.asarray(rho_lt)

    def get_lamda(self):
        """Return lamda value by array"""
        lamda_lt=[]
        for ind, problem in self.modeldict.items():

            tools.assert_optimized(self.optimization_status.loc[ind,'optimization_status'])
            lamda = list(problem.lamda[:].value)
            lamda_lt.append(lamda)

        return np.asarray(lamda_lt)




class DDF(DEA):
    def __init__(self,  data, sent, gy=[1], gx=[1], rts=RTS_VRS1, baseindex=None, refindex=None):
        """DEA: Directional distance function

        Args:
            data
            sent
            gy (list, optional): output directional vector. Defaults to [1].
            gx (list, optional): input directional vector. Defaults to [1].
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        self.y, self.x, self.yref, self.xref, self.gy, self.gx, self.data_index\
            = tools.assert_DDF(data, sent, gy, gx, baseindex, refindex)
        self.rts = rts

        self.modeldict = {}
        for i in self.data_index:
            self.ii = i
            # Initialize DEA model
            self.__model__ = ConcreteModel()
            self.__model__.R = Set(initialize=range(len(self.yref)))

            # Initialize sets
            # self.__model__.I = Set(initialize=range(len(self.y)))
            self.__model__.J = Set(initialize=range(len(self.x[0])))
            self.__model__.K = Set(initialize=range(len(self.y[0])))

            # Initialize variable
            self.__model__.rho = Var( doc='directional distance')

            self.__model__.lamda = Var( self.__model__.R, bounds=(0.0, None), doc='intensity variables')

            # Setup the objective function and constraints
            self.__model__.objective = Objective(rule=self._DEA__objective_rule(), sense=maximize, doc='objective function')
            # self.__model__.objective.pprint()
            self.__model__.input = Constraint( self.__model__.J, rule=self.__input_rule(), doc='input constraint')
            # self.__model__.input.pprint()


            self.__model__.output = Constraint(self.__model__.K, rule=self.__output_rule(), doc='output constraint')

            # self.__model__.output.pprint()

            if self.rts == RTS_VRS1:
                self.__model__.vrs = Constraint( rule=self.__vrs_rule(), doc='various return to scale rule')

            # Optimize model
            self.optimization_status = 0
            self.problem_status = 0
            self.modeldict[i] = self.__model__

    def __input_rule(self):
        """Return the proper input constraint"""
        def input_rule(model,  j):
            return sum(model.lamda[ r] * self.xref[r][j] for r in model.R) <= \
                    self.x[self.ii][j] - model.rho*self.gx[j]*self.x[self.ii][j]

        return input_rule

    def __output_rule(self):
        """Return the proper output constraint"""
        def output_rule(model, k):
            return sum(model.lamda[ r] * self.yref[r][k] for r in model.R) >= \
                        self.y[self.ii][k] + model.rho*self.gy[k]*self.y[self.ii][k]

        return output_rule



    def __vrs_rule(self):
        """Return the VRS constraint"""
        def vrs_rule(model):
            return sum(model.lamda[ r] for r in model.R) == 1
        return vrs_rule





class DEADUAL(DEA):
    def __init__(self, data, sent, muhat, epsilonhat, gy=[1], gx=[0], rts=RTS_VRS1, baseindex=None, refindex=None):
        """DEA: Multiplier problem

        Args:
            data
            sent
            orient (String): ORIENT_IO (input orientation) or ORIENT_OO (output orientation)
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """

        self.y, self.x, self.yref, self.xref, self.gy, self.gx,self.data_index \
            = tools.assert_DEA(data, sent, gy, gx,baseindex,refindex)

        self.rts = rts
        self.muhat = muhat
        self.epsilonhat = epsilonhat

        self.modeldict = {}
        for i in self.data_index:
            self.epsilonhatii = self.epsilonhat[i]
            self.ii = i
            # Initialize DEA model
            self.__model__ = ConcreteModel()
            # Initialize sets
            self.__model__.R = Set(initialize=range(len(self.yref)))
            # self.__model__.I = Set(initialize=range(len(self.y)))
            self.__model__.J = Set(initialize=range(len(self.x[0])))
            self.__model__.K = Set(initialize=range(len(self.y[0])))

            # Initialize variable
            self.__model__.delta = Var( self.__model__.J, bounds=(0.0, None), doc='multiplier x')
            self.__model__.gamma = Var( self.__model__.K, bounds=(0.0, None), doc='multiplier y')
            if self.rts == RTS_VRS1:
                self.__model__.alpha = Var( doc='variable return to scale')

            # Setup the objective function and constraints
            if sum(self.gx)>=1:
                self.__model__.objective = Objective(rule=self.__objective_rule(), sense=minimize, doc='objective function')
            elif sum(self.gy)>=1:
                self.__model__.objective = Objective(rule=self.__objective_rule(), sense=minimize, doc='objective function')
            else:
                raise ValueError("gx and gy must either be 1 or 0.")
            # self.__model__.objective.pprint()
            self.__model__.first = Constraint( self.__model__.R, rule=self.__first_rule(), doc='technology constraint')
            # self.__model__.first.pprint()

            self.__model__.second = Constraint(rule=self.__second_rule(), doc='normalization constraint')
            # self.__model__.second.pprint()
            # Optimize model
            self.optimization_status = 0
            self.problem_status = 0
            self.modeldict[i] = self.__model__







    def __objective_rule(self):
        """Return the proper objective function"""
        if sum(self.gx) >= 1:

            def objective_rule(model):
                if self.rts == RTS_VRS1:
                    return sum(model.delta[ j] * self.x[self.ii][j] * (1-self.gx[j])  for j in model.J) - \
                        sum(model.gamma[ k] * self.y[self.ii][k] * (1-self.gy[k])  for k in model.K) +\
                        model.alpha
                elif self.rts == RTS_CRS:
                    return sum(model.delta[ j] * self.x[self.ii][j] * (1-self.gx[j])  for j in model.J) - \
                        sum(model.gamma[ k] * self.y[self.ii][k] * (1-self.gy[k]) for k in model.K)
            return objective_rule

        elif sum(self.gy) >= 1:

            def objective_rule(model):
                if self.rts == RTS_VRS1:
                    return sum(model.delta[ j] * self.x[self.ii][j] * (1-self.gx[j])  for j in model.J) - \
                        sum(model.gamma[ k] * self.y[self.ii][k] * (1-self.gy[k]) for k in model.K) +\
                        model.alpha
                elif self.rts == RTS_CRS:
                    return sum(model.delta[ j] * self.x[self.ii][j] * (1-self.gx[j])  for j in model.J) - \
                        sum(model.gamma[ k] * self.y[self.ii][k] * (1-self.gy[k]) for k in model.K)

            return objective_rule

    def __first_rule(self):
        """Return the proper technology constraint"""
        if sum(self.gx) >= 1:
            if self.rts == RTS_VRS1:
                def first_rule(model, r):
                    return sum(model.delta[ j] * self.xref[r][j] for j in model.J) \
                        - sum(model.gamma[ k] * self.yref[r][k] for k in model.K) \
                        + model.alpha >= self.muhat+ self.epsilonhatii
                return first_rule
            elif self.rts == RTS_CRS:
                def first_rule(model, r):
                    return sum(model.delta[j] * self.xref[r][j] for j in model.J) \
                        - sum(model.gamma[ k] * self.yref[r][k] for k in model.K) >= self.muhat+ self.epsilonhatii
                return first_rule
        elif sum(self.gy) >= 1:
            if self.rts == RTS_VRS1:
                def first_rule(model,  r):
                    return sum(model.delta[j] * self.xref[r][j] for j in model.J) - \
                        sum(model.gamma[ k] * self.yref[r][k] for k in model.K) \
                        + model.alpha >= self.muhat+ self.epsilonhatii
                return first_rule
            elif self.rts == RTS_CRS:
                def first_rule(model, r):
                    return sum(model.delta[ j] * self.xref[r][j] for j in model.J) - \
                        sum(model.gamma[ k] * self.yref[r][k] for k in model.K) >= self.muhat+ self.epsilonhatii
                return first_rule

    def __second_rule(self):
        """Return the proper normalization constraint"""
        if sum(self.gx) >= 1:
            def second_rule(model):
                return sum(model.delta[ j]  *self.gx[j] for j in model.J) <= 1
            return second_rule
        elif sum(self.gy) >= 1:

            def second_rule(model):
                return sum(model.gamma[ k]  *self.gy[k] for k in model.K) == 1
            return second_rule


    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Optimize the function by requested method

        Args:
            email (string): The email address for remote optimization. It will optimize locally if OPT_LOCAL is given.
            solver (string): The solver chosen for optimization. It will optimize with default solver if OPT_DEFAULT is given.
        """
        # TODO(error/warning handling): Check problem status after optimization
        self.optimization_status =pd.DataFrame()
        for ind, problem in self.modeldict.items():
            self.optimization_status.loc[ind,'optimization_status'] = tools.optimize_model(problem, email, CET_ADDI, solver)






    def display_delta(self):
        """Display delta value"""
        for ind, problem in self.modeldict.items():
            tools.assert_optimized(self.optimization_status.loc[ind,'optimization_status'])
            problem.delta.display()

    def display_gamma(self):
        """Display gamma value"""
        for ind, problem in self.modeldict.items():
            tools.assert_optimized(self.optimization_status.loc[ind,'optimization_status'])
            problem.gamma.display()

    def display_alpha(self):
         """Display omega value"""
         for ind, problem in self.modeldict.items():
             tools.assert_optimized(self.optimization_status)
             tools.assert_various_return_to_scale_alpha(self.rts)
             problem.alpha.display()


    def get_delta(self):
        """Return delta value by array"""
        delta_lt = []
        for ind, problem in self.modeldict.items():
            tools.assert_optimized(self.optimization_status.loc[ind,'optimization_status'])
            delta = list(problem.delta[:].value)
            delta_lt.append(delta)
        return np.asarray(delta_lt)


    def get_gamma(self):
        """Return nu value by array"""
        gamma_lt = []
        for ind, problem in self.modeldict.items():
            tools.assert_optimized(self.optimization_status.loc[ind,'optimization_status'])
            gamma = list(problem.gamma[:].value)
            gamma_lt.append(gamma)
        return np.asarray(gamma_lt)

    def get_alpha(self):
        """Return omega value by array"""

        alpha_lt = []
        for ind, problem in self.modeldict.items():
            tools.assert_optimized(self.optimization_status.loc[ind,'optimization_status'])
            tools.assert_various_return_to_scale_alpha(self.rts)

            alpha = list(problem.alpha[:].value)
            alpha_lt.append(alpha)
        return np.asarray(alpha_lt)

    def get_obj(self):
        """Return obj value by array"""
        obj_lt = []
        for ind, problem in self.modeldict.items():
            tools.assert_optimized(self.optimization_status.loc[ind,'optimization_status'])
            obj = problem.objective()
            obj_lt.append(obj)
        if sum(self.gx) >= 1:
            return -np.asarray(obj_lt)
        return np.asarray(obj_lt)


class DDFDUAL(DEADUAL):

    def __init__(self,  data, sent, muhat, epsilonhat, gy=[1], gx=[1], rts=RTS_VRS1, baseindex=None, refindex=None):
        """DEA: Directional distance function

        Args:
            data
            sent
            gy (list, optional): output directional vector. Defaults to [1].
            gx (list, optional): input directional vector. Defaults to [1].
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2009,2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """

        self.y, self.x, self.yref, self.xref, self.gy, self.gx, self.data_index\
            = tools.assert_DDF(data, sent, gy, gx, baseindex, refindex)
        self.rts = rts
        self.muhat = muhat
        self.epsilonhat = epsilonhat

        self.modeldict = {}
        for i in self.data_index:
            self.epsilonhatii = self.epsilonhat[i]
            self.ii = i
            # Initialize DEA model
            self.__model__ = ConcreteModel()

            # Initialize sets
            self.__model__.R = Set(initialize=range(len(self.yref)))
            # self.__model__.I = Set(initialize=range(len(self.y)))
            self.__model__.J = Set(initialize=range(len(self.x[0])))
            self.__model__.K = Set(initialize=range(len(self.y[0])))

            # Initialize variable
            self.__model__.delta = Var( self.__model__.J, bounds=(0.0, None), doc='multiplier x')
            self.__model__.gamma = Var( self.__model__.K, bounds=(0.0, None), doc='multiplier y')

            if self.rts == RTS_VRS1:
                self.__model__.alpha = Var( doc='variable return to scale')

            # Setup the objective function and constraints
            self.__model__.objective = Objective(
                rule=self.__objective_rule(), sense=minimize, doc='objective function')
            self.__model__.first = Constraint(
                self.__model__.R, rule=self.__first_rule(), doc='technology constraint')
            self.__model__.second = Constraint(
                 rule=self.__second_rule(), doc='normalization constraint')


            # Optimize model
            self.optimization_status = 0
            self.problem_status = 0
            self.modeldict[i] = self.__model__

    def __objective_rule(self):
        """Return the proper objective function"""
        def objective_rule(model):
            if self.rts == RTS_VRS1:
                return sum(model.delta[ j] * self.x[self.ii][j] for j in model.J) - \
                    sum(model.gamma[ k] * self.y[self.ii][k] for k in model.K) + \
                    model.alpha
            elif self.rts == RTS_CRS:
                return sum(model.delta[ j] * self.x[self.ii][j] for j in model.J) - \
                    sum(model.gamma[ k] * self.y[self.ii][k] for k in model.K)
        return objective_rule


    def __first_rule(self):
        """Return the proper technology constraint"""
        if self.rts == RTS_VRS1:
            def first_rule(model,  r):
                return sum(model.delta[ j] * self.xref[r][j]-(self.muhat +self.epsilonhatii)*self.gx[j] for j in model.J) - \
                    sum(model.gamma[ k] * self.yref[r][k]+(self.muhat + self.epsilonhatii)*self.gy[k] for k in model.K) +\
                    model.alpha >=0
            return first_rule
        elif self.rts == RTS_CRS:
            def first_rule(model,  r):
                return sum(model.delta[ j] * self.xref[r][j]- (self.muhat +self.epsilonhatii)*self.gx[j] for j in model.J) - \
                    sum(model.gamma[ k] * self.yref[r][k]+(self.muhat + self.epsilonhatii)*self.gy[k] for k in model.K) >=  0
            return first_rule


    def __second_rule(self):
        """Return the proper normalization constraint"""
        def second_rule(model):
            return sum(model.delta[ j] *self.gx[j]  for j in model.J) + \
                sum(model.gamma[ k] *self.gy[k]   for k in model.K) == 1
        return second_rule

    def get_obj(self):
        """Return obj value by array"""
        obj_lt = []
        for ind, problem in self.modeldict.items():
            tools.assert_optimized(self.optimization_status.loc[ind,'optimization_status'])

            # print("aaaaaaaa",self.problem_status)
            obj = problem.objective()
            obj_lt.append(obj)

        return np.asarray(obj_lt)
