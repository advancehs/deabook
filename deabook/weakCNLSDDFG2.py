# import dependencies
import numpy as np
import pandas as pd
from .utils import   weakCNLSDDFZG1, weakCNLSDDFZG2, sweet, tools
from .constant import CET_ADDI, FUN_PROD, FUN_COST, OPT_DEFAULT, RTS_CRS, RTS_VRS,OPT_LOCAL
import time
from . import weakCNLS




class weakCNLSDDFG2():
    """Convex Nonparametric Least Square with weak disposability (weakCNLSDDF) and Genetic algorithm
    """
    def __init__(self, data,sent, z=None, gy=[1], gx=[1], gb=[1], fun=FUN_PROD, rts=RTS_VRS, \
                 baseindex=None, refindex=None):
        """weakCNLSG model

        Args:
            data (pandas.DataFrame): input pandas.
            sent (str): inputvars=outputvars: unoutputvars. e.g.: "K L = Y : CO2"
            z (ndarray, optional): Contextual variable(s). Defaults to None.
            gy (list, optional): output directional vector. Defaults to [1].
            gx (list, optional): input directional vector. Defaults to [1].
            gb (list, optional): undesirable output directional vector. Defaults to [1].
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS.
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        self.outputvars,self.inputvars,self.unoutputvars,self.zvars,self.gy, self.gx, self.gb \
            = tools.assert_valid_yxbz(sent,z,gy,gx,gb)
        self.y, self.x, self.b, self.z, self.yref, self.xref, self.bref, self.zref,self.referenceflag\
            = tools.assert_valid_yxbz2(baseindex,refindex,data,\
                                       self.outputvars,self.inputvars,self.unoutputvars,self.zvars)

        self.xcol = self.x.columns
        self.ycol = self.y.columns
        self.bcol = self.b.columns
        self.zcol = self.z.columns if type(z) != type(None) else None

        self.fun = fun
        self.rts = rts

        self.data = data
        self.sent = sent
        self.baseindex = baseindex
        self.refindex = refindex


        # self.cutactive = sweet.sweet(pd.concat([self.x,self.b],axis=1))
        self.cutactive = sweet.sweet(self.x)
        # self.cutactive_ref = sweet.sweetref(pd.concat([self.x,self.b],axis=1),pd.concat([self.xref,self.bref],axis=1))


        # active (added) violated concavity constraint by iterative procedure
        self.active = pd.DataFrame(np.zeros((self.x.shape[0], self.x.shape[0])),index=self.x.index,columns = self.x.index)
        # violated concavity constraint
        self.active2 = pd.DataFrame(np.zeros((self.x.shape[0], self.x.shape[0])),index=self.x.index,columns = self.x.index)
        # print("active2",self.active2)
        # active (added) violated concavity constraint for weak disposbility constrains by iterative procedure
        self.activeweak = pd.DataFrame(np.zeros((self.x.shape[0], self.x.shape[0])),index=self.x.index,columns = self.x.index)
        # violated concavity constraint for weak disposbility constrains
        self.activeweak2 = pd.DataFrame(np.zeros((self.x.shape[0], self.x.shape[0])),index=self.x.index,columns = self.x.index)

        if self.referenceflag:
            self.activeref = pd.DataFrame(np.zeros((self.x.shape[0], self.xref.shape[0])),index=self.x.index,columns = self.xref.index)
            self.activeref2 = pd.DataFrame(np.zeros((self.x.shape[0], self.xref.shape[0])),index=self.x.index,columns = self.xref.index)
            self.activeweakref = pd.DataFrame(np.zeros((self.x.shape[0], self.xref.shape[0])),index=self.x.index,columns = self.xref.index)
            self.activeweakref2 = pd.DataFrame(np.zeros((self.x.shape[0], self.xref.shape[0])),index=self.x.index,columns = self.xref.index)
        else:
            self.activeref = None
            self.activeref2 = None
            self.activeweakref = None
            self.activeweakref2 = None

        print("xcol,ycol,bcol are:",self.x.columns,self.y.columns,self.b.columns)

        print("gx,gy,gb are:",self.gx,self.gy,self.gb)
        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0



    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Optimize the function by requested method"""
        # TODO(error/warning handling): Check problem status after optimization
        self.t0 = time.time()

        model1 = weakCNLSDDFZG1.weakCNLSDDFZG1(
            self.data, self.sent, self.z, self.cutactive, self.gy, self.gx, self.gb, \
                    self.fun, self.rts, self.baseindex, self.refindex)

        model1.optimize(email, solver)
        self.alpha = model1.get_alpha()
        self.beta = model1.get_beta()
        self.gamma = model1.get_gamma()
        self.delta = model1.get_delta()
        self.__model__ = model1.__model__

        self.count = 0
        if self.referenceflag:
            while self.__convergence_test() > 0.0001 \
                    or self.__convergence_test_weak() > 0.0001 \
                    or self.__convergence_test_ref() > 0.0001 \
                    or self.__convergence_test_ref_weak() > 0.0001 :

                model2 = weakCNLSDDFZG2.weakCNLSDDFZG2(
                    self.data, self.sent, self.z,self.cutactive, self.active,self.activeweak,
                    self.activeref, self.activeweakref, self.gy, self.gx, self.gb,  self.fun, self.rts, \
                    self.baseindex, self.refindex)

                model2.optimize(email, solver)
                self.alpha = model2.get_alpha()
                self.beta = model2.get_beta()
                self.gamma = model2.get_gamma()
                self.delta = model2.get_delta()


                self.__model__ = model2.__model__
                self.count += 1
                # if self.count > 20:
                #     break

            self.optimization_status = 1
            self.tt = time.time() - self.t0

        else:
            while (self.__convergence_test() > 0.0000000001) :
                while (self.__convergence_test_weak() > 0.0000000000001):
                    model2 = weakCNLSDDFZG2.weakCNLSDDFZG2(
                        self.data, self.sent, self.z, self.cutactive, self.active, self.activeweak,
                        self.activeref, self.activeweakref, self.gy, self.gx, self.gb, self.fun, self.rts, \
                        self.baseindex, self.refindex)

                    model2.optimize(email, solver)
                    self.alpha2 = self.alpha
                    self.alpha = model2.get_alpha()

                    self.beta = model2.get_beta()
                    self.gamma = model2.get_gamma()
                    self.delta = model2.get_delta()

                    self.__model__ = model2.__model__

                self.count += 1

            self.optimization_status = 1
            self.tt = time.time() - self.t0

    def __convergence_test(self):
        activetmp1 = 0.0
        if self.rts == RTS_VRS and self.fun == FUN_PROD:
        # go into the loop
            for i in self.x.index:
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in self.x.index:
                    self.active2.loc[i, j] = self.alpha.loc[i] + np.sum(np.array(self.beta.loc[i, :]) * np.array(self.x.loc[i, :])) \
                                        + np.sum(np.array(self.delta.loc[i, :]) * np.array(self.b.loc[i, :])) \
                                         - np.sum(np.array(self.gamma.loc[i, :]) * np.array(self.y.loc[i, :])) \
                                         - self.alpha.loc[j] - np.sum(np.array(self.beta.loc[j, :]) * np.array(self.x.loc[i, :])) \
                                          - np.sum(np.array(self.delta.loc[j, :]) * np.array(self.b.loc[j, :])) \
                                             + np.sum(np.array(self.gamma.loc[j, :]) * np.array(self.y.loc[i, :]))
                    if self.active2.loc[i, j] > activetmp:
                        activetmp = self.active2.loc[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in self.x.index:
                    if self.active2.loc[i, j] >= activetmp and activetmp > 0:
                        self.active.loc[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_VRS and self.fun == FUN_COST:
        # go into the loop
            for i in self.x.index:
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in self.x.index:
                    self.active2.loc[i, j] = -self.alpha.loc[i] - np.sum(np.array(self.beta.loc[i, :]) * np.array(self.x.loc[i, :])) \
                                       - np.sum(np.array(self.delta.loc[i, :]) * np.array(self.b.loc[i, :])) \
                                        + np.sum(np.array(self.gamma.loc[i, :]) * np.array(self.y.loc[i, :])) \
                                         + self.alpha.loc[j] + np.sum(np.array(self.beta.loc[j, :]) * np.array(self.x.loc[i, :])) \
                                          + np.sum(np.array(self.delta.loc[j, :]) * np.array(self.b.loc[j, :])) \
                                             - np.sum(np.array(self.gamma.loc[j, :]) * np.array(self.y.loc[i, :]))

                    if self.active2.loc[i, j] > activetmp:
                        activetmp = self.active2.loc[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in self.x.index:
                    if self.active2.loc[i, j] >= activetmp and activetmp > 0:
                        self.active.loc[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_CRS and self.fun == FUN_PROD:
        # go into the loop
            for i in self.x.index:
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in self.x.index:
                    self.active2.loc[i, j] =  np.sum(np.array(self.beta.loc[i, :]) * np.array(self.x.loc[i, :])) \
                                        + np.sum(np.array(self.delta.loc[i, :]) * np.array(self.b.loc[i, :])) \
                                         - np.sum(np.array(self.gamma.loc[i, :]) * np.array(self.y.loc[i, :])) \
                                          - np.sum(np.array(self.beta.loc[j, :]) * np.array(self.x.loc[i, :])) \
                                          - np.sum(np.array(self.delta.loc[j, :]) * np.array(self.b.loc[j, :])) \
                                             + np.sum(np.array(self.gamma.loc[j, :]) * np.array(self.y.loc[i, :]))

                    if self.active2.loc[i, j] > activetmp:
                        activetmp = self.active2.loc[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in self.x.index:
                    if self.active2.lco[i, j] >= activetmp and activetmp > 0:
                        self.active.loc[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_CRS and self.fun == FUN_COST:
        # go into the loop
            for i in self.x.index:
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in self.x.index:
                    self.active2.loc[i, j] = - np.sum(np.array(self.beta.loc[i, :]) * np.array(self.x.loc[i, :])) \
                                        - np.sum(np.array(self.delta.loc[i, :]) * np.array(self.b.loc[i, :])) \
                                         + np.sum(np.array(self.gamma.loc[i, :]) * np.array(self.y.loc[i, :])) \
                                          + np.sum(np.array(self.beta.loc[j, :]) * np.array(self.x.loc[i, :])) \
                                          + np.sum(np.array(self.delta.loc[j, :]) * np.array(self.b.loc[j, :])) \
                                         - np.sum(np.array(self.gamma.loc[j, :]) * np.array(self.y.loc[i, :]))
                    if self.active2.loc[i, j] > activetmp:
                        activetmp = self.active2.loc[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in self.x.index:
                    if self.active2.loc[i, j] >= activetmp and activetmp > 0:
                        self.active.loc[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp
        return activetmp

    def __convergence_test_weak(self):

        activetmp1 = 0.0
        if self.rts == RTS_VRS and self.fun == FUN_PROD:
        # go into the loop
            for i in self.x.index:

                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in self.x.index:
                    self.activeweak2.loc[i, j] = - self.alpha.loc[j] - np.sum(np.array(self.beta.loc[j, :]) * np.array(self.x.loc[i, :]))

                    if self.activeweak2.loc[i, j] > activetmp:
                        activetmp = self.activeweak2.loc[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in self.x.index:
                    if self.activeweak2.loc[i, j] >= activetmp and activetmp > 0:
                        self.activeweak.loc[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_VRS and self.fun == FUN_COST:
        # go into the loop
            for i in self.x.index:
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in self.x.index:
                    self.activeweak2.loc[i, j] =  self.alpha.loc[j] + np.sum(np.array(self.beta.loc[j, :]) * np.array(self.x.loc[i, :]))

                    if self.activeweak2.loc[i, j] > activetmp:
                        activetmp = self.activeweak2.loc[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in self.x.index:
                    if self.activeweak2.loc[i, j] >= activetmp and activetmp > 0:
                        self.activeweak.loc[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_CRS and self.fun == FUN_PROD:
        # go into the loop
            for i in self.x.index:
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in self.x.index:
                    self.activeweak2.loc[i, j] =  - np.sum(np.array(self.beta.loc[j, :]) * np.array(self.x.loc[i, :]))

                    if self.activeweak2.loc[i, j] > activetmp:
                        activetmp = self.activeweak2.loc[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in self.x.index:
                    if self.activeweak2.loc[i, j] >= activetmp and activetmp > 0:
                        self.activeweak.loc[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_CRS and self.fun == FUN_COST:
        # go into the loop
            for i in self.x.index:
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in self.x.index:
                    self.activeweak2.loc[i, j] =  np.sum(np.array(self.beta.loc[j, :]) * np.array(self.x.loc[i, :]))
                    if self.activeweak2.loc[i, j] > activetmp:
                        activetmp = self.activeweak2.loc[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in self.x.index:
                    if self.activeweak2.loc[i, j] >= activetmp and activetmp > 0:
                        self.activeweak.loc[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp
        return activetmp

    def __convergence_test_ref(self):
        activetmp1 = 0.0
        if self.rts == RTS_VRS and self.fun == FUN_PROD:
        # go into the loop
            for i in self.x.index:
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in self.xref.index:
                    self.activeref2.loc[i, j] = 0 \
                                         - self.alpha.loc[i] - np.sum(np.array(self.beta.loc[i, :]) * np.array(self.xref.loc[j, :])) \
                                    - np.sum(np.array(self.delta.loc[i, :]) * np.array(self.bref.loc[j, :])) \
                                     + np.sum(np.array(self.gamma.loc[i, :]) * np.array(self.yref.loc[j, :]))

                    if self.activeref2.loc[i, j] > activetmp:
                        activetmp = self.activeref2.loc[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in self.xref.index:
                    if self.activeref2.loc[i, j] >= activetmp and activetmp > 0:
                        self.activeref.loc[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_VRS and self.fun == FUN_COST:
        # go into the loop
            for i in self.x.index:
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in self.xref.index:
                    self.activeref2.loc[i, j] = 0 \
                                         + self.alpha.loc[i] +np.sum(np.array(self.beta.loc[i, :]) * np.array(self.xref.loc[j, :])) \
                                    + np.sum(np.array(self.delta.loc[i, :]) * np.array(self.bref.loc[j, :])) \
                                     - np.sum(np.array(self.gamma.loc[i, :]) * np.array(self.yref.loc[j, :]))

                    if self.activeref2.loc[i, j] > activetmp:
                        activetmp = self.activeref2.loc[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in self.xref.index:
                    if self.activeref2.loc[i, j] >= activetmp and activetmp > 0:
                        self.activeref.loc[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_CRS and self.fun == FUN_PROD:
        # go into the loop
            for i in self.x.index:
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in self.xref.index:
                    self.activeref2.loc[i, j] = 0 \
                                     - np.sum(np.array(self.beta.loc[i, :]) * np.array(self.xref.loc[j, :])) \
                                    - np.sum(np.array(self.delta.loc[i, :]) * np.array(self.bref.loc[j, :])) \
                                     + np.sum(np.array(self.gamma.loc[i, :]) * np.array(self.yref.loc[j, :]))

                    if self.activeref2.loc[i, j] > activetmp:
                        activetmp = self.activeref2.loc[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in self.xref.index:
                    if self.activeref2.lco[i, j] >= activetmp and activetmp > 0:
                        self.activeref.loc[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_CRS and self.fun == FUN_COST:
        # go into the loop
            for i in self.x.index:
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in self.xref.index:
                    self.activeref2.loc[i, j] = 0 \
                                     + np.sum(np.array(self.beta.loc[i, :]) * np.array(self.xref.loc[j, :])) \
                                    + np.sum(np.array(self.delta.loc[i, :]) * np.array(self.bref.loc[j, :])) \
                                     - np.sum(np.array(self.gamma.loc[i, :]) * np.array(self.yref.loc[j, :]))
                    if self.activeref2.loc[i, j] > activetmp:
                        activetmp = self.activeref2.loc[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in self.xref.index:
                    if self.activeref2.loc[i, j] >= activetmp and activetmp > 0:
                        self.activeref.loc[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp
        return activetmp


    def __convergence_test_ref_weak(self):
        activetmp1 = 0.0
        if self.rts == RTS_VRS and self.fun == FUN_PROD:
        # go into the loop
            for i in self.x.index:
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in self.xref.index:
                    self.activeweakref2.loc[i, j] = - self.alpha.loc[i] - np.sum(np.array(self.beta.loc[i, :]) * np.array(self.xref.loc[j, :]))

                    if self.activeweakref2.loc[i, j] > activetmp:
                        activetmp = self.activeweakref2.loc[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in self.xref.index:
                    if self.activeweakref2.loc[i, j] >= activetmp and activetmp > 0:
                        self.activeweakref.loc[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_VRS and self.fun == FUN_COST:
        # go into the loop
            for i in self.x.index:
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in self.xref.index:
                    self.activeweakref2.loc[i, j] = + self.alpha.loc[i] + np.sum(np.array(self.beta.loc[i, :]) * np.array(self.xref.loc[j, :]))

                    if self.activeweakref2.loc[i, j] > activetmp:
                        activetmp = self.activeweakref2.loc[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in self.xref.index:
                    if self.activeweakref2.loc[i, j] >= activetmp and activetmp > 0:
                        self.activeweakref.loc[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_CRS and self.fun == FUN_PROD:
        # go into the loop
            for i in self.x.index:
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in self.xref.index:
                    self.activeweakref2.loc[i, j] =  - np.sum(np.array(self.beta.loc[i, :]) * np.array(self.xref.loc[j, :]))

                    if self.activeweakref2.loc[i, j] > activetmp:
                        activetmp = self.activeweakref2.loc[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in self.xref.index:
                    if self.activeweakref2.loc[i, j] >= activetmp and activetmp > 0:
                        self.activeweakref.loc[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_CRS and self.fun == FUN_COST:
        # go into the loop
            for i in self.x.index:
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in self.xref.index:
                    self.activeweakref2.loc[i, j] =   np.sum(np.array(self.beta.loc[i, :]) * np.array(self.xref.loc[j, :]))
                    if self.activeweakref2.loc[i, j] > activetmp:
                        activetmp = self.activeweakref2.loc[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in self.xref.index:
                    if self.activeweakref2.loc[i, j] >= activetmp and activetmp > 0:
                        self.activeweakref.loc[i, j] = 1
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
        """Display residual value"""
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
        residual = list(self.__model__.epsilon[:].value)
        return np.asarray(residual)

    def get_lamda(self):
        """Return beta value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_contextual_variable(self.z)
        lamda = list(self.__model__.lamda[:].value)
        return np.asarray(lamda)

    def get_frontier(self):
        """Return estimated frontier value by array"""
        raise ValueError("DDF has no frontier.")

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
