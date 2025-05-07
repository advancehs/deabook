from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, maximize, Constraint, Reals, PositiveReals
import numpy as np
import pandas as pd
from .constant import CET_ADDI, ORIENT_IO, ORIENT_OO, RTS_VRS, RTS_CRS, OPT_DEFAULT, OPT_LOCAL
from .utils import tools
import ast


class DDFDUAL():

    def __init__(self, data,sent = "inputvar=outputvar:unoutputvar",  gy=[1], gx=[1], gb=[1], rts=RTS_VRS, baseindex=None,refindex=None):
        """DDFDUAL: Dual of Directional distance function

        Args:
            data (pandas.DataFrame): input pandas.
            sent (str): inputvars=outputvars[: unoutputvars]. e.g.: "K L = Y : CO2"
            gy (pandas.Series): output directional vector. Defaults to [1].
            gx (pandas.Series): input directional vector. Defaults to [1].
            gb (pandas.Series): undesirable output directional vector. Defaults to None.
            rts (String, optional): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # Initialize DEA model
        self.outputvars, self.inputvars, self.unoutputvars,  self.gy, self.gx, self.gb \
            = tools.assert_valid_yxb(sent, gy, gx, gb)
        self.y, self.x, self.b, self.yref, self.xref, self.bref, \
            = tools.assert_valid_yxb2(baseindex, refindex, data, self.outputvars, \
                                       self.inputvars, self.unoutputvars,)

        self.rts = rts
        self.xcol = self.x.columns
        self.ycol = self.y.columns
        self.bcol = self.b.columns

        print("xcol,ycol,bcol are:",self.x.columns,self.y.columns,self.b.columns)

        print("gx,gy,gb are:",self.gx,self.gy,self.gb)

        self.I = self.x.index          ## I 是 被评价决策单元的索引
        self.__modeldict = {}
        for i in self.I:
            # print(i)
            self.I0 = i                                                 ## I 是 被评价决策单元的数量

            self.__model__ = ConcreteModel()
            # Initialize sets
            self.__model__.I2 = Set(initialize=self.xref.index)      ## I2 是 参考决策单元的数量
            self.__model__.K = Set(initialize=range(len(self.x.iloc[0])))          ## K 是投入个数
            self.__model__.L = Set(initialize=range(len(self.y.iloc[0])))           ## L 是产出个数 被评价单元和参考单元的K，L一样
            self.__model__.J = Set(initialize=range(len(self.b.iloc[0])))   ## B 是 非期望产出个数


            # Initialize variable

            self.__model__.spx = Var(self.__model__.K,initialize=1,bounds=(0.0, None), within=Reals,doc='shadow price of x')
            self.__model__.spy = Var(self.__model__.L, initialize=1,bounds=(0.0, None),within=Reals, doc='shadow price of y')

            self.__model__.spb = Var(self.__model__.J,bounds=(1e-6, None),within=PositiveReals, doc='shadow price of b')
            if self.rts == RTS_VRS:
                self.__model__.spalpha = Var(Set(initialize=range(1)),  within=Reals,doc='shadow price of 1')

            # Setup the objective function and constraints
            self.__model__.objective = Objective(rule=self.__objective_rule(), sense=minimize, doc='objective function')
            self.__model__.first = Constraint(self.__model__.I2,  rule=self.__first_rule(), doc='first constraint')
            self.__model__.second = Constraint(self.__model__.I2,  rule=self.__second_rule(), doc='second constraint')
            self.__model__.third = Constraint(                    rule=self.__third_rule(), doc='third constraint')


            self.__modeldict[i] = self.__model__

        # Optimize model
    def __objective_rule(self):
        """Return the proper objective function"""
        def objective_rule(model):

                return sum(model.spx[k]*self.x.loc[self.I0,self.xcol[k]] for k in model.K
                    ) - sum(model.spy[l]*self.y.loc[self.I0,self.ycol[l]] for l in model.L
                    ) + sum(model.spb[j]*self.b.loc[self.I0,self.bcol[j]] for j in model.J
                    ) + (model.spalpha[0]*1 if self.rts == RTS_VRS else 0)
        return objective_rule

    def __first_rule(self):
        """Return the proper first constraint"""
        def first_rule(model, i2):

            return sum(model.spx[k] * self.xref.loc[i2,self.xcol[k]] for k in model.K
                ) - sum(model.spy[l] * self.yref.loc[i2,self.ycol[l]] for l in model.L
                ) + sum(model.spb[j] * self.bref.loc[i2,self.bcol[j]] for j in model.J
                ) + (model.spalpha[0]*1 if self.rts == RTS_VRS else 0)   >=0
        return first_rule

    def __second_rule(self):
        """Return the proper second constraint"""
        def second_rule(model, i2):

            return sum(model.spx[k] * self.xref.loc[i2,self.xcol[k]] for k in model.K
                ) + (model.spalpha[0]*1 if self.rts == RTS_VRS else 0)   >=0
        return second_rule

    def __third_rule(self):
        """Return the proper third constraint"""
        def third_rule(model):
            return sum(model.spx[ k] * self.gx.loc[self.I0,self.xcol[k]] for k in model.K) \
                + sum(model.spy[ l] * self.gy.loc[self.I0,self.ycol[l]] for l in model.L) \
                + sum(model.spb[ j] * self.gb.loc[self.I0,self.bcol[j]] for j in model.J) == 1

        return third_rule


    def optimize(self,  solver=OPT_DEFAULT):
        """Optimize the function by requested method

        Args:
            solver (string): The solver chosen for optimization. It will optimize with default solver if OPT_DEFAULT is given.
        """
        # TODO(error/warning handling): Check problem status after optimization

        data2,obj,spx,spy,spb,  = pd.DataFrame,{},{},{},{},
        for ind, problem in self.__modeldict.items():
            _, _ = tools.optimize_model2(problem, ind, solver)

            obj[ind]= problem.objective()
            spx[ind]= np.asarray(list(problem.spx[:].value))
            spy[ind]= np.asarray(list(problem.spy[:].value))
            spb[ind]= np.asarray(list(problem.spb[:].value))

        obj = pd.DataFrame(obj,index=["obj"]).T
        spx = pd.DataFrame(spx).T
        spx.columns = spx.columns.map(lambda x : "Input"+ str(x)+"'s shadow price" )
        spy = pd.DataFrame(spy).T
        spy.columns = spy.columns.map(lambda y : "Output"+ str(y)+"'s shadow price" )
        spb = pd.DataFrame(spb).T
        spb.columns = spb.columns.map(lambda b : "Undesirable Output"+ str(b)+"'s shadow price" )
        sp=pd.concat([spx,spy],axis=1)
        sp=pd.concat([sp,spb],axis=1)
        data3 = pd.concat([obj,sp],axis=1)
        return data3


    def info(self, dmu = "all"):
        """Show the infomation of the lp model

        Args:
            dmu (string): The solver chosen for optimization. Default is "all".
        """
        if dmu =="all":
            for ind, problem in list(self.__modeldict.items()):
                print(ind,"\n",problem.pprint())

        print(self.__modeldict[int(dmu)].pprint())

