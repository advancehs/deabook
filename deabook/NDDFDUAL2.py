from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, maximize, Constraint, Reals, PositiveReals
import numpy as np
import pandas as pd
from .constant import CET_ADDI, ORIENT_IO, ORIENT_OO, RTS_VRS, RTS_CRS, OPT_DEFAULT, OPT_LOCAL
from .utils import tools
import ast


class NDDFDUAL2():

    def __init__(self, data,sent = "inputvar=outputvar",  gy=[1], gx=[1], gb=None, weight =None,\
                 rts=RTS_VRS, baseindex=None,refindex=None):
        """DEA: Directional distance function

        Args:
            data (pandas.DataFrame): input pandas.
            sent (str): inputvars=outputvars[: unoutputvars]. e.g.: "K L = Y : CO2"
            gy (list, optional): output directional vector. Defaults to [1].
            gx (list, optional): input directional vector. Defaults to [1].
            gb (list, optional): undesirable output directional vector. Defaults to None.
            weight(list, optional): weght matrix
            rts (String): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale)
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

        if type(weight) != type(None):
            self.weight= weight
        else:
            self.weight=[]
            if type(self.b) != type(None):
                fenmu = 1*int(self.gx.iloc[0,0]!=0) + 1*int(self.gy.iloc[0,0]!=0) + 1*int(self.gb.iloc[0,0]!=0)
                print(fenmu)
                for _ in range(len(self.x.iloc[0])):
                    self.weight.append(1/fenmu/len(self.x.iloc[0]))
                for _ in range(len(self.y.iloc[0])):
                    self.weight.append(1/fenmu/len(self.y.iloc[0]))
                for _ in range(len(self.b.iloc[0])):
                    self.weight.append(1/fenmu/len(self.b.iloc[0]))
            else:
                fenmu =1*int(self.gx.iloc[0,0]!=0) + 1*int(self.gy.iloc[0,0]!=0)

                for _ in range(len(self.x.iloc[0])):
                    self.weight.append(1/fenmu/len(self.x.iloc[0]))
                for _ in range(len(self.y.iloc[0])):
                    self.weight.append(1/fenmu/len(self.y.iloc[0]))

        self.iweight = self.weight[0:len(self.x.iloc[0])]
        self.oweight = self.weight[len(self.x.iloc[0]):len(self.x.iloc[0])+len(self.y.iloc[0])]
        if type(self.b) != type(None):
            self.bweight = self.weight[len(self.x.iloc[0])+len(self.y.iloc[0]):len(self.x.iloc[0])+len(self.y.iloc[0])+len(self.b.iloc[0])]


        self.xcol = self.x.columns
        self.ycol = self.y.columns
        self.bcol = self.b.columns if type(self.b) != type(None) else None

        print(self.iweight,self.oweight,self.bweight)
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
            if type(self.b) != type(None):
                self.__model__.J = Set(initialize=range(len(self.b.iloc[0])))   ## B 是 非期望产出个数

            # Initialize variable
            self.__model__.px = Var(self.__model__.K,initialize=1,bounds=(0.0, None), within=Reals,doc='shadow price of x')
            self.__model__.py = Var(self.__model__.L, initialize=1,bounds=(0.0, None),within=Reals, doc='shadow price of y')
            if type(self.b) != type(None):
                self.__model__.pb = Var(self.__model__.J,bounds=(0.0, None),within=Reals, doc='shadow price of b')
            if self.rts == RTS_VRS:
                self.__model__.pomega = Var(Set(initialize=range(1)),  within=Reals,doc='shadow price of 1')

            # Setup the objective function and constraints
            self.__model__.objective = Objective(rule=self.__objective_rule(), sense=minimize, doc='objective function')
            self.__model__.first = Constraint(self.__model__.I2,  rule=self.__first_rule(), doc='first constraint')
            self.__model__.second = Constraint(self.__model__.K,  rule=self.__second_rule(), doc='second constraint')
            self.__model__.third = Constraint(self.__model__.L,  rule=self.__third_rule(), doc='third constraint')

            if type(self.b) != type(None):
                self.__model__.forth = Constraint(self.__model__.J,  rule=self.__forth_rule(), doc='forth constraint')

            self.__modeldict[i] = self.__model__

        # Optimize model
    def __objective_rule(self):
        """Return the proper objective function"""
        def objective_rule(model):
            if type(self.b) != type(None):
                return sum(model.px[k]*self.x.loc[self.I0,self.xcol[k]] for k in model.K
                    ) - sum(model.py[l]*self.y.loc[self.I0,self.ycol[l]] for l in model.L
                    ) + sum(model.pb[j]*self.b.loc[self.I0,self.bcol[j]] for j in model.J
                    ) + (model.pomega[0]*1 if self.rts == RTS_VRS else 0)
            else:
                return sum(model.px[k]*self.x.loc[self.I0,self.xcol[k]] for k in model.K
                    ) - sum(model.py[l]*self.y.loc[self.I0,self.ycol[l]] for l in model.L
                    ) + (model.pomega[0]*1 if self.rts == RTS_VRS else 0 )
        return objective_rule

    def __first_rule(self):
        """Return the proper first constraint"""
        def first_rule(model, i2):
            if type(self.b) != type(None):
                return sum(model.px[k] * self.xref.loc[i2,self.xcol[k]] for k in model.K
                    ) - sum(model.py[l] * self.yref.loc[i2,self.ycol[l]] for l in model.L
                    ) + sum(model.pb[j] * self.bref.loc[i2,self.bcol[j]] for j in model.J
                    ) + (model.pomega[0]*1 if self.rts == RTS_VRS else 0)   >=0
            else:
                return sum(model.px[k] *self.xref.loc[i2,self.xcol[k]]   for k in model.K
                    ) - sum(model.py[l] *self.yref.loc[i2,self.ycol[l]]  for l in model.L
                    ) + (model.pomega[0]*1 if self.rts == RTS_VRS else 0 )  >=0
        return first_rule

    def __second_rule(self):
        """Return the proper second constraint"""
        def second_rule(model, k):
            if self.gx.loc[self.I0,self.xcol[k]]==0:
                return Constraint.Skip
            return  self.gx.loc[self.I0,self.xcol[k]] * model.px[k] >= self.iweight[k]
        return second_rule

    def __third_rule(self):
        """Return the proper third constraint"""
        def third_rule(model, l):
            if self.gy.loc[self.I0,self.ycol[l]]==0:
                return Constraint.Skip
            return self.gy.loc[self.I0,self.ycol[l]]*model.py[l] >= self.oweight[l]

        return third_rule

    def __forth_rule(self):
        """Return the proper forth constraint"""
        def forth_rule(model, j):
            if self.gb.loc[self.I0,self.bcol[j]]==0:
                return Constraint.Skip
            return  self.gb.loc[self.I0,self.bcol[j]]*model.pb[j] >= self.bweight[j]
        return forth_rule

    def optimize(self,  solver=OPT_DEFAULT):
        """Optimize the function by requested method

        Args:
            solver (string): The solver chosen for optimization. It will optimize with default solver if OPT_DEFAULT is given.
        """
        # TODO(error/warning handling): Check problem status after optimization

        data2,obj,px,py,pb,  = pd.DataFrame,{},{},{},{},
        for ind, problem in self.__modeldict.items():
            _, _ = tools.optimize_model2(problem, ind, solver)
            if type(self.b) != type(None):
                obj[ind]= problem.objective()
                px[ind]= np.asarray(list(problem.px[:].value))
                py[ind]= np.asarray(list(problem.py[:].value))
                pb[ind]= np.asarray(list(problem.pb[:].value))
                # pomega[ind]=
            else:
                obj[ind]= problem.objective()
                px[ind]= np.asarray(list(problem.px[:].value))
                py[ind]= np.asarray(list(problem.py[:].value))
                # pomega[ind]=
        obj = pd.DataFrame(obj,index=["obj"]).T
        px = pd.DataFrame(px).T
        px.columns = px.columns.map(lamda x : "Input"+ str(x)+"'s shadow price" )
        py = pd.DataFrame(py).T
        py.columns = py.columns.map(lamda y : "Output"+ str(y)+"'s shadow price" )
        pb = pd.DataFrame(pb).T
        pb.columns = pb.columns.map(lamda b : "Undesirable Output"+ str(b)+"'s shadow price" )
        p=pd.concat([px,py],axis=1)
        p=pd.concat([p,pb],axis=1)
        # data3 = pd.concat([data2,obj],axis=1)
        data3 = pd.concat([obj,p],axis=1)
        return data3


    def info(self, dmu = "all"):
        """Show the information of the lp model

        Args:
            dmu (string): The solver chosen for optimization. Default is "all".
        """
        if dmu =="all":
            for ind, problem in list(self.__modeldict.items()):
                print(ind,"\n",problem.pprint())

        print(self.__modeldict[int(dmu)].pprint())

