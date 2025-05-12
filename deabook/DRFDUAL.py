from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, maximize, Constraint, Reals,PositiveReals
import numpy as np
import pandas as pd
from .constant import LEFT, RIGHT, RTS_VRS, RTS_CRS, OPT_DEFAULT, OPT_LOCAL
from .utils import tools
import ast


class DRFDUAL():

    def __init__(self, data,year,sent = "inputvar=outputvar:unoutputvar", fenmu="unoutputvar", fenzi="inputvar", \
                      side=LEFT,  rts=RTS_VRS, baseindex=None,refindex=None):
        """DRFDUAL: Dual of Directional response function

        Args:
            data (pandas.DataFrame): input pandas.
            sent (str): inputvars=outputvars[: unoutputvars]. e.g.: "K L = Y : CO2"
            gy (list): output directional vector. Defaults to [1].
            gx (list): input directional vector. Defaults to [1].
            gb (list): undesirable output directional vector. Defaults to None.
            rts (String, optional): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # Initialize DEA model
        self.data=data
        self.year = year
        self.tlt=pd.Series(self.year).drop_duplicates().sort_values()
        self.outputvars,self.inputvars,self.unoutputvars ,self.obj_coeflt, self.rule4_coeflt,self.neg_obj \
            = tools.assert_valid_yxb_drf(sent,fenmu,fenzi)    ## 判断分母是x，b or y，是x，b的，目标要加负号

        self.rts = rts
        self.side = side
        self.baseindex = baseindex
        if type(baseindex) != type(None):
            self.varname1=self.baseindex.split('=')[0]
            print(self.baseindex)
            self.varvalue1=ast.literal_eval(self.baseindex.split('=')[1])
            self.y, self.x, self.b = self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.outputvars
                                        ], self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.inputvars
                                        ], self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.unoutputvars
                                        ]

        else:

            self.y, self.x, self.b = self.data.loc[:, self.outputvars
                                        ], self.data.loc[:, self.inputvars
                                        ], self.data.loc[:, self.unoutputvars
                                        ]


        # print(type(self.varname1),self.varvalue1,self.x,)
        self.refindex = refindex
        if type(refindex) != type(None):
            self.varname=self.refindex.split('=')[0]
            self.varvalue=ast.literal_eval(self.refindex.split('=')[1])

            self.yref, self.xref, self.bref = self.data.loc[self.data[self.varname].isin(self.varvalue), self.outputvars
                                                ], self.data.loc[self.data[self.varname].isin(self.varvalue), self.inputvars
                                                ], self.data.loc[self.data[self.varname].isin(self.varvalue), self.unoutputvars
                                                ]
        else:
            self.yref, self.xref, self.bref = self.data.loc[:, self.outputvars
                                        ], self.data.loc[:, self.inputvars
                                        ], self.data.loc[:, self.unoutputvars
                                        ]



        self.xcol = self.x.columns
        self.ycol = self.y.columns
        self.bcol = self.b.columns

        self.xobj_coef = self.obj_coeflt['xobj_coef']
        self.yobj_coef = self.obj_coeflt['yobj_coef']
        self.bobj_coef = self.obj_coeflt['bobj_coef']

        self.xrule4_coef = self.rule4_coeflt['xrule4_coef']
        self.yrule4_coef = self.rule4_coeflt['yrule4_coef']
        self.brule4_coef = self.rule4_coeflt['brule4_coef']

        print("xcol,ycol,bcol are:",self.x.columns,self.y.columns,self.b.columns)

        print("fenmu is:",self.xobj_coef,self.yobj_coef,self.bobj_coef)
        print("fenzi is:",self.xrule4_coef,self.yrule4_coef,self.brule4_coef)
        print("neg is:",self.neg_obj)

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
            xlb,xub={},{}
            for ii,j in enumerate(self.xobj_coef):
                if j>0:
                    xlb[ii]=None
                    xub[ii]=None
                else:
                    xlb[ii]=0
                    xub[ii]=None
            def xfb(__model__, i):
                return (xlb[i], xub[i])

            ylb,yub={},{}
            for ii,j in enumerate(self.yobj_coef):
                if j>0:
                    ylb[ii]=None
                    yub[ii]=None
                else:
                    ylb[ii]=0
                    yub[ii]=None
            def yfb(__model__, i):
                return (ylb[i], yub[i])

            blb,bub={},{}
            for ii,j in enumerate(self.bobj_coef):
                if j>0:
                    blb[ii]=None
                    bub[ii]=None
                else:
                    blb[ii]=0
                    bub[ii]=None
            def bfb(__model__, i):
                return (blb[i], bub[i])

            self.__model__.spx = Var(self.__model__.K,bounds=xfb, within=Reals,doc='shadow price of x')
            self.__model__.spy = Var(self.__model__.L,bounds=yfb,within=Reals, doc='shadow price of y')
            self.__model__.spb = Var(self.__model__.J,bounds=bfb ,within=Reals, doc='shadow price of b')



            if self.rts == RTS_VRS:
                self.__model__.spalpha = Var(Set(initialize=range(1)),  within=Reals,doc='shadow price of 1')

            # Setup the objective function and constraints
            if (self.side== RIGHT) :
                if self.neg_obj:
                    self.__model__.objective = Objective(rule=self.__objective_rule(), sense=maximize, doc='objective function')
                else:
                    self.__model__.objective = Objective(rule=self.__objective_rule(), sense=minimize, doc='objective function')

            elif (self.side== LEFT) :
                if self.neg_obj:
                    self.__model__.objective = Objective(rule=self.__objective_rule(), sense=minimize, doc='objective function')
                else:
                    self.__model__.objective = Objective(rule=self.__objective_rule(), sense=maximize, doc='objective function')

            self.__model__.first = Constraint(self.__model__.I2,  rule=self.__first_rule(), doc='first constraint')
            self.__model__.second = Constraint(self.__model__.I2,  rule=self.__second_rule(), doc='second constraint')
            self.__model__.third = Constraint(                    rule=self.__third_rule(), doc='inefficiency=0 constraint')
            self.__model__.forth = Constraint(                    rule=self.__forth_rule(), doc='fenmu=1 constraint')


            self.__modeldict[i] = self.__model__

        # Optimize model
    def __objective_rule(self):
        """Return the proper objective function"""
        if self.neg_obj:
            def objective_rule(model):

                    return -1 * (sum(model.spx[k]*self.xobj_coef[k] for k in model.K
                        ) - sum(model.spy[l]*self.yobj_coef[l] for l in model.L
                        ) + sum(model.spb[j]*self.bobj_coef[j] for j in model.J))
            return objective_rule
        else:
            def objective_rule(model):

                    return sum(model.spx[k]*self.xobj_coef[k] for k in model.K
                        ) - sum(model.spy[l]*self.yobj_coef[l] for l in model.L
                        ) + sum(model.spb[j]*self.bobj_coef[j] for j in model.J)
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

            return sum(model.spx[k] * self.x.loc[self.I0, self.xcol[k]] for k in model.K
                ) - sum(model.spy[l] * self.y.loc[self.I0, self.ycol[l]] for l in model.L
                ) + sum(model.spb[j] * self.b.loc[self.I0, self.bcol[j]] for j in model.J
                ) + (model.spalpha[0]*1 if self.rts == RTS_VRS else 0)  == 0

        return third_rule

    def __forth_rule(self):
        """Return the proper third constraint"""
        def forth_rule(model):
            return  sum(self.xrule4_coef[k]*model.spx[k] for k in model.K
                  )+sum(self.yrule4_coef[l]*model.spy[l] for l in model.L
                  )+sum(self.brule4_coef[j]*model.spb[j] for j in model.J) == 1

        return forth_rule

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

