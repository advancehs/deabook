from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize,maximize, Constraint,Reals
from pyomo.core.expr.numvalue import NumericValue
import pandas as pd
import numpy as np

from . import  weakCNLS
from .constant import CET_ADDI, FUN_COST, FUN_PROD, RTS_VRS, RTS_CRS,OPT_DEFAULT, OPT_LOCAL,LEFT,RIGHT
from .utils import tools
import ast


class weakCNLSDRF():

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

                                       ## I 是 被评价决策单元的数量

        self.__model__ = ConcreteModel()
        # Initialize sets
        self.__model__.I = Set(initialize=self.x.index)  ## I 是 被评价决策单元的数量

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

        self.__model__.beta = Var(self.__model__.I,self.__model__.K,bounds=(0,None), within=Reals,doc='shadow price of x')
        self.__model__.gamma = Var(self.__model__.I,self.__model__.L,bounds=(0,None),within=Reals, doc='shadow price of y')
        self.__model__.delta = Var(self.__model__.I,self.__model__.J,bounds=(0,None) ,within=Reals, doc='shadow price of b')



        if self.rts == RTS_VRS:
            self.__model__.alpha = Var(self.__model__.I,  within=Reals,doc='shadow price of 1')

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

        self.__model__.frontier_rule = Constraint(self.__model__.I,  rule=self.__frontier_rule(), doc='frontier constraint')
        self.__model__.afriat_rule = Constraint(self.__model__.I, self.__model__.I2,
                                          rule=self.__afriat_rule(), doc='afriat_rule constraint')
        self.__model__.disposability_rule = Constraint( self.__model__.I, self.__model__.I2,
                                          rule=self.__disposability_rule(), doc='disposability_rule constraint')
        self.__model__.translation_rule = Constraint( self.__model__.I, \
                                                      rule=self.__translation_rule(), doc='translation_rule constraint')



        # Optimize model
    def __objective_rule(self):
        """Return the proper objective function"""
        if self.neg_obj:
            def objective_rule(model):

                    return -1 * (sum( sum(model.beta[i,k]*self.xobj_coef[k] for k in model.K
                                    ) - sum(model.gamma[i,l]*self.yobj_coef[l] for l in model.L
                                    ) + sum(model.delta[i,j]*self.bobj_coef[j] for j in model.J) for i in model.I))
            return objective_rule
        else:
            def objective_rule(model):

                    return sum( sum(model.beta[i,k]*self.xobj_coef[k] for k in model.K
                                    ) - sum(model.gamma[i,l]*self.yobj_coef[l] for l in model.L
                                    ) + sum(model.delta[i,j]*self.bobj_coef[j] for j in model.J) for i in model.I)
            return objective_rule

    def __frontier_rule(self):
        """Return the proper first constraint"""
        def frontier_rule(model, i):

            return sum(model.beta[i,k] * self.xref.loc[i,self.xcol[k]] for k in model.K
                ) - sum(model.gamma[i,l] * self.yref.loc[i,self.ycol[l]] for l in model.L
                ) + sum(model.delta[i,j] * self.bref.loc[i,self.bcol[j]] for j in model.J
                ) + (model.alpha[0]*1 if self.rts == RTS_VRS else 0)   ==0
        return frontier_rule

    def __afriat_rule(self):
        """Return the proper second constraint"""

        __operator = NumericValue.__le__

        if self.rts == RTS_VRS:
            def afriat_rule(model, i, h):
                if i == h:
                    return Constraint.Skip
                return __operator(model.alpha[i] \
                                  + sum(model.beta[i, k] * self.x.loc[i,self.xcol[k]]  for k in model.K) \
                                  + sum(model.delta[i, j] * self.b.loc[i,self.bcol[j]] for j in model.J) \
                                  - sum(model.gamma[i, l] * self.y.loc[i,self.ycol[l]] for l in model.L),
                                  model.alpha[h]
                                  + sum(model.beta[h, k] * self.x.loc[i,self.xcol[k]]  for k in model.K) \
                                  + sum(model.delta[h, j] * self.b.loc[i,self.bcol[j]] for j in model.J) \
                                  - sum(model.gamma[h, l] * self.y.loc[i,self.ycol[l]] for l in model.L))
            return afriat_rule

        elif self.rts == RTS_CRS:
            def afriat_rule(model, i, h):
                if i == h:
                    return Constraint.Skip
                return __operator(  sum(model.beta[i, k] * self.x.loc[i,self.xcol[k]]  for k in model.K) \
                                  + sum(model.delta[i, j] * self.b.loc[i,self.bcol[j]] for j in model.J) \
                                  - sum(model.gamma[i, l] * self.y.loc[i,self.ycol[l]] for l in model.L),
                                    sum(model.beta[h, k] * self.x.loc[i,self.xcol[k]] for k in model.K) \
                                  + sum(model.delta[h, j] * self.b.loc[i,self.bcol[j]] for j in model.J) \
                                  - sum(model.gamma[h, l] * self.y.loc[i,self.ycol[l]] for l in model.L))
            return afriat_rule

        raise ValueError("Undefined model parameters.")
    def __disposability_rule(self):
        """Return the proper third constraint"""
        if self.rts == RTS_VRS:
            def disposability_rule(model, i, h):
                if i == h:
                    return Constraint.Skip
                return model.alpha[i] + sum(model.beta[i, k] * self.x.loc[h,self.xcol[k]] for k in model.K) >= 0
            return disposability_rule

        elif self.rts == RTS_CRS:
            def disposability_rule(model, i, h):
                if i == h:
                    return Constraint.Skip
                return sum(model.beta[i, k] * self.x.loc[h,self.xcol[k]] for k in model.K) >= 0
            return disposability_rule

        raise ValueError("Undefined model parameters.")

    def __translation_rule(self):
        """Return the proper third constraint"""
        def translation_rule(model, i):
            return sum(model.beta[i, k] * self.x.loc[i,self.xcol[k]] for k in model.K) \
                + sum(model.gamma[i, l] * self.y.loc[i,self.ycol[l]] for l in model.L) \
                + sum(model.delta[i, j] * self.b.loc[i,self.bcol[j]] for j in model.J) == 1

        return translation_rule

    def optimize(self,  solver=OPT_DEFAULT):
        """Optimize the function by requested method

        Args:
            solver (string): The solver chosen for optimization. It will optimize with default solver if OPT_DEFAULT is given.
        """
        # TODO(error/warning handling): Check problem status after optimization

        data2,obj,spx,spy,spb,  = pd.DataFrame,{},{},{},{},

        _, _ = tools.optimize_model3(self.__model__, solver)

        obj= self.__model__.objective()
        gamma = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.gamma),
                                                           list(self.__model__.gamma[:, :].value))])
        gamma = pd.DataFrame(gamma, columns=['Name', 'Key', 'Value'])
        gamma = gamma.pivot(index='Name', columns='Key', values='Value')


        delta = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.delta),
                                                           list(self.__model__.delta[:, :].value))])
        delta = pd.DataFrame(delta, columns=['Name', 'Key', 'Value'])
        delta = delta.pivot(index='Name', columns='Key', values='Value')

        beta = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.beta),
                                                           list(self.__model__.beta[:, :].value))])
        beta = pd.DataFrame(beta, columns=['Name', 'Key', 'Value'])
        beta = beta.pivot(index='Name', columns='Key', values='Value')
        #
        print(obj)
        # obj = pd.DataFrame(obj,index=["obj"]).T
        beta.columns = beta.columns.map(lamda x : "Input"+ str(x)+"'s shadow price" )
        gamma.columns = gamma.columns.map(lamda y : "Output"+ str(y)+"'s shadow price" )
        delta.columns = delta.columns.map(lamda b : "Undesirable"+ str(b)+"'s shadow price" )

        # spy = pd.DataFrame(spy).T
        # spy.columns = spy.columns.map(lamda y : "Output"+ str(y)+"'s shadow price" )
        # spb = pd.DataFrame(spb).T
        # spb.columns = spb.columns.map(lamda b : "Undesirable Output"+ str(b)+"'s shadow price" )
        sp=pd.concat([beta,gamma],axis=1)
        sp=pd.concat([sp,delta],axis=1)
        # data3 = pd.concat([obj,sp],axis=1)
        return sp


    def info(self, dmu = "all"):
        """Show the information of the lp model

        Args:
            dmu (string): The solver chosen for optimization. Default is "all".
        """

        print(self.__model__.delta.pprint())


