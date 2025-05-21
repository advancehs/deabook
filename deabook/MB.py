"""Main module."""
# import dependencies

from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, maximize, Constraint, Reals
import numpy as np
import pandas as pd
from .constant import CET_ADDI, RTS_VRS1, RTS_CRS,OPT_DEFAULT, OPT_LOCAL
from .utils import tools,MB_
from .DEA import DEA
import ast


class MB():
    def __init__(self, data,sent = "inputvar_np + inputvar_p =outputvar_np + outputvar_p:unoutputvar",  \
                 sx=[[1,1,1],[1,1,1]], sy=[[1],[1]], level=5 ,rts=RTS_VRS1, baseindex=None,refindex=None):
        """DEA: Directional distance function

        Args:
            data (pandas.DataFrame): input pandas.
            sent (str): inputvars=outputvars[: unoutputvars]. e.g.: "K L+ E = Y : CO2" 使用“+”区分污染和废污染投入（或产出）
            sx (list): 投入包含污染物质系数. Defaults to [[1,1,1],[1,1,1]].
            sy (list, optional): 期望产出包含污染物质系数. Defaults to [[1],[1]].
            level(int, optional): 返回求的层级。1：只变量化b；2：再变量化x_p；3：再变量化x_np；\
                                                4：再变量化y_p（没有则变量化y_np）；5：再变量化y_np
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # Initialize DEA model
        self.data=data
        # self.year = year
        self.sent = sent
        # self.tlt=pd.Series(self.year).drop_duplicates().sort_values()
        self.rts = rts
        self.baseindex = baseindex
        self.refindex = refindex

        self.inputvars_np, self.inputvars_p,self.outputvars_np,\
            self.outputvars_p,self.unoutputvars, self.sx, self.sy,self.level = tools.split_MB(self.sent, sx, sy,level)

    def optimize(self, solver=OPT_DEFAULT,dmu="1"):

        if type(self.outputvars_np) != type(None):
            if type(self.outputvars_p) != type(None):

                model1111 = MB_.MB1111(  ## K L + E = Y + Y2 : CO2
                               self.data, self.inputvars_np, self.inputvars_p,self.outputvars_np,self.outputvars_p,\
                               self.unoutputvars, self.sx, self.sy,self.rts,self.level,self.baseindex,self.refindex)
                data3 = model1111.optimize(solver=solver)
                info = model1111.info(dmu=dmu)
            elif type(self.outputvars_p) == type(None):
                model1110 = MB_.MB1110(  ## K L + E = Y  : CO2
                               self.data, self.inputvars_np, self.inputvars_p,self.outputvars_np,\
                                self.unoutputvars, self.sx, self.sy,self.rts,self.level,self.baseindex,self.refindex)
                data3 = model1110.optimize(solver=solver)
                info = model1110.info(dmu=dmu)
        elif type(self.outputvars_np) == type(None):
            if type(self.outputvars_p) != type(None):
                model1101 = MB_.MB1101(  ## K L + E = + Y2 : CO2
                               self.data, self.inputvars_np, self.inputvars_p,self.outputvars_p,\
                                self.unoutputvars, self.sx, self.sy,self.rts,self.level,self.baseindex,self.refindex)
                data3 = model1101.optimize(solver=solver)
                info = model1101.info(dmu=dmu)
            elif type(self.outputvars_p) == type(None):
                model1100 = MB_.MB1100(  ## K L + E =   : CO2
                               self.data, self.inputvars_np, self.inputvars_p,\
                                self.unoutputvars, self.sx, self.sy,self.rts,self.level,self.baseindex,self.refindex)
                data3 = model1100.optimize(solver=solver)
                info = model1100.info(dmu=dmu)
        return data3,info




    def info(self, dmu="all"):
        if type(self.outputvars_np) != type(None):
            if type(self.outputvars_p) != type(None):
                model1111 = MB_.MB1111(
                               self.data, self.inputvars_np, self.inputvars_p,self.outputvars_np,\
                               self.outputvars_p,self.unoutputvars, self.sx, self.sy, \
                                self.rts,self.level,self.baseindex,self.refindex)
                return model1111.info(dmu=dmu)
            elif type(self.outputvars_p) == type(None):
                model1110 = MB_.MB1110(  ## K L + E = Y  : CO2
                               self.data, self.inputvars_np, self.inputvars_p,self.outputvars_np,\
                                self.unoutputvars, self.sx, self.sy,self.rts,self.level,self.baseindex,self.refindex)
                return model1110.info(dmu=dmu)
        elif type(self.outputvars_np) == type(None):
            if type(self.outputvars_p) != type(None):
                model1101 = MB_.MB1101(  ## K L + E = + Y2 : CO2
                               self.data, self.inputvars_np, self.inputvars_p,self.outputvars_p,\
                                self.unoutputvars, self.sx, self.sy,self.rts,self.level,self.baseindex,self.refindex)
                return model1101.info(dmu=dmu)
            elif type(self.outputvars_p) == type(None):
                model1100 = MB_.MB1100(  ## K L + E =   : CO2
                               self.data, self.inputvars_np, self.inputvars_p,\
                                self.unoutputvars, self.sx, self.sy,self.rts,self.level,self.baseindex,self.refindex)
                return model1100.info(dmu=dmu)



class MBx(MB):
    def __init__(self, data,year,sent = "inputvar=outputvar",  sx=[[1,1,1],[1,1,1]], sy=[[1],[1]], rts=RTS_VRS1, baseindex=None,refindex=None):
        """DEA: Directional distance function

        Args:
            data (pandas.DataFrame): input pandas.
            sent (str): inputvars=outputvars[: unoutputvars]. e.g.: "K L = Y : CO2"
            sx (list): 投入包含污染物质系数. Defaults to [[1,1,1],[1,1,1]].
            sy (list, optional): 期望产出包含污染物质系数. Defaults to [[1],[1]].
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # Initialize DEA model
        self.data=data
        self.year = year
        self.sent = sent
        self.tlt=pd.Series(self.year).drop_duplicates().sort_values()
        self.inputvars = self.sent.split('=')[0].strip(' ').split(' ')
        try:
            self.outputvars = self.sent.split('=')[1]   .split(':')[0].strip(' ').split(' ')
            self.unoutputvars = self.sent.split('=')[1]   .split(':')[1].strip(' ').split(' ')
        except:
            self.outputvars = self.sent.split('=')[1]    .strip(' ').split(' ')
            self.unoutputvars=None
        self.sx, self.sy = sx, sy
        self.rts = rts
        # print(self.sx, self.sy)

        self.baseindex = baseindex
        if type(baseindex) != type(None):
            self.varname1=self.baseindex.split('=')[0]
            # print(self.baseindex)
            self.varvalue1=ast.literal_eval(self.baseindex.split('=')[1])
            self.y, self.x, self.b = self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.outputvars
                                        ], self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.inputvars
                                        ], self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.unoutputvars
                                        ]if type(self.unoutputvars) != type(None) else None

        else:

            self.y, self.x, self.b = self.data.loc[:, self.outputvars
                                        ], self.data.loc[:, self.inputvars
                                        ], self.data.loc[:, self.unoutputvars
                                        ] if type(self.unoutputvars) != type(None) else None


        # print(self.b)
        self.refindex = refindex
        if type(refindex) != type(None):
            self.varname=self.refindex.split('=')[0]
            self.varvalue=ast.literal_eval(self.refindex.split('=')[1])

            self.yref, self.xref, self.bref = self.data.loc[self.data[self.varname].isin(self.varvalue), self.outputvars
                                                ], self.data.loc[self.data[self.varname].isin(self.varvalue), self.inputvars
                                                ], self.data.loc[self.data[self.varname].isin(self.varvalue), self.unoutputvars
                                                ] if type(self.unoutputvars) != type(None) else None
        else:
            self.yref, self.xref, self.bref = self.data.loc[:, self.outputvars
                                        ], self.data.loc[:, self.inputvars
                                        ], self.data.loc[:, self.unoutputvars
                                        ] if type(self.unoutputvars) != type(None) else None

        self.xcol = self.x.columns
        self.ycol = self.y.columns
        self.bcol = self.b.columns if type(self.unoutputvars) != type(None) else None

        # print(self.xcol)

        self.I = self.x.index          ## I 是 被评价决策单元的索引
        self.__modeldict = {}
        for i in self.I:
            # print(i)
            self.I0 = i                                                 ## I 是 被评价决策单元的数量

            self.__model__ = ConcreteModel()
            # Initialize sets
            self.__model__.I2 = Set(initialize= self.xref.index)                     ## I2 是 参考决策单元的数量
            self.__model__.K = Set(initialize=range(len(self.x.iloc[0])))          ## K 是投入个数
            self.__model__.L = Set(initialize=range(len(self.y.iloc[0])))          ## L 是产出个数 被评价单元和参考单元的K，L一样

            if type(self.b) != type(None):
                self.__model__.B = Set(initialize=range(len(self.b.iloc[0])))      ## B 是 非期望产出个数

            # Initialize variable

            self.__model__.objx = Var(self.__model__.K,bounds=(0.0, None), within=Reals,doc='object x')

            self.__model__.thetax = Var(self.__model__.K,bounds=(0.0, None), doc='slack x')
            self.__model__.thetay = Var(self.__model__.L,bounds=(0.0, None), doc='slack y')
            if type(self.b) != type(None):
                self.__model__.thetab = Var(self.__model__.B,bounds=(0.0, None), doc='slack b')
                self.__model__.theta = Var(self.__model__.B,bounds=(0.0, None), within=Reals,doc='object b')

            self.__model__.lamda = Var(self.__model__.I2, bounds=(0.0, None),within=Reals, doc='intensity variables')


            # Setup the objective function and constraints
            self.__model__.objective = Objective(rule=self.__objective_rule(), sense=minimize, doc='objective function')
            self.__model__.input = Constraint(self.__model__.K,  rule=self.__input_rule(), doc='input constraint')
            self.__model__.output = Constraint(self.__model__.L,  rule=self.__output_rule(), doc='output constraint')


            if type(self.b) != type(None):
                self.__model__.undesirable_output = Constraint(self.__model__.B, rule=self.__undesirable_output_rule(), \
                                                               doc='undesirable output constraint')
                self.__model__.mb = Constraint(self.__model__.B,  rule=self.__mb_rule(), \
                                                                doc='material balance constraint')
            if self.rts == RTS_VRS1:
                self.__model__.vrs = Constraint(rule=self.__vrs_rule(), doc='various return to scale rule')

            self.__modeldict[i] = self.__model__

        # Optimize model
    def __objective_rule(self):
        """Return the proper objective function"""
        def objective_rule(model):
            return sum(model.theta[b]*1 for b in model.B)
        return objective_rule

    def __input_rule(self):
        """Return the proper input constraint"""
        def input_rule(model, k):
            return sum(model.lamda[i2] * self.xref.loc[i2,self.xcol[k]] for i2 in model.I2
                    ) + model.thetax[k] == model.objx[k]
        return input_rule

    def __output_rule(self):
        """Return the proper output constraint"""
        def output_rule(model, l):
            return sum(model.lamda[i2] * self.yref.loc[i2,self.ycol[l]] for i2 in model.I2
                    ) - model.thetay[l] == self.y.loc[self.I0,self.ycol[l]]
        return output_rule

    def __undesirable_output_rule(self):
        """Return the proper undesirable output constraint"""
        def undesirable_output_rule(model, b):
            return sum(model.lamda[i2] * self.bref.loc[i2,self.bcol[b]] for i2 in model.I2
                    ) + model.thetab[b] == model.theta[b]*1
        return undesirable_output_rule

    def __mb_rule(self):
        """Return the proper undesirable output constraint"""
        def mb_rule(model, b):
            return sum(self.sx[b][k] * (model.thetax[k] + self.x.loc[self.I0,self.xcol[k]] - model.objx[k] ) for k in model.K) \
                    + sum(self.sy[b][l] * model.thetay[l] for l in model.L) \
                    == self.b.loc[self.I0,self.bcol[b]] - model.theta[b]
        return mb_rule

    def __vrs_rule(self):
        def vrs_rule(model):
            return sum(model.lamda[ i2] for i2 in model.I2) == 1

        return vrs_rule

    def optimize(self,  solver=OPT_DEFAULT):
        """Optimize the function by requested method

        Args:
            solver (string): The solver chosen for optimization. It will optimize with default solver if OPT_DEFAULT is given.
        """
        # TODO(error/warning handling): Check problem status after optimization

        data2,obj,theta,thetax,thetay,thetab,lamda = pd.DataFrame(),{},{},{},{},{},{}
        for ind, problem in self.__modeldict.items():
            _, data2.loc[ind,"optimization_status"] = tools.optimize_model(problem, ind, solver)

            if type(self.b) != type(None):
                obj[ind] = problem.objective()
                theta[ind], = np.asarray(list(problem.theta[:].value))
                thetax[ind]= np.asarray(list(problem.thetax[:].value))
                thetay[ind]= np.asarray(list(problem.thetay[:].value))
                thetab[ind]= np.asarray(list(problem.thetab[:].value))
                lamda[ind]= np.asarray(list(problem.lamda[:].value))
            else:
                obj[ind] = problem.objective()
                theta[ind], = np.asarray(list(problem.theta[:].value))
                thetax[ind]= np.asarray(list(problem.thetax[:].value))
                thetay[ind]= np.asarray(list(problem.thetay[:].value))
                lamda[ind]= np.asarray(list(problem.lamda[:].value))

                # print(list(problem.thetax[:].value ),list(problem.t[:].value ))
        obj = pd.DataFrame(obj,index=["obj"]).T
        theta = pd.DataFrame(theta,index=["var Undesirable"]).T
        thetax = pd.DataFrame(thetax).T
        thetax.columns = thetax.columns.map(lambda x : "Input"+ str(x)+"'s slack" )
        thetay = pd.DataFrame(thetay).T
        thetay.columns = thetay.columns.map(lambda y : "Output"+ str(y)+"'s slack" )

        theta_=pd.concat([theta,thetax],axis=1)
        theta_=pd.concat([theta_,thetay],axis=1)
        if type(self.b) != type(None):
            thetab = pd.DataFrame(thetab).T
            thetab.columns = thetab.columns.map(lambda b : "Undesirable Output"+ str(b)+"'s slack" )
            theta_ = pd.concat([theta_,thetab],axis=1)

        lamda =pd.DataFrame(lamda) .T
        lamda.columns = lamda.columns.map(lambda x : "lamda"+ str(x) )
        data3 = pd.concat([data2,obj],axis=1)
        data3 = pd.concat([data3,theta_],axis=1)
        # data3 = pd.concat([data3,lamda],axis=1)
        return data3


    def info(self, dmu = "all"):
        """Show the infomation of the lp model

        Args:
            dmu (string): The solver chosen for optimization. Default is "all".
        """
        if dmu =="all":
            for ind, problem in self.__modeldict.items():
                # print(ind,"\n",problem.pprint())
                pass
        print(self.__modeldict[int(dmu)].pprint())


class MBx2(MB):
    def __init__(self, data,year,sent = "inputvar=outputvar",  sx=[[1,1,1],[1,1,1]], sy=[[1],[1]], rts=RTS_VRS1, baseindex=None,refindex=None):
        """DEA: Directional distance function

        Args:
            data (pandas.DataFrame): input pandas.
            sent (str): inputvars=outputvars[: unoutputvars]. e.g.: "K L = Y : CO2"
            sx (list): 投入包含污染物质系数. Defaults to [[1,1,1],[1,1,1]].
            sy (list, optional): 期望产出包含污染物质系数. Defaults to [[1],[1]].
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # Initialize DEA model
        self.data=data
        self.year = year
        self.sent = sent
        self.tlt=pd.Series(self.year).drop_duplicates().sort_values()
        self.inputvars = self.sent.split('=')[0].strip(' ').split(' ')
        try:
            self.outputvars = self.sent.split('=')[1]   .split(':')[0].strip(' ').split(' ')
            self.unoutputvars = self.sent.split('=')[1]   .split(':')[1].strip(' ').split(' ')
        except:
            self.outputvars = self.sent.split('=')[1]    .strip(' ').split(' ')
            self.unoutputvars=None
        self.sx, self.sy = sx, sy
        self.rts = rts
        # print(self.sx, self.sy)

        self.baseindex = baseindex
        if type(baseindex) != type(None):
            self.varname1=self.baseindex.split('=')[0]
            # print(self.baseindex)
            self.varvalue1=ast.literal_eval(self.baseindex.split('=')[1])
            self.y, self.x, self.b = self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.outputvars
                                        ], self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.inputvars
                                        ], self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.unoutputvars
                                        ]if type(self.unoutputvars) != type(None) else None

        else:

            self.y, self.x, self.b = self.data.loc[:, self.outputvars
                                        ], self.data.loc[:, self.inputvars
                                        ], self.data.loc[:, self.unoutputvars
                                        ] if type(self.unoutputvars) != type(None) else None


        # print(self.b)
        self.refindex = refindex
        if type(refindex) != type(None):
            self.varname=self.refindex.split('=')[0]
            self.varvalue=ast.literal_eval(self.refindex.split('=')[1])

            self.yref, self.xref, self.bref = self.data.loc[self.data[self.varname].isin(self.varvalue), self.outputvars
                                                ], self.data.loc[self.data[self.varname].isin(self.varvalue), self.inputvars
                                                ], self.data.loc[self.data[self.varname].isin(self.varvalue), self.unoutputvars
                                                ] if type(self.unoutputvars) != type(None) else None
        else:
            self.yref, self.xref, self.bref = self.data.loc[:, self.outputvars
                                        ], self.data.loc[:, self.inputvars
                                        ], self.data.loc[:, self.unoutputvars
                                        ] if type(self.unoutputvars) != type(None) else None

        self.xcol = self.x.columns
        self.ycol = self.y.columns
        self.bcol = self.b.columns if type(self.unoutputvars) != type(None) else None

        # print(self.xcol)

        self.I = self.x.index          ## I 是 被评价决策单元的索引
        self.__modeldict = {}
        for i in self.I:
            # print(i)
            self.I0 = i                                                 ## I 是 被评价决策单元的数量

            self.__model__ = ConcreteModel()
            # Initialize sets
            self.__model__.I2 = Set(initialize= self.xref.index)                     ## I2 是 参考决策单元的数量
            self.__model__.K = Set(initialize=range(len(self.x.iloc[0])))          ## K 是投入个数
            self.__model__.L = Set(initialize=range(len(self.y.iloc[0])))          ## L 是产出个数 被评价单元和参考单元的K，L一样

            if type(self.b) != type(None):
                self.__model__.B = Set(initialize=range(len(self.b.iloc[0])))      ## B 是 非期望产出个数

            # Initialize variable

            self.__model__.objx = Var(self.__model__.K,bounds=(0.0, None), within=Reals,doc='object x')

            self.__model__.thetax = Var(self.__model__.K,bounds=(0.0, None), doc='slack x')
            self.__model__.thetay = Var(self.__model__.L,bounds=(0.0, None), doc='slack y')
            if type(self.b) != type(None):
                self.__model__.thetab = Var(self.__model__.B,bounds=(0.0, None), doc='slack b')
                self.__model__.theta = Var(self.__model__.B,bounds=(0.0, None), within=Reals,doc='object b')

            self.__model__.lamda = Var(self.__model__.I2, bounds=(0.0, None),within=Reals, doc='intensity variables')


            # Setup the objective function and constraints
            self.__model__.objective = Objective(rule=self.__objective_rule(), sense=minimize, doc='objective function')
            self.__model__.input = Constraint(self.__model__.K,  rule=self.__input_rule(), doc='input constraint')
            self.__model__.output = Constraint(self.__model__.L,  rule=self.__output_rule(), doc='output constraint')


            if type(self.b) != type(None):
                self.__model__.undesirable_output = Constraint(self.__model__.B, rule=self.__undesirable_output_rule(), \
                                                               doc='undesirable output constraint')
                self.__model__.mb = Constraint(self.__model__.B,  rule=self.__mb_rule(), \
                                                                doc='material balance constraint')
            if self.rts == RTS_VRS1:
                self.__model__.vrs = Constraint(rule=self.__vrs_rule(), doc='various return to scale rule')

            self.__modeldict[i] = self.__model__

        # Optimize model
    def __objective_rule(self):
        """Return the proper objective function"""
        def objective_rule(model):
            return sum(model.theta[b]*1 for b in model.B)
        return objective_rule

    def __input_rule(self):
        """Return the proper input constraint"""
        def input_rule(model, k):
            return sum(model.lamda[i2] * self.xref.loc[i2,self.xcol[k]] for i2 in model.I2
                    ) + model.thetax[k] == model.objx[k]
        return input_rule

    def __output_rule(self):
        """Return the proper output constraint"""
        def output_rule(model, l):
            return sum(model.lamda[i2] * self.yref.loc[i2,self.ycol[l]] for i2 in model.I2
                    ) - model.thetay[l] == self.y.loc[self.I0,self.ycol[l]]
        return output_rule

    def __undesirable_output_rule(self):
        """Return the proper undesirable output constraint"""
        def undesirable_output_rule(model, b):
            return sum(model.lamda[i2] * self.bref.loc[i2,self.bcol[b]] for i2 in model.I2
                    ) + model.thetab[b] == model.theta[b]*1
        return undesirable_output_rule

    def __mb_rule(self):
        """Return the proper undesirable output constraint"""
        def mb_rule(model, b):
            return sum(self.sx[b][k] * (model.thetax[k] + self.x.loc[self.I0,self.xcol[k]] - model.objx[k] ) for k in model.K) \
                    + sum(self.sy[b][l] * model.thetay[l] for l in model.L) \
                    == self.b.loc[self.I0,self.bcol[b]] - model.theta[b]
        return mb_rule

    def __vrs_rule(self):
        def vrs_rule(model):
            return sum(model.lamda[ i2] for i2 in model.I2) == 1

        return vrs_rule

    def optimize(self,  solver=OPT_DEFAULT):
        """Optimize the function by requested method

        Args:
            solver (string): The solver chosen for optimization. It will optimize with default solver if OPT_DEFAULT is given.
        """
        # TODO(error/warning handling): Check problem status after optimization

        data2,obj,theta,thetax,thetay,thetab,lamda = pd.DataFrame(),{},{},{},{},{},{}
        for ind, problem in self.__modeldict.items():
            _, data2.loc[ind,"optimization_status"] = tools.optimize_model(problem, ind, solver)

            if type(self.b) != type(None):
                obj[ind] = problem.objective()
                theta[ind], = np.asarray(list(problem.theta[:].value))
                thetax[ind]= np.asarray(list(problem.thetax[:].value))
                thetay[ind]= np.asarray(list(problem.thetay[:].value))
                thetab[ind]= np.asarray(list(problem.thetab[:].value))
                lamda[ind]= np.asarray(list(problem.lamda[:].value))
            else:
                obj[ind] = problem.objective()
                theta[ind], = np.asarray(list(problem.theta[:].value))
                thetax[ind]= np.asarray(list(problem.thetax[:].value))
                thetay[ind]= np.asarray(list(problem.thetay[:].value))
                lamda[ind]= np.asarray(list(problem.lamda[:].value))

                # print(list(problem.thetax[:].value ),list(problem.t[:].value ))
        obj = pd.DataFrame(obj,index=["obj"]).T
        theta = pd.DataFrame(theta,index=["var Undesirable"]).T
        thetax = pd.DataFrame(thetax).T
        thetax.columns = thetax.columns.map(lambda x : "Input"+ str(x)+"'s slack" )
        thetay = pd.DataFrame(thetay).T
        thetay.columns = thetay.columns.map(lambda y : "Output"+ str(y)+"'s slack" )

        theta_=pd.concat([theta,thetax],axis=1)
        theta_=pd.concat([theta_,thetay],axis=1)
        if type(self.b) != type(None):
            thetab = pd.DataFrame(thetab).T
            thetab.columns = thetab.columns.map(lambda b : "Undesirable Output"+ str(b)+"'s slack" )
            theta_ = pd.concat([theta_,thetab],axis=1)

        lamda =pd.DataFrame(lamda) .T
        lamda.columns = lamda.columns.map(lambda x : "lamda"+ str(x) )
        data3 = pd.concat([data2,obj],axis=1)
        data3 = pd.concat([data3,theta_],axis=1)
        # data3 = pd.concat([data3,lamda],axis=1)
        return data3


    def info(self, dmu = "all"):
        """Show the infomation of the lp model

        Args:
            dmu (string): The solver chosen for optimization. Default is "all".
        """
        if dmu =="all":
            for ind, problem in self.__modeldict.items():
                print(ind,"\n",problem.pprint())

        print(self.__modeldict[int(dmu)].pprint())


class MBxy(MB):
    def __init__(self, data,year,sent = "inputvar=outputvar",  sx=[[1,1,1],[1,1,1]], sy=[[1],[1]], rts=RTS_VRS1, baseindex=None,refindex=None):
        """DEA: Directional distance function

        Args:
            data (pandas.DataFrame): input pandas.
            sent (str): inputvars=outputvars[: unoutputvars]. e.g.: "K L = Y : CO2"
            sx (list): 投入包含污染物质系数. Defaults to [[1,1,1],[1,1,1]].
            sy (list, optional): 期望产出包含污染物质系数. Defaults to [[1],[1]].
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # Initialize DEA model
        self.data=data
        self.year = year
        self.sent = sent
        self.tlt=pd.Series(self.year).drop_duplicates().sort_values()
        self.inputvars = self.sent.split('=')[0].strip(' ').split(' ')
        try:
            self.outputvars = self.sent.split('=')[1]   .split(':')[0].strip(' ').split(' ')
            self.unoutputvars = self.sent.split('=')[1]   .split(':')[1].strip(' ').split(' ')
        except:
            self.outputvars = self.sent.split('=')[1]    .strip(' ').split(' ')
            self.unoutputvars=None
        self.sx, self.sy = sx, sy
        self.rts = rts
        # print(self.sx, self.sy)

        self.baseindex = baseindex
        if type(baseindex) != type(None):
            self.varname1=self.baseindex.split('=')[0]
            # print(self.baseindex)
            self.varvalue1=ast.literal_eval(self.baseindex.split('=')[1])
            self.y, self.x, self.b = self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.outputvars
                                        ], self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.inputvars
                                        ], self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.unoutputvars
                                        ]if type(self.unoutputvars) != type(None) else None

        else:

            self.y, self.x, self.b = self.data.loc[:, self.outputvars
                                        ], self.data.loc[:, self.inputvars
                                        ], self.data.loc[:, self.unoutputvars
                                        ] if type(self.unoutputvars) != type(None) else None


        # print(self.b)
        self.refindex = refindex
        if type(refindex) != type(None):
            self.varname=self.refindex.split('=')[0]
            self.varvalue=ast.literal_eval(self.refindex.split('=')[1])

            self.yref, self.xref, self.bref = self.data.loc[self.data[self.varname].isin(self.varvalue), self.outputvars
                                                ], self.data.loc[self.data[self.varname].isin(self.varvalue), self.inputvars
                                                ], self.data.loc[self.data[self.varname].isin(self.varvalue), self.unoutputvars
                                                ] if type(self.unoutputvars) != type(None) else None
        else:
            self.yref, self.xref, self.bref = self.data.loc[:, self.outputvars
                                        ], self.data.loc[:, self.inputvars
                                        ], self.data.loc[:, self.unoutputvars
                                        ] if type(self.unoutputvars) != type(None) else None

        self.xcol = self.x.columns
        self.ycol = self.y.columns
        self.bcol = self.b.columns if type(self.unoutputvars) != type(None) else None

        # print(self.xcol)

        self.I = self.x.index          ## I 是 被评价决策单元的索引
        self.__modeldict = {}
        for i in self.I:
            # print(i)
            self.I0 = i                                                 ## I 是 被评价决策单元的数量

            self.__model__ = ConcreteModel()
            # Initialize sets
            self.__model__.I2 = Set(initialize= self.xref.index)                     ## I2 是 参考决策单元的数量
            self.__model__.K = Set(initialize=range(len(self.x.iloc[0])))          ## K 是投入个数
            self.__model__.L = Set(initialize=range(len(self.y.iloc[0])))          ## L 是产出个数 被评价单元和参考单元的K，L一样

            if type(self.b) != type(None):
                self.__model__.B = Set(initialize=range(len(self.b.iloc[0])))      ## B 是 非期望产出个数

            # Initialize variable

            self.__model__.objx = Var(self.__model__.K,bounds=(0.0, None), within=Reals,doc='object x')
            self.__model__.objy = Var(self.__model__.L,bounds=(0.0, None), within=Reals,doc='object y')

            self.__model__.thetax = Var(self.__model__.K,bounds=(0.0, None), doc='slack x')
            self.__model__.thetay = Var(self.__model__.L,bounds=(0.0, None), doc='slack y')
            if type(self.b) != type(None):
                self.__model__.thetab = Var(self.__model__.B,bounds=(0.0, None), doc='slack b')
                self.__model__.theta = Var(self.__model__.B,bounds=(0.0, None), within=Reals,doc='object b')

            self.__model__.lamda = Var(self.__model__.I2, bounds=(0.0, None),within=Reals, doc='intensity variables')


            # Setup the objective function and constraints
            self.__model__.objective = Objective(rule=self.__objective_rule(), sense=minimize, doc='objective function')
            self.__model__.input = Constraint(self.__model__.K,  rule=self.__input_rule(), doc='input constraint')
            self.__model__.output = Constraint(self.__model__.L,  rule=self.__output_rule(), doc='output constraint')


            if type(self.b) != type(None):
                self.__model__.undesirable_output = Constraint(self.__model__.B, rule=self.__undesirable_output_rule(), \
                                                               doc='undesirable output constraint')
                self.__model__.mb = Constraint(self.__model__.B,  rule=self.__mb_rule(), \
                                                                doc='material balance constraint')
            if self.rts == RTS_VRS1:
                self.__model__.vrs = Constraint(rule=self.__vrs_rule(), doc='various return to scale rule')

            self.__modeldict[i] = self.__model__

        # Optimize model
    def __objective_rule(self):
        """Return the proper objective function"""
        def objective_rule(model):
            return sum(model.theta[b]*1 for b in model.B)
        return objective_rule

    def __input_rule(self):
        """Return the proper input constraint"""
        def input_rule(model, k):
            return sum(model.lamda[i2] * self.xref.loc[i2,self.xcol[k]] for i2 in model.I2
                    ) + model.thetax[k] == model.objx[k]
        return input_rule

    def __output_rule(self):
        """Return the proper output constraint"""
        def output_rule(model, l):
            return sum(model.lamda[i2] * self.yref.loc[i2,self.ycol[l]] for i2 in model.I2
                    ) - model.thetay[l] == model.objy[l]
        return output_rule

    def __undesirable_output_rule(self):
        """Return the proper undesirable output constraint"""
        def undesirable_output_rule(model, b):
            return sum(model.lamda[i2] * self.bref.loc[i2,self.bcol[b]] for i2 in model.I2
                    ) + model.thetab[b] == model.theta[b]*1
        return undesirable_output_rule

    def __mb_rule(self):
        """Return the proper undesirable output constraint"""
        def mb_rule(model, b):
            return sum(self.sx[b][k] * (model.thetax[k] + self.x.loc[self.I0,self.xcol[k]] - model.objx[k]) for k in model.K) \
                 + sum(self.sy[b][l] * (model.thetay[l] + model.objy[l] - self.y.loc[self.I0,self.ycol[l]]) for l in model.L) \
                == model.thetab[b]+self.b.loc[self.I0,self.bcol[b]] - model.theta[b]
        return mb_rule


    def __vrs_rule(self):
        def vrs_rule(model):
            return sum(model.lamda[i2] for i2 in model.I2) == 1

        return vrs_rule

    def optimize(self,  solver=OPT_DEFAULT):
        """Optimize the function by requested method

        Args:
            solver (string): The solver chosen for optimization. It will optimize with default solver if OPT_DEFAULT is given.
        """
        # TODO(error/warning handling): Check problem status after optimization

        data2,obj,theta,thetax,thetay,thetab,lamda = pd.DataFrame(),{},{},{},{},{},{}
        for ind, problem in self.__modeldict.items():
            _, data2.loc[ind,"optimization_status"] = tools.optimize_model(problem, ind, solver)

            if type(self.b) != type(None):
                obj[ind] = problem.objective()
                theta[ind], = np.asarray(list(problem.theta[:].value))
                thetax[ind]= np.asarray(list(problem.thetax[:].value))
                thetay[ind]= np.asarray(list(problem.thetay[:].value))
                thetab[ind]= np.asarray(list(problem.thetab[:].value))
                lamda[ind]= np.asarray(list(problem.lamda[:].value))
            else:
                obj[ind] = problem.objective()
                theta[ind], = np.asarray(list(problem.theta[:].value))
                thetax[ind]= np.asarray(list(problem.thetax[:].value))
                thetay[ind]= np.asarray(list(problem.thetay[:].value))
                lamda[ind]= np.asarray(list(problem.lamda[:].value))

                # print(list(problem.thetax[:].value ),list(problem.t[:].value ))
        obj = pd.DataFrame(obj,index=["obj"]).T
        theta = pd.DataFrame(theta,index=["var Undesirable"]).T
        thetax = pd.DataFrame(thetax).T
        thetax.columns = thetax.columns.map(lambda x : "Input"+ str(x)+"'s slack" )
        thetay = pd.DataFrame(thetay).T
        thetay.columns = thetay.columns.map(lambda y : "Output"+ str(y)+"'s slack" )

        theta_=pd.concat([theta,thetax],axis=1)
        theta_=pd.concat([theta_,thetay],axis=1)
        if type(self.b) != type(None):
            thetab = pd.DataFrame(thetab).T
            thetab.columns = thetab.columns.map(lambda b : "Undesirable Output"+ str(b)+"'s slack" )
            theta_ = pd.concat([theta_,thetab],axis=1)

        lamda =pd.DataFrame(lamda) .T
        lamda.columns = lamda.columns.map(lambda x : "lamda"+ str(x) )
        data3 = pd.concat([data2,obj],axis=1)
        data3 = pd.concat([data3,theta_],axis=1)
        # data3 = pd.concat([data3,lamda],axis=1)
        return data3


    def info(self, dmu = "all"):
        """Show the infomation of the lp model

        Args:
            dmu (string): The solver chosen for optimization. Default is "all".
        """
        if dmu =="all":
            for ind, problem in self.__modeldict.items():
                print(ind,"\n",problem.pprint())

        print(self.__modeldict[int(dmu)].pprint())



class MBxy2(MB):
    def __init__(self, data,year,sent = "inputvar=outputvar",  sx=[[1,1,1],[1,1,1]], sy=[[1],[1]], rts=RTS_VRS1, baseindex=None,refindex=None):
        """DEA: Directional distance function

        Args:
            data (pandas.DataFrame): input pandas.
            sent (str): inputvars=outputvars[: unoutputvars]. e.g.: "K L = Y : CO2"
            sx (list): 投入包含污染物质系数. Defaults to [[1,1,1],[1,1,1]].
            sy (list, optional): 期望产出包含污染物质系数. Defaults to [[1],[1]].
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # Initialize DEA model
        self.data=data
        self.year = year
        self.sent = sent
        self.tlt=pd.Series(self.year).drop_duplicates().sort_values()
        self.inputvars = self.sent.split('=')[0].strip(' ').split(' ')
        try:
            self.outputvars = self.sent.split('=')[1]   .split(':')[0].strip(' ').split(' ')
            self.unoutputvars = self.sent.split('=')[1]   .split(':')[1].strip(' ').split(' ')
        except:
            self.outputvars = self.sent.split('=')[1]    .strip(' ').split(' ')
            self.unoutputvars=None
        self.sx, self.sy = sx, sy
        self.rts = rts
        # print(self.sx, self.sy)

        self.baseindex = baseindex
        if type(baseindex) != type(None):
            self.varname1=self.baseindex.split('=')[0]
            # print(self.baseindex)
            self.varvalue1=ast.literal_eval(self.baseindex.split('=')[1])
            self.y, self.x, self.b = self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.outputvars
                                        ], self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.inputvars
                                        ], self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.unoutputvars
                                        ]if type(self.unoutputvars) != type(None) else None

        else:

            self.y, self.x, self.b = self.data.loc[:, self.outputvars
                                        ], self.data.loc[:, self.inputvars
                                        ], self.data.loc[:, self.unoutputvars
                                        ] if type(self.unoutputvars) != type(None) else None


        # print(self.b)
        self.refindex = refindex
        if type(refindex) != type(None):
            self.varname=self.refindex.split('=')[0]
            self.varvalue=ast.literal_eval(self.refindex.split('=')[1])

            self.yref, self.xref, self.bref = self.data.loc[self.data[self.varname].isin(self.varvalue), self.outputvars
                                                ], self.data.loc[self.data[self.varname].isin(self.varvalue), self.inputvars
                                                ], self.data.loc[self.data[self.varname].isin(self.varvalue), self.unoutputvars
                                                ] if type(self.unoutputvars) != type(None) else None
        else:
            self.yref, self.xref, self.bref = self.data.loc[:, self.outputvars
                                        ], self.data.loc[:, self.inputvars
                                        ], self.data.loc[:, self.unoutputvars
                                        ] if type(self.unoutputvars) != type(None) else None

        self.xcol = self.x.columns
        self.ycol = self.y.columns
        self.bcol = self.b.columns if type(self.unoutputvars) != type(None) else None

        # print(self.xcol)

        self.I = self.x.index          ## I 是 被评价决策单元的索引
        self.__modeldict = {}
        for i in self.I:
            # print(i)
            self.I0 = i                                                 ## I 是 被评价决策单元的数量

            self.__model__ = ConcreteModel()
            # Initialize sets
            self.__model__.I2 = Set(initialize= self.xref.index)                     ## I2 是 参考决策单元的数量
            self.__model__.K = Set(initialize=range(len(self.x.iloc[0])))          ## K 是投入个数
            self.__model__.L = Set(initialize=range(len(self.y.iloc[0])))          ## L 是产出个数 被评价单元和参考单元的K，L一样

            if type(self.b) != type(None):
                self.__model__.B = Set(initialize=range(len(self.b.iloc[0])))      ## B 是 非期望产出个数

            # Initialize variable

            self.__model__.objx = Var(self.__model__.K,bounds=(0.0, None), within=Reals,doc='object x')
            self.__model__.objy = Var(self.__model__.L,bounds=(0.0, None), within=Reals,doc='object y')

            self.__model__.thetax = Var(self.__model__.K,bounds=(0.0, None), doc='slack x')
            self.__model__.thetay = Var(self.__model__.L,bounds=(0.0, None), doc='slack y')
            if type(self.b) != type(None):
                self.__model__.thetab = Var(self.__model__.B,bounds=(0.0, None), doc='slack b')

            self.__model__.lamda = Var(self.__model__.I2, bounds=(0.0, None),within=Reals, doc='intensity variables')


            # Setup the objective function and constraints
            self.__model__.objective = Objective(rule=self.__objective_rule(), sense=maximize, doc='objective function')
            self.__model__.input = Constraint(self.__model__.K,  rule=self.__input_rule(), doc='input constraint')
            self.__model__.output = Constraint(self.__model__.L,  rule=self.__output_rule(), doc='output constraint')


            if type(self.b) != type(None):
                self.__model__.undesirable_output = Constraint(self.__model__.B, rule=self.__undesirable_output_rule(), \
                                                               doc='undesirable output constraint')
                self.__model__.mb = Constraint(self.__model__.B,  rule=self.__mb_rule(), \
                                                                doc='material balance constraint')
            if self.rts == RTS_VRS1:
                self.__model__.vrs = Constraint(rule=self.__vrs_rule(), doc='various return to scale rule')

            self.__modeldict[i] = self.__model__

        # Optimize model
    def __objective_rule(self):
        """Return the proper objective function"""
        def objective_rule(model):
            return sum(model.thetab[b] for b in model.B)
        return objective_rule

    def __input_rule(self):
        """Return the proper input constraint"""
        def input_rule(model, k):
            return sum(model.lamda[i2] * self.xref.loc[i2,self.xcol[k]] for i2 in model.I2
                    ) + model.thetax[k] == model.objx[k]
        return input_rule

    def __output_rule(self):
        """Return the proper output constraint"""
        def output_rule(model, l):
            return sum(model.lamda[i2] * self.yref.loc[i2,self.ycol[l]] for i2 in model.I2
                    ) - model.thetay[l] == model.objy[l]
        return output_rule

    def __undesirable_output_rule(self):
        """Return the proper undesirable output constraint"""
        def undesirable_output_rule(model, b):
            return sum(model.lamda[i2] * self.bref.loc[i2,self.bcol[b]] for i2 in model.I2
                    ) + model.thetab[b] == self.b.loc[self.I0,self.bcol[b]]
        return undesirable_output_rule

    def __mb_rule(self):
        """Return the proper undesirable output constraint"""
        def mb_rule(model, b):
            return sum(self.sx[b][k] * (model.thetax[k] + self.x.loc[self.I0,self.xcol[k]] - model.objx[k]) for k in model.K) \
                 + sum(self.sy[b][l] * (model.thetay[l] + model.objy[l] - self.y.loc[self.I0,self.ycol[l]]) for l in model.L) \
                == model.thetab[b]
        return mb_rule


    def __vrs_rule(self):
        def vrs_rule(model):
            return sum(model.lamda[i2] for i2 in model.I2) == 1

        return vrs_rule

    def optimize(self,  solver=OPT_DEFAULT):
        """Optimize the function by requested method

        Args:
            solver (string): The solver chosen for optimization. It will optimize with default solver if OPT_DEFAULT is given.
        """
        # TODO(error/warning handling): Check problem status after optimization

        data2,obj,thetax,thetay,thetab,lamda = pd.DataFrame(),{},{},{},{},{}
        for ind, problem in self.__modeldict.items():
            _, data2.loc[ind,"optimization_status"] = tools.optimize_model(problem, ind, solver)

            if type(self.b) != type(None):
                obj[ind] = problem.objective()
                # theta[ind], = np.asarray(list(problem.theta[:].value))
                thetax[ind]= np.asarray(list(problem.thetax[:].value))
                thetay[ind]= np.asarray(list(problem.thetay[:].value))
                thetab[ind]= np.asarray(list(problem.thetab[:].value))
                lamda[ind]= np.asarray(list(problem.lamda[:].value))
            else:
                obj[ind] = problem.objective()
                # theta[ind], = np.asarray(list(problem.theta[:].value))
                thetax[ind]= np.asarray(list(problem.thetax[:].value))
                thetay[ind]= np.asarray(list(problem.thetay[:].value))
                lamda[ind]= np.asarray(list(problem.lamda[:].value))

                # print(list(problem.thetax[:].value ),list(problem.t[:].value ))
        obj = pd.DataFrame(obj,index=["obj"]).T
        # theta = pd.DataFrame(theta,index=["var Undesirable"]).T
        thetax = pd.DataFrame(thetax).T
        thetax.columns = thetax.columns.map(lambda x : "Input"+ str(x)+"'s slack" )
        thetay = pd.DataFrame(thetay).T
        thetay.columns = thetay.columns.map(lambda y : "Output"+ str(y)+"'s slack" )

        # theta_=pd.concat([theta,thetax],axis=1)
        theta_=pd.concat([thetax,thetay],axis=1)
        if type(self.b) != type(None):
            thetab = pd.DataFrame(thetab).T
            thetab.columns = thetab.columns.map(lambda b : "Undesirable Output"+ str(b)+"'s slack" )
            theta_ = pd.concat([theta_,thetab],axis=1)

        lamda =pd.DataFrame(lamda) .T
        lamda.columns = lamda.columns.map(lambda x : "lamda"+ str(x) )
        data3 = pd.concat([data2,obj],axis=1)
        data3 = pd.concat([data3,theta_],axis=1)
        # data3 = pd.concat([data3,lamda],axis=1)
        return data3


    def info(self, dmu = "all"):
        """Show the infomation of the lp model

        Args:
            dmu (string): The solver chosen for optimization. Default is "all".
        """
        if dmu =="all":
            for ind, problem in self.__modeldict.items():
                print(ind,"\n",problem.pprint())

        print(self.__modeldict[int(dmu)].pprint())

class MB2(MB):
    def __init__(self, data,year,sent = "inputvar=outputvar",  sx=[[1,1,1],[1,1,1]], sy=[[1],[1]], rts=RTS_VRS1, baseindex=None,refindex=None):
        """DEA: Directional distance function

        Args:
            data (pandas.DataFrame): input pandas.
            sent (str): inputvars=outputvars[: unoutputvars]. e.g.: "K L = Y : CO2"
            sx (list): 投入包含污染物质系数. Defaults to [[1,1,1],[1,1,1]].
            sy (list, optional): 期望产出包含污染物质系数. Defaults to [[1],[1]].
            rts (String): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale)
            baseindex (String, optional): estimate index. Defaults to None. e.g.: "Year=[2010]"
            refindex (String, optional): reference index. Defaults to None. e.g.: "Year=[2010]"
        """
        # Initialize DEA model
        self.data=data
        self.year = year
        self.sent = sent
        self.tlt=pd.Series(self.year).drop_duplicates().sort_values()
        self.inputvars = self.sent.split('=')[0].strip(' ').split(' ')
        try:
            self.outputvars = self.sent.split('=')[1]   .split(':')[0].strip(' ').split(' ')
            self.unoutputvars = self.sent.split('=')[1]   .split(':')[1].strip(' ').split(' ')
        except:
            self.outputvars = self.sent.split('=')[1]    .strip(' ').split(' ')
            self.unoutputvars=None
        self.sx, self.sy = sx, sy
        self.rts = rts
        # print(self.sx, self.sy)

        self.baseindex = baseindex
        if type(baseindex) != type(None):
            self.varname1=self.baseindex.split('=')[0]
            # print(self.baseindex)
            self.varvalue1=ast.literal_eval(self.baseindex.split('=')[1])
            self.y, self.x, self.b = self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.outputvars
                                        ], self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.inputvars
                                        ], self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.unoutputvars
                                        ]if type(self.unoutputvars) != type(None) else None

        else:

            self.y, self.x, self.b = self.data.loc[:, self.outputvars
                                        ], self.data.loc[:, self.inputvars
                                        ], self.data.loc[:, self.unoutputvars
                                        ] if type(self.unoutputvars) != type(None) else None


        # print(self.b)
        self.refindex = refindex
        if type(refindex) != type(None):
            self.varname=self.refindex.split('=')[0]
            self.varvalue=ast.literal_eval(self.refindex.split('=')[1])

            self.yref, self.xref, self.bref = self.data.loc[self.data[self.varname].isin(self.varvalue), self.outputvars
                                                ], self.data.loc[self.data[self.varname].isin(self.varvalue), self.inputvars
                                                ], self.data.loc[self.data[self.varname].isin(self.varvalue), self.unoutputvars
                                                ] if type(self.unoutputvars) != type(None) else None
        else:
            self.yref, self.xref, self.bref = self.data.loc[:, self.outputvars
                                        ], self.data.loc[:, self.inputvars
                                        ], self.data.loc[:, self.unoutputvars
                                        ] if type(self.unoutputvars) != type(None) else None

        self.xcol = self.x.columns
        self.ycol = self.y.columns
        self.bcol = self.b.columns if type(self.unoutputvars) != type(None) else None

        # print(self.xcol)

        self.I = self.x.index          ## I 是 被评价决策单元的索引
        self.__modeldict = {}
        for i in self.I:
            # print(i)
            self.I0 = i                                                 ## I 是 被评价决策单元的数量

            self.__model__ = ConcreteModel()
            # Initialize sets
            self.__model__.I2 = Set(initialize= self.xref.index)                     ## I2 是 参考决策单元的数量
            self.__model__.K = Set(initialize=range(len(self.x.iloc[0])))          ## K 是投入个数
            self.__model__.L = Set(initialize=range(len(self.y.iloc[0])))          ## L 是产出个数 被评价单元和参考单元的K，L一样

            if type(self.b) != type(None):
                self.__model__.B = Set(initialize=range(len(self.b.iloc[0])))      ## B 是 非期望产出个数

            # Initialize variable
            self.__model__.thetax = Var(self.__model__.K,bounds=(0.0, None), doc='slack x')
            self.__model__.thetay = Var(self.__model__.L,bounds=(0.0, None), doc='slack y')
            if type(self.b) != type(None):
                self.__model__.thetab = Var(self.__model__.B,bounds=(0.0, None), doc='slack b')

            self.__model__.lamda = Var(self.__model__.I2, bounds=(0.0, None),within=Reals, doc='intensity variables')


            # Setup the objective function and constraints
            self.__model__.objective = Objective(rule=self.__objective_rule(), sense=maximize, doc='objective function')
            self.__model__.input = Constraint(self.__model__.K,  rule=self.__input_rule(), doc='input constraint')
            self.__model__.output = Constraint(self.__model__.L,  rule=self.__output_rule(), doc='output constraint')


            if type(self.b) != type(None):
                self.__model__.undesirable_output = Constraint(self.__model__.B, rule=self.__undesirable_output_rule(), \
                                                               doc='undesirable output constraint')
                self.__model__.mb = Constraint(self.__model__.B,  rule=self.__mb_rule(), \
                                                                doc='material balance constraint')
            if self.rts == RTS_VRS1:
                self.__model__.vrs = Constraint(rule=self.__vrs_rule(), doc='various return to scale rule')

            self.__modeldict[i] = self.__model__

        # Optimize model
    def __objective_rule(self):
        """Return the proper objective function"""
        def objective_rule(model):
            return sum(model.thetab[b]*1 for b in model.B)
        return objective_rule

    def __input_rule(self):
        """Return the proper input constraint"""
        def input_rule(model, k):
            return sum(model.lamda[i2] * self.xref.loc[i2,self.xcol[k]] for i2 in model.I2
                    ) + model.thetax[k] == self.x.loc[self.I0,self.xcol[k]]
        return input_rule

    def __output_rule(self):
        """Return the proper output constraint"""
        def output_rule(model, l):
            return sum(model.lamda[i2] * self.yref.loc[i2,self.ycol[l]] for i2 in model.I2
                    ) - model.thetay[l] == self.y.loc[self.I0,self.ycol[l]]
        return output_rule

    def __undesirable_output_rule(self):
        """Return the proper undesirable output constraint"""
        def undesirable_output_rule(model, b):
            return sum(model.lamda[i2] * self.bref.loc[i2,self.bcol[b]] for i2 in model.I2
                    ) + model.thetab[b] == self.b.loc[self.I0,self.bcol[b]]
        return undesirable_output_rule

    def __mb_rule(self):
        """Return the proper undesirable output constraint"""
        def mb_rule(model, b):
            return sum(self.sx[b][k] * model.thetax[k] for k in model.K) \
                    + sum(self.sy[b][l] * model.thetay[l] for l in model.L) \
                    == model.thetab[b]
        return mb_rule

    def __vrs_rule(self):
        def vrs_rule(model):
            return sum(model.lamda[ i2] for i2 in model.I2) == 1

        return vrs_rule

    def optimize(self,  solver=OPT_DEFAULT):
        """Optimize the function by requested method

        Args:
            solver (string): The solver chosen for optimization. It will optimize with default solver if OPT_DEFAULT is given.
        """
        # TODO(error/warning handling): Check problem status after optimization

        data2, obj, thetax, thetay, thetab, lamda = pd.DataFrame(), {}, {}, {}, {}, {}
        for ind, problem in self.__modeldict.items():
            _, data2.loc[ind, "optimization_status"] = tools.optimize_model(problem, ind, solver)

            if type(self.b) != type(None):
                obj[ind] = problem.objective()
                thetax[ind] = np.asarray(list(problem.thetax[:].value))
                thetay[ind] = np.asarray(list(problem.thetay[:].value))
                thetab[ind] = np.asarray(list(problem.thetab[:].value))
                lamda[ind] = np.asarray(list(problem.lamda[:].value))
            else:
                obj[ind] = problem.objective()
                thetax[ind] = np.asarray(list(problem.thetax[:].value))
                thetay[ind] = np.asarray(list(problem.thetay[:].value))
                lamda[ind] = np.asarray(list(problem.lamda[:].value))

                # print(list(problem.thetax[:].value ),list(problem.t[:].value ))
        obj = pd.DataFrame(obj, index=["obj"]).T
        thetax = pd.DataFrame(thetax).T
        thetax.columns = thetax.columns.map(lambda x: "Input" + str(x) + "'s slack")
        thetay = pd.DataFrame(thetay).T
        thetay.columns = thetay.columns.map(lambda y: "Output" + str(y) + "'s slack")

        theta_ = pd.concat([thetax, thetay], axis=1)
        if type(self.b) != type(None):
            thetab = pd.DataFrame(thetab).T
            thetab.columns = thetab.columns.map(lambda b: "Undesirable Output" + str(b) + "'s slack")
            theta_ = pd.concat([theta_, thetab], axis=1)

        lamda = pd.DataFrame(lamda).T
        lamda.columns = lamda.columns.map(lambda x: "lamda" + str(x))
        data3 = pd.concat([data2, obj], axis=1)
        data3 = pd.concat([data3, theta_], axis=1)
        # data3 = pd.concat([data3,lamda],axis=1)
        return data3

    def info(self, dmu = "all"):
        """Show the infomation of the lp model

        Args:
            dmu (string): The solver chosen for optimization. Default is "all".
        """
        if dmu =="all":
            for ind, problem in self.__modeldict.items():
                print(ind,"\n",problem.pprint())

        print(self.__modeldict[int(dmu)].pprint())


