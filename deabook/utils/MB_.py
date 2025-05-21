# import dependencies

from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, maximize, Constraint, Reals,PositiveReals
import numpy as np
import pandas as pd
from ..constant import CET_ADDI, CET_MULT, FUN_PROD, FUN_COST, RTS_CRS, RTS_VRS1, OPT_DEFAULT, OPT_LOCAL
import ast
from . import tools

class MB1111():
    """initial Group-VC-added CNLSZ (CNLSZ+G) model
    """

    def __init__(self, data, inputvars_np, inputvars_p,outputvars_np,\
                         outputvars_p, unoutputvars, sx, sy,  rts,level,baseindex,refindex):
        """CNLSZ+G model

        Args:
            y (float): output variable.
            x (float): input variables.
            z (float, optional): Contextual variable(s). Defaults to None.
            cutactive (float): active concavity constraint.
            cet (String, optional): CET_ADDI (additive composite error term) or CET_MULT (multiplicative composite error term). Defaults to CET_ADDI.
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS1.
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.data, self.inputvars_np, self.inputvars_p, self.outputvars_np, self.outputvars_p, \
                        self.unoutputvars, self.sx, self.sy, self.rts, \
                        self.level, self.baseindex, self.refindex = \
            data, inputvars_np, inputvars_p,outputvars_np,outputvars_p, \
                        unoutputvars, sx, sy,  rts,\
                        level,baseindex,refindex

        self.sx_np = np.array(self.sx)[:,0:len(self.inputvars_np)]
        self.sx_p = np.array(self.sx)[:,len(self.inputvars_np):]
        self.sy_np = np.array(self.sy)[:,0:len(self.outputvars_np)]
        self.sy_p = np.array(self.sy)[:,len(self.outputvars_np):]

        # print(self.inputvars_np, self.inputvars_p,self.outputvars_np,\
        #     self.outputvars_p,self.unoutputvars,self.sx, self.sy)


        if type(baseindex) != type(None):
            self.varname1=self.baseindex.split('=')[0]
            # print(self.baseindex)
            self.varvalue1=ast.literal_eval(self.baseindex.split('=')[1])
            self.y_p,self.y_np, self.x_p,self.x_np, self.b = \
                                            self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.outputvars_p
                                        ], self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.outputvars_np
                                        ],self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.inputvars_p
                                        ], self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.inputvars_np
                                        ], self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.unoutputvars
                                        ]

        else:
            self.y_p,self.y_np, self.x_p,self.x_np, self.b  = \
                                            self.data.loc[:, self.outputvars_p
                                        ], self.data.loc[:, self.outputvars_np
                                        ] , self.data.loc[:, self.inputvars_p
                                        ], self.data.loc[:, self.inputvars_np
                                        ], self.data.loc[:, self.unoutputvars
                                        ]


        if type(refindex) != type(None):
            self.varname=self.refindex.split('=')[0]
            self.varvalue=ast.literal_eval(self.refindex.split('=')[1])

            self.yref_p,self.yref_np, self.xref_p, self.xref_np, self.bref = \
                                         self.data.loc[self.data[self.varname].isin(self.varvalue), self.outputvars_p
                                    ], self.data.loc[self.data[self.varname].isin(self.varvalue), self.outputvars_np
                                    ] , self.data.loc[self.data[self.varname].isin(self.varvalue), self.inputvars_p
                                    ], self.data.loc[self.data[self.varname].isin(self.varvalue), self.inputvars_np
                                    ], self.data.loc[self.data[self.varname].isin(self.varvalue), self.unoutputvars ]
        else:
            self.yref_p, self.yref_np,self.xref_p, self.xref_np, self.bref = \
                                        self.data.loc[:, self.outputvars_p
                                    ], self.data.loc[:, self.outputvars_np
                                    ], self.data.loc[:, self.inputvars_p
                                    ], self.data.loc[:, self.inputvars_np
                                    ], self.data.loc[:, self.unoutputvars ]


        self.xcol_p = self.x_p.columns
        self.xcol_np = self.x_np.columns
        self.ycol_p = self.y_p.columns
        self.ycol_np = self.y_np.columns
        self.bcol = self.b.columns

        self.I = self.x_p.index          ## I 是 被评价决策单元的索引

        self.__modeldict = {}
        for i in self.I:
            # print(i)
            self.I0 = i                                                 ## I 是 被评价决策单元的数量

            self.__model__ = ConcreteModel()

            self.__model__.I2 = Set(initialize=self.xref_p.index)  ## I2 是 参考决策单元的数量

            self.__model__.Knp = Set(initialize=range(len(self.x_np.iloc[0])))  ## K 是投入个数
            self.__model__.Kp = Set(initialize=range(len(self.x_p.iloc[0])))  ## K 是投入个数
            self.__model__.Lnp = Set(initialize=range(len(self.y_np.iloc[0])))  ## L 是产出个数 被评价单元和参考单元的K，L一样

            self.__model__.Lp = Set(initialize=range(len(self.y_p.iloc[0])))  ## L 是产出个数 被评价单元和参考单元的K，L一样

            self.__model__.B = Set(initialize=range(len(self.b.iloc[0])))  ## B 是 非期望产出个数

            # Initialize variable
            self.__model__.thetax_p = Var(self.__model__.Kp, bounds=(0.0, None), doc='slack x_p')
            self.__model__.thetax_np = Var(self.__model__.Knp, bounds=(0.0, None), doc='slack x_np')
            self.__model__.thetay_p = Var(self.__model__.Lp, bounds=(0.0, None), doc='slack y_p')
            self.__model__.thetay_np = Var(self.__model__.Lnp, bounds=(0.0, None), doc='slack y_np')
            self.__model__.thetab = Var(self.__model__.B, bounds=(0.0, None), doc='slack b')

            self.__model__.objb = Var(self.__model__.B, bounds=(0.0, None), within=Reals, doc='object b')
            if self.level >= 2:
                self.__model__.objx_p = Var(self.__model__.Kp, bounds=(0.0, None), within=Reals, doc='object x_p')
            if self.level >= 3:
                self.__model__.objx_np = Var(self.__model__.Knp, bounds=(0.0, None), within=Reals, doc='object x_np')
            if (self.level >= 4) and (type(self.y_p) != type(None)):
                self.__model__.objy_p = Var(self.__model__.Lp, bounds=(0.0, None), within=Reals, doc='object y_p')
            if self.level >= 5:
                self.__model__.objy_np = Var(self.__model__.Lnp, bounds=(0.0, None), within=Reals, doc='object y_np')

            self.__model__.lamda = Var(self.__model__.I2, bounds=(0.0, None), within=Reals, doc='intensity variables')

            # Setup the objective function and constraints
            self.__model__.objective = Objective(rule=self.__objective_rule(), sense=minimize, doc='objective function')

            self.__model__.input_np = Constraint(self.__model__.Knp, rule=self.__input_np_rule(),
                                                 doc='input_np constraint')
            self.__model__.input_p = Constraint(self.__model__.Kp, rule=self.__input_p_rule(), doc='input_p constraint')
            self.__model__.output_np = Constraint(self.__model__.Lnp, rule=self.__output_np_rule(),
                                                  doc='output_np constraint')

            self.__model__.output_p = Constraint(self.__model__.Lp, rule=self.__output_p_rule(),
                                                     doc='output_p constraint')

            self.__model__.undesirable_output = Constraint(self.__model__.B, rule=self.__undesirable_output_rule(), \
                                                           doc='undesirable output constraint')
            self.__model__.mb = Constraint(self.__model__.B, rule=self.__mb_rule(), \
                                           doc='material balance constraint')
            if self.rts == RTS_VRS1:
                self.__model__.vrs = Constraint(rule=self.__vrs_rule(), doc='various return to scale rule')

            self.__modeldict[i] = self.__model__

        # Optimize model

    def __objective_rule(self):
        """Return the proper objective function"""

        def objective_rule(model):
            return sum(model.objb[b] * 1 for b in model.B)

        return objective_rule

    def __input_p_rule(self):
        """Return the proper input constraint"""
        if self.level < 2: # level = 1
            def input_p_rule(model, kp):
                return sum(model.lamda[i2] * self.xref_p.loc[i2, self.xcol_p[kp]] for i2 in model.I2
                           ) + model.thetax_p[kp] == self.x_p.loc[self.I0, self.xcol_p[kp]]
            return input_p_rule

        else:
            def input_p_rule(model, kp):
                return sum(model.lamda[i2] * self.xref_p.loc[i2, self.xcol_p[kp]] for i2 in model.I2
                           ) + model.thetax_p[kp] == model.objx_p[kp]
            return input_p_rule

    def __input_np_rule(self):
        """Return the proper input constraint"""
        if self.level < 3:
            def input_np_rule(model, knp):
                return sum(model.lamda[i2] * self.xref_np.loc[i2, self.xcol_np[knp]] for i2 in model.I2
                           ) + model.thetax_np[knp] == self.x_np.loc[self.I0, self.xcol_np[knp]]

            return input_np_rule
        else:
            def input_np_rule(model, knp):
                return sum(model.lamda[i2] * self.xref_np.loc[i2, self.xcol_np[knp]] for i2 in model.I2
                           ) + model.thetax_np[knp] == model.objx_np[knp]

            return input_np_rule

    def __output_p_rule(self):
        """Return the proper output constraint"""
        if self.level < 4:
            def output_p_rule(model, lp):
                return sum(model.lamda[i2] * self.yref_p.loc[i2, self.ycol_p[lp]] for i2 in model.I2
                           ) - model.thetay_p[lp] == self.y_p.loc[self.I0, self.ycol_p[lp]]
            return output_p_rule
        else:
            def output_p_rule(model, lp):
                return sum(model.lamda[i2] * self.yref_p.loc[i2, self.ycol_p[lp]] for i2 in model.I2
                           ) - model.thetay_p[lp] == model.objy_p[lp]
            return output_p_rule

    def __output_np_rule(self):
        """Return the proper output constraint"""
        if self.level < 5:
            def output_np_rule(model, lnp):
                return sum(model.lamda[i2] * self.yref_np.loc[i2, self.ycol_np[lnp]] for i2 in model.I2
                           ) - model.thetay_np[lnp] == self.y_np.loc[self.I0, self.ycol_np[lnp]]

            return output_np_rule
        else:
            def output_np_rule(model, lnp):
                return sum(model.lamda[i2] * self.yref_np.loc[i2, self.ycol_np[lnp]] for i2 in model.I2
                           ) - model.thetay_np[lnp] == model.objy_np[lnp]

            return output_np_rule

    def __undesirable_output_rule(self):
        """Return the proper undesirable output constraint"""

        def undesirable_output_rule(model, b):
            return sum(model.lamda[i2] * self.bref.loc[i2, self.bcol[b]] for i2 in model.I2
                       ) + model.thetab[b] == model.objb[b]
        return undesirable_output_rule

    def __mb_rule(self):
        """Return the proper undesirable output constraint"""
        if self.level ==1:
            def mb_rule(model, b):
                return sum(self.sx_p[b][kp] * model.thetax[kp] for kp in model.Kp) \
                       + sum(self.sy_p[b][lp] * model.thetay[lp] for lp in model.Lp)\
                    == self.b.loc[self.I0, self.bcol[b]] - model.objb[b]
            return mb_rule

        elif self.level ==2:
            def mb_rule(model, b):
                return sum(
                    self.sx_p[b][kp] * (self.x_p.loc[self.I0, self.xcol_p[kp]] - model.objx_p[kp]) for kp in model.Kp) \
                       + sum(self.sy_p[b][lp] * model.thetay[lp] for lp in model.Lp) \
                        == self.b.loc[self.I0, self.bcol[b]] - model.objb[b]
            return mb_rule

        elif self.level ==3:
            def mb_rule(model, b):
                return sum(
                    self.sx_p[b][kp] * (self.x_p.loc[self.I0, self.xcol_p[kp]] - model.objx_p[kp]) for kp in model.Kp) \
                       + sum(self.sy_p[b][lp] * model.thetay[lp] for lp in model.Lp) \
                        == self.b.loc[self.I0, self.bcol[b]] - model.objb[b]
            return mb_rule

        elif self.level ==4:
            def mb_rule(model, b):
                return sum(
                    self.sx_p[b][kp] * (self.x_p.loc[self.I0, self.xcol_p[kp]] - model.objx_p[kp]) for kp in model.Kp) \
                       + sum(self.sy_p[b][lp] * (model.objy_p[lp] - self.y_p.loc[self.I0, self.ycol_p[lp]]) for lp in
                             model.Lp) == self.b.loc[self.I0, self.bcol[b]] - model.objb[b]
            return mb_rule

        elif self.level ==5:
            def mb_rule(model, b):
                return sum(
                    self.sx_p[b][kp] * (self.x_p.loc[self.I0, self.xcol_p[kp]] - model.objx_p[kp]) for kp in model.Kp) \
                       + sum(self.sy_p[b][lp] * (model.objy_p[lp] - self.y_p.loc[self.I0, self.ycol_p[lp]]) for lp in
                             model.Lp) == self.b.loc[self.I0, self.bcol[b]] - model.objb[b]
            return mb_rule

    def __vrs_rule(self):
        def vrs_rule(model):
            return sum(model.lamda[i2] for i2 in model.I2) == 1

        return vrs_rule

    def optimize(self, solver=OPT_DEFAULT):
        """Optimize the function by requested method

        Args:
            solver (string): The solver chosen for optimization. It will optimize with default solver if OPT_DEFAULT is given.
        """
        # TODO(error/warning handling): Check problem status after optimization

        data2, obj, objb, = pd.DataFrame(), {}, {}
        for ind, problem in self.__modeldict.items():
            _, data2.loc[ind, "optimization_status"] = tools.optimize_model4(problem, ind, solver)

            if type(self.b) != type(None):
                obj[ind] = problem.objective()
                objb[ind], = np.asarray(list(problem.objb[:].value))

            else:
                obj[ind] = problem.objective()
                objb[ind], = np.asarray(list(problem.theta[:].value))


                # print(list(problem.thetax[:].value ),list(problem.t[:].value ))
        obj = pd.DataFrame(obj, index=["obj"]).T
        objb = pd.DataFrame(objb, index=["best of Undesirable"]).T

        theta_ = pd.concat([obj, objb], axis=1)
        data3 = pd.concat([data2,theta_],axis=1)
        return data3

    def info(self, dmu="all"):
        """Show the infomation of the lp model

        Args:
            dmu (string): The solver chosen for optimization. Default is "all".
        """
        if dmu == "all":
            for ind, problem in self.__modeldict.items():
                # print(ind, "\n", problem.pprint())
                pass
        # print(self.__modeldict[int(dmu)].pprint())


class MB1110():
    """ 所有期望产出不含污染物质
    """

    def __init__(self, data, inputvars_np, inputvars_p,outputvars_np,\
                          unoutputvars, sx, sy,  rts,level,baseindex,refindex):
        """CNLSZ+G model

        Args:
            y (float): output variable.
            x (float): input variables.
            z (float, optional): Contextual variable(s). Defaults to None.
            cutactive (float): active concavity constraint.
            cet (String, optional): CET_ADDI (additive composite error term) or CET_MULT (multiplicative composite error term). Defaults to CET_ADDI.
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS1.
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.data, self.inputvars_np, self.inputvars_p, self.outputvars_np, \
                        self.unoutputvars, self.sx, self.sy, self.rts, \
                        self.level, self.baseindex, self.refindex = \
            data, inputvars_np, inputvars_p,outputvars_np, \
                        unoutputvars, sx, sy,  rts,\
                        level,baseindex,refindex

        self.sx_np = np.array(self.sx)[:,0:len(self.inputvars_np)]
        self.sx_p = np.array(self.sx)[:,len(self.inputvars_np):]
        self.sy_np = np.array(self.sy)[:,0:len(self.outputvars_np)]

        # print(self.inputvars_np, self.inputvars_p,self.outputvars_np,\
            # self.unoutputvars,self.sx, self.sy)


        if type(baseindex) != type(None):
            self.varname1=self.baseindex.split('=')[0]
            # print(self.baseindex)
            self.varvalue1=ast.literal_eval(self.baseindex.split('=')[1])
            self.y_np, self.x_p,self.x_np, self.b = \
                                           self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.outputvars_np
                                        ],self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.inputvars_p
                                        ], self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.inputvars_np
                                        ], self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.unoutputvars
                                        ]

        else:
            self.y_np, self.x_p,self.x_np, self.b  = \
                                           self.data.loc[:, self.outputvars_np
                                        ] , self.data.loc[:, self.inputvars_p
                                        ], self.data.loc[:, self.inputvars_np
                                        ], self.data.loc[:, self.unoutputvars
                                        ]


        if type(refindex) != type(None):
            self.varname=self.refindex.split('=')[0]
            self.varvalue=ast.literal_eval(self.refindex.split('=')[1])

            self.yref_np, self.xref_p, self.xref_np, self.bref = \
                                        self.data.loc[self.data[self.varname].isin(self.varvalue), self.outputvars_np
                                    ] , self.data.loc[self.data[self.varname].isin(self.varvalue), self.inputvars_p
                                    ], self.data.loc[self.data[self.varname].isin(self.varvalue), self.inputvars_np
                                    ], self.data.loc[self.data[self.varname].isin(self.varvalue), self.unoutputvars ]
        else:
            self.yref_np,self.xref_p, self.xref_np, self.bref = \
                                         self.data.loc[:, self.outputvars_np
                                    ], self.data.loc[:, self.inputvars_p
                                    ], self.data.loc[:, self.inputvars_np
                                    ], self.data.loc[:, self.unoutputvars ]


        self.xcol_p = self.x_p.columns
        self.xcol_np = self.x_np.columns
        # self.ycol_p = self.y_p.columns
        self.ycol_np = self.y_np.columns
        self.bcol = self.b.columns

        self.I = self.x_p.index          ## I 是 被评价决策单元的索引

        self.__modeldict = {}
        for i in self.I:
            # print(i)
            self.I0 = i                                                 ## I 是 被评价决策单元的数量

            self.__model__ = ConcreteModel()

            self.__model__.I2 = Set(initialize=self.xref_p.index)  ## I2 是 参考决策单元的数量

            self.__model__.Knp = Set(initialize=range(len(self.x_np.iloc[0])))  ## K 是投入个数
            self.__model__.Kp = Set(initialize=range(len(self.x_p.iloc[0])))  ## K 是投入个数
            self.__model__.Lnp = Set(initialize=range(len(self.y_np.iloc[0])))  ## L 是产出个数 被评价单元和参考单元的K，L一样

            # self.__model__.Lp = Set(initialize=range(len(self.y_p.iloc[0])))  ## L 是产出个数 被评价单元和参考单元的K，L一样

            self.__model__.B = Set(initialize=range(len(self.b.iloc[0])))  ## B 是 非期望产出个数

            # Initialize variable
            self.__model__.thetax_p = Var(self.__model__.Kp, bounds=(0.0, None), doc='slack x_p')
            self.__model__.thetax_np = Var(self.__model__.Knp, bounds=(0.0, None), doc='slack x_np')
            # self.__model__.thetay_p = Var(self.__model__.Lp, bounds=(0.0, None), doc='slack y_p')
            self.__model__.thetay_np = Var(self.__model__.Lnp, bounds=(0.0, None), doc='slack y_np')
            self.__model__.thetab = Var(self.__model__.B, bounds=(0.0, None), doc='slack b')

            self.__model__.objb = Var(self.__model__.B, bounds=(0.0, None), within=Reals, doc='object b')
            if self.level >= 2:
                self.__model__.objx_p = Var(self.__model__.Kp, bounds=(0.0, None), within=Reals, doc='object x_p')
            if self.level >= 3:
                self.__model__.objx_np = Var(self.__model__.Knp, bounds=(0.0, None), within=Reals, doc='object x_np')
            if (self.level >= 4) :
                self.__model__.objy_np = Var(self.__model__.Lnp, bounds=(0.0, None), within=Reals, doc='object y_np')

            self.__model__.lamda = Var(self.__model__.I2, bounds=(0.0, None), within=Reals, doc='intensity variables')

            # Setup the objective function and constraints
            self.__model__.objective = Objective(rule=self.__objective_rule(), sense=minimize, doc='objective function')

            self.__model__.input_np = Constraint(self.__model__.Knp, rule=self.__input_np_rule(),
                                                 doc='input_np constraint')
            self.__model__.input_p = Constraint(self.__model__.Kp, rule=self.__input_p_rule(), doc='input_p constraint')
            self.__model__.output_np = Constraint(self.__model__.Lnp, rule=self.__output_np_rule(),
                                                  doc='output_np constraint')



            self.__model__.undesirable_output = Constraint(self.__model__.B, rule=self.__undesirable_output_rule(), \
                                                           doc='undesirable output constraint')
            self.__model__.mb = Constraint(self.__model__.B, rule=self.__mb_rule(), \
                                           doc='material balance constraint')
            if self.rts == RTS_VRS1:
                self.__model__.vrs = Constraint(rule=self.__vrs_rule(), doc='various return to scale rule')

            self.__modeldict[i] = self.__model__

        # Optimize model

    def __objective_rule(self):
        """Return the proper objective function"""

        def objective_rule(model):
            return sum(model.objb[b] * 1 for b in model.B)

        return objective_rule

    def __input_p_rule(self):
        """Return the proper input constraint"""
        if self.level < 2: # level = 1
            def input_p_rule(model, kp):
                return sum(model.lamda[i2] * self.xref_p.loc[i2, self.xcol_p[kp]] for i2 in model.I2
                           ) + model.thetax_p[kp] == self.x_p.loc[self.I0, self.xcol_p[kp]]
            return input_p_rule

        else:
            def input_p_rule(model, kp):
                return sum(model.lamda[i2] * self.xref_p.loc[i2, self.xcol_p[kp]] for i2 in model.I2
                           ) + model.thetax_p[kp] == model.objx_p[kp]
            return input_p_rule

    def __input_np_rule(self):
        """Return the proper input constraint"""
        if self.level < 3:
            def input_np_rule(model, knp):
                return sum(model.lamda[i2] * self.xref_np.loc[i2, self.xcol_np[knp]] for i2 in model.I2
                           ) + model.thetax_np[knp] == self.x_np.loc[self.I0, self.xcol_np[knp]]

            return input_np_rule
        else:
            def input_np_rule(model, knp):
                return sum(model.lamda[i2] * self.xref_np.loc[i2, self.xcol_np[knp]] for i2 in model.I2
                           ) + model.thetax_np[knp] == model.objx_np[knp]

            return input_np_rule


    def __output_np_rule(self):
        """Return the proper output constraint"""
        if self.level < 4:
            def output_np_rule(model, lnp):
                return sum(model.lamda[i2] * self.yref_np.loc[i2, self.ycol_np[lnp]] for i2 in model.I2
                           ) - model.thetay_np[lnp] == self.y_np.loc[self.I0, self.ycol_np[lnp]]
            return output_np_rule

        else:
            def output_np_rule(model, lnp):
                return sum(model.lamda[i2] * self.yref_np.loc[i2, self.ycol_np[lnp]] for i2 in model.I2
                           ) - model.thetay_np[lnp] == model.objy_np[lnp]
            return output_np_rule

    def __undesirable_output_rule(self):
        """Return the proper undesirable output constraint"""
        def undesirable_output_rule(model, b):
            return sum(model.lamda[i2] * self.bref.loc[i2, self.bcol[b]] for i2 in model.I2
                       ) + model.thetab[b] == model.objb[b]
        return undesirable_output_rule

    def __mb_rule(self):
        """Return the proper undesirable output constraint"""
        if self.level ==1:
            def mb_rule(model, b):
                return sum(self.sx_p[b][kp] * model.thetax_p[kp] for kp in model.Kp) \
                               == self.b.loc[self.I0, self.bcol[b]] - model.objb[b]
            return mb_rule

        elif self.level ==2:
            def mb_rule(model, b):
                return sum(
                    self.sx_p[b][kp] * (self.x_p.loc[self.I0, self.xcol_p[kp]] - model.objx_p[kp]) for kp in model.Kp) \
                               == self.b.loc[self.I0, self.bcol[b]] - model.objb[b]
            return mb_rule

        elif self.level ==3:
            def mb_rule(model, b):
                return sum(
                    self.sx_p[b][kp] * (self.x_p.loc[self.I0, self.xcol_p[kp]] - model.objx_p[kp]) for kp in model.Kp) \
                               == self.b.loc[self.I0, self.bcol[b]] - model.objb[b]
            return mb_rule

        elif self.level ==4:
            def mb_rule(model, b):
                return sum(
                    self.sx_p[b][kp] * (self.x_p.loc[self.I0, self.xcol_p[kp]] - model.objx_p[kp]) for kp in model.Kp) \
                                == self.b.loc[self.I0, self.bcol[b]] - model.objb[b]
            return mb_rule


    def __vrs_rule(self):
        def vrs_rule(model):
            return sum(model.lamda[i2] for i2 in model.I2) == 1

        return vrs_rule

    def optimize(self, solver=OPT_DEFAULT):
        """Optimize the function by requested method

        Args:
            solver (string): The solver chosen for optimization. It will optimize with default solver if OPT_DEFAULT is given.
        """
        # TODO(error/warning handling): Check problem status after optimization

        data2, obj, objb, = pd.DataFrame(), {}, {}
        for ind, problem in self.__modeldict.items():
            _, data2.loc[ind, "optimization_status"] = tools.optimize_model4(problem, ind, solver)

            if type(self.b) != type(None):
                obj[ind] = problem.objective()
                objb[ind], = np.asarray(list(problem.objb[:].value))

            else:
                obj[ind] = problem.objective()
                objb[ind], = np.asarray(list(problem.theta[:].value))


                # print(list(problem.thetax[:].value ),list(problem.t[:].value ))
        obj = pd.DataFrame(obj, index=["obj"]).T
        objb = pd.DataFrame(objb, index=["best of Undesirable"]).T

        theta_ = pd.concat([obj, objb], axis=1)
        data3 = pd.concat([data2,theta_],axis=1)
        return data3

    def info(self, dmu="all"):
        """Show the infomation of the lp model

        Args:
            dmu (string): The solver chosen for optimization. Default is "all".
        """
        if dmu == "all":
            for ind, problem in self.__modeldict.items():
                # print(ind, "\n", problem.pprint())
                pass
        # print(self.__modeldict[int(dmu)].pprint())


class MB1101():
    """ 期望产出都包含污染物质
    """

    def __init__(self, data, inputvars_np, inputvars_p,outputvars_p,\
                        unoutputvars, sx, sy,  rts,level,baseindex,refindex):
        """CNLSZ+G model

        Args:
            y (float): output variable.
            x (float): input variables.
            z (float, optional): Contextual variable(s). Defaults to None.
            cutactive (float): active concavity constraint.
            cet (String, optional): CET_ADDI (additive composite error term) or CET_MULT (multiplicative composite error term). Defaults to CET_ADDI.
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS1.
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.data, self.inputvars_np, self.inputvars_p, self.outputvars_p, \
                        self.unoutputvars, self.sx, self.sy, self.rts, \
                        self.level, self.baseindex, self.refindex = \
            data, inputvars_np, inputvars_p,outputvars_p, \
                        unoutputvars, sx, sy,  rts,\
                        level,baseindex,refindex

        self.sx_np = np.array(self.sx)[:,0:len(self.inputvars_np)]
        self.sx_p = np.array(self.sx)[:,len(self.inputvars_np):]

        self.sy_p = np.array(self.sy)[:,0:len(self.outputvars_p)]

        # print(self.inputvars_np, self.inputvars_p,self.outputvars_p,self.unoutputvars,self.sx, self.sy)
        # print("sssssssssss",self.sy_p)

        if type(baseindex) != type(None):
            self.varname1=self.baseindex.split('=')[0]
            # print(self.baseindex)
            self.varvalue1=ast.literal_eval(self.baseindex.split('=')[1])
            self.y_p, self.x_p,self.x_np, self.b = \
                                            self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.outputvars_p
                                        ],self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.inputvars_p
                                        ], self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.inputvars_np
                                        ], self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.unoutputvars
                                        ]

        else:
            self.y_p, self.x_p,self.x_np, self.b  = \
                                            self.data.loc[:, self.outputvars_p
                                        ], self.data.loc[:, self.inputvars_p
                                        ], self.data.loc[:, self.inputvars_np
                                        ], self.data.loc[:, self.unoutputvars
                                        ]


        if type(refindex) != type(None):
            self.varname=self.refindex.split('=')[0]
            self.varvalue=ast.literal_eval(self.refindex.split('=')[1])

            self.yref_p, self.xref_p, self.xref_np, self.bref = \
                                         self.data.loc[self.data[self.varname].isin(self.varvalue), self.outputvars_p
                                    ], self.data.loc[self.data[self.varname].isin(self.varvalue), self.inputvars_p
                                    ], self.data.loc[self.data[self.varname].isin(self.varvalue), self.inputvars_np
                                    ], self.data.loc[self.data[self.varname].isin(self.varvalue), self.unoutputvars ]
        else:
            self.yref_p,self.xref_p, self.xref_np, self.bref = \
                                        self.data.loc[:, self.outputvars_p
                                    ], self.data.loc[:, self.inputvars_p
                                    ], self.data.loc[:, self.inputvars_np
                                    ], self.data.loc[:, self.unoutputvars ]


        self.xcol_p = self.x_p.columns
        self.xcol_np = self.x_np.columns
        self.ycol_p = self.y_p.columns
        self.bcol = self.b.columns

        self.I = self.x_p.index          ## I 是 被评价决策单元的索引

        self.__modeldict = {}
        for i in self.I:
            # print(i)
            self.I0 = i                                                 ## I 是 被评价决策单元的数量

            self.__model__ = ConcreteModel()

            self.__model__.I2 = Set(initialize=self.xref_p.index)  ## I2 是 参考决策单元的数量

            self.__model__.Knp = Set(initialize=range(len(self.x_np.iloc[0])))  ## K 是投入个数
            self.__model__.Kp = Set(initialize=range(len(self.x_p.iloc[0])))  ## K 是投入个数

            self.__model__.Lp = Set(initialize=range(len(self.y_p.iloc[0])))  ## L 是产出个数 被评价单元和参考单元的K，L一样

            self.__model__.B = Set(initialize=range(len(self.b.iloc[0])))  ## B 是 非期望产出个数

            # Initialize variable
            self.__model__.thetax_p = Var(self.__model__.Kp, bounds=(0.0, None), doc='slack x_p')
            self.__model__.thetax_np = Var(self.__model__.Knp, bounds=(0.0, None), doc='slack x_np')
            self.__model__.thetay_p = Var(self.__model__.Lp, bounds=(0.0, None), doc='slack y_p')
            self.__model__.thetab = Var(self.__model__.B, bounds=(0.0, None), doc='slack b')

            self.__model__.objb = Var(self.__model__.B, bounds=(0.0, None), within=Reals, doc='object b')
            if self.level >= 2:
                self.__model__.objx_p = Var(self.__model__.Kp, bounds=(0.0, None), within=Reals, doc='object x_p')
            if self.level >= 3:
                self.__model__.objx_np = Var(self.__model__.Knp, bounds=(0.0, None), within=Reals, doc='object x_np')
            if (self.level >= 4) and (type(self.y_p) != type(None)):
                self.__model__.objy_p = Var(self.__model__.Lp, bounds=(0.0, None), within=Reals, doc='object y_p')

            self.__model__.lamda = Var(self.__model__.I2, bounds=(0.0, None), within=Reals, doc='intensity variables')

            # Setup the objective function and constraints
            self.__model__.objective = Objective(rule=self.__objective_rule(), sense=minimize, doc='objective function')

            self.__model__.input_np = Constraint(self.__model__.Knp, rule=self.__input_np_rule(),
                                                 doc='input_np constraint')
            self.__model__.input_p = Constraint(self.__model__.Kp, rule=self.__input_p_rule(), doc='input_p constraint')

            self.__model__.output_p = Constraint(self.__model__.Lp, rule=self.__output_p_rule(),
                                                     doc='output_p constraint')

            self.__model__.undesirable_output = Constraint(self.__model__.B, rule=self.__undesirable_output_rule(), \
                                                           doc='undesirable output constraint')
            self.__model__.mb = Constraint(self.__model__.B, rule=self.__mb_rule(), \
                                           doc='material balance constraint')
            if self.rts == RTS_VRS1:
                self.__model__.vrs = Constraint(rule=self.__vrs_rule(), doc='various return to scale rule')

            self.__modeldict[i] = self.__model__

        # Optimize model

    def __objective_rule(self):
        """Return the proper objective function"""

        def objective_rule(model):
            return sum(model.objb[b] * 1 for b in model.B)

        return objective_rule

    def __input_p_rule(self):
        """Return the proper input constraint"""
        if self.level < 2: # level = 1
            def input_p_rule(model, kp):
                return sum(model.lamda[i2] * self.xref_p.loc[i2, self.xcol_p[kp]] for i2 in model.I2
                           ) + model.thetax_p[kp] == self.x_p.loc[self.I0, self.xcol_p[kp]]
            return input_p_rule

        else:
            def input_p_rule(model, kp):
                return sum(model.lamda[i2] * self.xref_p.loc[i2, self.xcol_p[kp]] for i2 in model.I2
                           ) + model.thetax_p[kp] == model.objx_p[kp]
            return input_p_rule

    def __input_np_rule(self):
        """Return the proper input constraint"""
        if self.level < 3:
            def input_np_rule(model, knp):
                return sum(model.lamda[i2] * self.xref_np.loc[i2, self.xcol_np[knp]] for i2 in model.I2
                           ) + model.thetax_np[knp] == self.x_np.loc[self.I0, self.xcol_np[knp]]

            return input_np_rule
        else:
            def input_np_rule(model, knp):
                return sum(model.lamda[i2] * self.xref_np.loc[i2, self.xcol_np[knp]] for i2 in model.I2
                           ) + model.thetax_np[knp] == model.objx_np[knp]

            return input_np_rule

    def __output_p_rule(self):
        """Return the proper output constraint"""
        if self.level < 4:
            def output_p_rule(model, lp):
                return sum(model.lamda[i2] * self.yref_p.loc[i2, self.ycol_p[lp]] for i2 in model.I2
                           ) - model.thetay_p[lp] == self.y_p.loc[self.I0, self.ycol_p[lp]]
            return output_p_rule
        else:
            def output_p_rule(model, lp):
                return sum(model.lamda[i2] * self.yref_p.loc[i2, self.ycol_p[lp]] for i2 in model.I2
                           ) - model.thetay_p[lp] == model.objy_p[lp]
            return output_p_rule


    def __undesirable_output_rule(self):
        """Return the proper undesirable output constraint"""

        def undesirable_output_rule(model, b):
            return sum(model.lamda[i2] * self.bref.loc[i2, self.bcol[b]] for i2 in model.I2
                       ) + model.thetab[b] == model.objb[b]
        return undesirable_output_rule

    def __mb_rule(self):
        """Return the proper undesirable output constraint"""
        if self.level ==1:
            def mb_rule(model, b):
                return sum(self.sx_p[b][kp] * model.thetax[kp] for kp in model.Kp) \
                       + sum(self.sy_p[b][lp] * model.thetay[lp] for lp in model.Lp) if type(self.sy_p) != type(
                    None) else 0 \
                               == self.b.loc[self.I0, self.bcol[b]] - model.objb[b]
            return mb_rule

        elif self.level ==2:
            def mb_rule(model, b):
                return sum(
                    self.sx_p[b][kp] * (self.x_p.loc[self.I0, self.xcol_p[kp]] - model.objx_p[kp]) for kp in model.Kp) \
                       + sum(self.sy_p[b][lp] * model.thetay[lp] for lp in model.Lp) if type(self.sy_p) != type(
                    None) else 0 \
                               == self.b.loc[self.I0, self.bcol[b]] - model.objb[b]
            return mb_rule

        elif self.level ==3:
            def mb_rule(model, b):
                return sum(
                    self.sx_p[b][kp] * (self.x_p.loc[self.I0, self.xcol_p[kp]] - model.objx_p[kp]) for kp in model.Kp) \
                       + sum(self.sy_p[b][lp] * model.thetay[lp] for lp in model.Lp) if type(self.sy_p) != type(
                    None) else 0 \
                               == self.b.loc[self.I0, self.bcol[b]] - model.objb[b]
            return mb_rule

        elif self.level ==4:
            def mb_rule(model, b):
                return sum(
                    self.sx_p[b][kp] * (self.x_p.loc[self.I0, self.xcol_p[kp]] - model.objx_p[kp]) for kp in model.Kp) \
              + sum(self.sy_p[b][lp] * (model.objy_p[lp] - self.y_p.loc[self.I0, self.ycol_p[lp]]) for lp in model.Lp) \
                                     == self.b.loc[self.I0, self.bcol[b]] - model.objb[b]
            return mb_rule



    def __vrs_rule(self):
        def vrs_rule(model):
            return sum(model.lamda[i2] for i2 in model.I2) == 1

        return vrs_rule

    def optimize(self, solver=OPT_DEFAULT):
        """Optimize the function by requested method

        Args:
            solver (string): The solver chosen for optimization. It will optimize with default solver if OPT_DEFAULT is given.
        """
        # TODO(error/warning handling): Check problem status after optimization

        data2, obj, objb, = pd.DataFrame(), {}, {}
        for ind, problem in self.__modeldict.items():
            _, data2.loc[ind, "optimization_status"] = tools.optimize_model4(problem, ind, solver)

            if type(self.b) != type(None):
                obj[ind] = problem.objective()
                objb[ind], = np.asarray(list(problem.objb[:].value))

            else:
                obj[ind] = problem.objective()
                objb[ind], = np.asarray(list(problem.theta[:].value))


                # print(list(problem.thetax[:].value ),list(problem.t[:].value ))
        obj = pd.DataFrame(obj, index=["obj"]).T
        objb = pd.DataFrame(objb, index=["best of Undesirable"]).T

        theta_ = pd.concat([obj, objb], axis=1)
        data3 = pd.concat([data2,theta_],axis=1)
        return data3

    def info(self, dmu="all"):
        """Show the infomation of the lp model

        Args:
            dmu (string): The solver chosen for optimization. Default is "all".
        """
        if dmu == "all":
            for ind, problem in self.__modeldict.items():
                # print(ind, "\n", problem.pprint())
                pass
        # print(self.__modeldict[int(dmu)].pprint())



class MB1100():
    """ 没有期望产出
    """

    def __init__(self, data, inputvars_np, inputvars_p,\
                        unoutputvars, sx, sy,  rts,level,baseindex,refindex):
        """CNLSZ+G model

        Args:
            y (float): output variable.
            x (float): input variables.
            z (float, optional): Contextual variable(s). Defaults to None.
            cutactive (float): active concavity constraint.
            cet (String, optional): CET_ADDI (additive composite error term) or CET_MULT (multiplicative composite error term). Defaults to CET_ADDI.
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS1 (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS1.
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.data, self.inputvars_np, self.inputvars_p, \
                        self.unoutputvars, self.sx, self.sy, self.rts, \
                        self.level, self.baseindex, self.refindex = \
            data, inputvars_np, inputvars_p, \
                        unoutputvars, sx, sy,  rts,\
                        level,baseindex,refindex

        self.sx_np = np.array(self.sx)[:,0:len(self.inputvars_np)]
        self.sx_p = np.array(self.sx)[:,len(self.inputvars_np):]

        # print(self.inputvars_np, self.inputvars_p,self.unoutputvars,self.sx, self.sy)


        if type(baseindex) != type(None):
            self.varname1=self.baseindex.split('=')[0]
            # print(self.baseindex)
            self.varvalue1=ast.literal_eval(self.baseindex.split('=')[1])
            self.x_p,self.x_np, self.b = \
                                        self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.inputvars_p
                                        ], self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.inputvars_np
                                        ], self.data.loc[self.data[self.varname1].isin(self.varvalue1), self.unoutputvars
                                        ]

        else:
            self.x_p,self.x_np, self.b  = \
                                        self.data.loc[:, self.inputvars_p
                                        ], self.data.loc[:, self.inputvars_np
                                        ], self.data.loc[:, self.unoutputvars
                                        ]


        if type(refindex) != type(None):
            self.varname=self.refindex.split('=')[0]
            self.varvalue=ast.literal_eval(self.refindex.split('=')[1])

            self.xref_p, self.xref_np, self.bref = \
                                        self.data.loc[self.data[self.varname].isin(self.varvalue), self.inputvars_p
                                    ], self.data.loc[self.data[self.varname].isin(self.varvalue), self.inputvars_np
                                    ], self.data.loc[self.data[self.varname].isin(self.varvalue), self.unoutputvars ]
        else:
            self.xref_p, self.xref_np, self.bref = \
                                     self.data.loc[:, self.inputvars_p
                                    ], self.data.loc[:, self.inputvars_np
                                    ], self.data.loc[:, self.unoutputvars ]


        self.xcol_p = self.x_p.columns
        self.xcol_np = self.x_np.columns
        self.bcol = self.b.columns

        self.I = self.x_p.index          ## I 是 被评价决策单元的索引

        self.__modeldict = {}
        for i in self.I:
            # print(i)
            self.I0 = i                                                 ## I 是 被评价决策单元的数量

            self.__model__ = ConcreteModel()

            self.__model__.I2 = Set(initialize=self.xref_p.index)  ## I2 是 参考决策单元的数量

            self.__model__.Knp = Set(initialize=range(len(self.x_np.iloc[0])))  ## K 是投入个数
            self.__model__.Kp = Set(initialize=range(len(self.x_p.iloc[0])))  ## K 是投入个数


            self.__model__.B = Set(initialize=range(len(self.b.iloc[0])))  ## B 是 非期望产出个数

            # Initialize variable
            self.__model__.thetax_p = Var(self.__model__.Kp, bounds=(0.0, None), doc='slack x_p')
            self.__model__.thetax_np = Var(self.__model__.Knp, bounds=(0.0, None), doc='slack x_np')
            self.__model__.thetab = Var(self.__model__.B, bounds=(0.0, None), doc='slack b')

            self.__model__.objb = Var(self.__model__.B, bounds=(0.0, None), within=Reals, doc='object b')
            if self.level >= 2:
                self.__model__.objx_p = Var(self.__model__.Kp, bounds=(0.0, None), within=Reals, doc='object x_p')
            if self.level >= 3:
                self.__model__.objx_np = Var(self.__model__.Knp, bounds=(0.0, None), within=Reals, doc='object x_np')

            self.__model__.lamda = Var(self.__model__.I2, bounds=(0.0, None), within=Reals, doc='intensity variables')

            # Setup the objective function and constraints
            self.__model__.objective = Objective(rule=self.__objective_rule(), sense=minimize, doc='objective function')

            self.__model__.input_np = Constraint(self.__model__.Knp, rule=self.__input_np_rule(),
                                                 doc='input_np constraint')
            self.__model__.input_p = Constraint(self.__model__.Kp, rule=self.__input_p_rule(), doc='input_p constraint')



            self.__model__.undesirable_output = Constraint(self.__model__.B, rule=self.__undesirable_output_rule(), \
                                                           doc='undesirable output constraint')
            self.__model__.mb = Constraint(self.__model__.B, rule=self.__mb_rule(), \
                                           doc='material balance constraint')
            if self.rts == RTS_VRS1:
                self.__model__.vrs = Constraint(rule=self.__vrs_rule(), doc='various return to scale rule')

            self.__modeldict[i] = self.__model__

        # Optimize model

    def __objective_rule(self):
        """Return the proper objective function"""

        def objective_rule(model):
            return sum(model.objb[b] * 1 for b in model.B)

        return objective_rule

    def __input_p_rule(self):
        """Return the proper input constraint"""
        if self.level < 2: # level = 1
            def input_p_rule(model, kp):
                return sum(model.lamda[i2] * self.xref_p.loc[i2, self.xcol_p[kp]] for i2 in model.I2
                           ) + model.thetax_p[kp] == self.x_p.loc[self.I0, self.xcol_p[kp]]
            return input_p_rule

        else:
            def input_p_rule(model, kp):
                return sum(model.lamda[i2] * self.xref_p.loc[i2, self.xcol_p[kp]] for i2 in model.I2
                           ) + model.thetax_p[kp] == model.objx_p[kp]
            return input_p_rule

    def __input_np_rule(self):
        """Return the proper input constraint"""
        if self.level < 3:
            def input_np_rule(model, knp):
                return sum(model.lamda[i2] * self.xref_np.loc[i2, self.xcol_np[knp]] for i2 in model.I2
                           ) + model.thetax_np[knp] == self.x_np.loc[self.I0, self.xcol_np[knp]]

            return input_np_rule
        else:
            def input_np_rule(model, knp):
                return sum(model.lamda[i2] * self.xref_np.loc[i2, self.xcol_np[knp]] for i2 in model.I2
                           ) + model.thetax_np[knp] == model.objx_np[knp]

            return input_np_rule



    def __undesirable_output_rule(self):
        """Return the proper undesirable output constraint"""

        def undesirable_output_rule(model, b):
            return sum(model.lamda[i2] * self.bref.loc[i2, self.bcol[b]] for i2 in model.I2
                       ) + model.thetab[b] == model.objb[b]
        return undesirable_output_rule

    def __mb_rule(self):
        """Return the proper undesirable output constraint"""
        if self.level ==1:
            def mb_rule(model, b):
                return sum(self.sx_p[b][kp] * model.thetax[kp] for kp in model.Kp) \
                      == self.b.loc[self.I0, self.bcol[b]] - model.objb[b]
            return mb_rule

        elif self.level ==2:
            def mb_rule(model, b):
                return sum(
                    self.sx_p[b][kp] * (self.x_p.loc[self.I0, self.xcol_p[kp]] - model.objx_p[kp]) for kp in model.Kp) \
                        == self.b.loc[self.I0, self.bcol[b]] - model.objb[b]
            return mb_rule

        elif self.level ==3:
            def mb_rule(model, b):
                return sum(
                    self.sx_p[b][kp] * (self.x_p.loc[self.I0, self.xcol_p[kp]] - model.objx_p[kp]) for kp in model.Kp) \
                        == self.b.loc[self.I0, self.bcol[b]] - model.objb[b]
            return mb_rule





    def __vrs_rule(self):
        def vrs_rule(model):
            return sum(model.lamda[i2] for i2 in model.I2) == 1

        return vrs_rule

    def optimize(self, solver=OPT_DEFAULT):
        """Optimize the function by requested method

        Args:
            solver (string): The solver chosen for optimization. It will optimize with default solver if OPT_DEFAULT is given.
        """
        # TODO(error/warning handling): Check problem status after optimization

        data2, obj, objb, = pd.DataFrame(), {}, {}
        for ind, problem in self.__modeldict.items():
            _, data2.loc[ind, "optimization_status"] = tools.optimize_model4(problem, ind, solver)

            if type(self.b) != type(None):
                obj[ind] = problem.objective()
                objb[ind], = np.asarray(list(problem.objb[:].value))

            else:
                obj[ind] = problem.objective()
                objb[ind], = np.asarray(list(problem.theta[:].value))


                # print(list(problem.thetax[:].value ),list(problem.t[:].value ))
        obj = pd.DataFrame(obj, index=["obj"]).T
        objb = pd.DataFrame(objb, index=["best of Undesirable"]).T

        theta_ = pd.concat([obj, objb], axis=1)
        data3 = pd.concat([data2,theta_],axis=1)
        return data3

    def info(self, dmu="all"):
        """Show the infomation of the lp model

        Args:
            dmu (string): The solver chosen for optimization. Default is "all".
        """
        if dmu == "all":
            for ind, problem in self.__modeldict.items():
                # print(ind, "\n", problem.pprint())
                pass
        # print(self.__modeldict[int(dmu)].pprint())


