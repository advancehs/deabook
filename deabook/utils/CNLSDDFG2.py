# import dependencies
from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, Constraint, log
from pyomo.core.expr.numvalue import NumericValue
import numpy as np
import pandas as pd
from ..constant import CET_ADDI, FUN_PROD, FUN_COST, RTS_CRS, RTS_VRS1, OPT_DEFAULT, OPT_LOCAL
from .tools import optimize_model
from . import CNLSDDFZG2, CNLSDDFG1

class CNLSDDFG2(CNLSDDFZG2.CNLSDDFZG2, CNLSDDFG1.CNLSDDFG1):
    """initial Group-VC-added CNLSDDF (CNLSDDF+G) model
    """

    def __init__(self, y, x, b, cutactive, active, gy=[1], gx=[1], gb=[1],  fun=FUN_PROD, rts=RTS_VRS1):
        """CNLSDDF+G model

        Args:
            y (float): output variable.
            x (float): input variables.
            b (float, optional): undersiable variables.
            cutactive (float, optional): active concavity constraint.
            active (float or ndarray ): violated concavity constraint.
            gy (list, optional): output directional vector. Defaults to [1].
            gx (list, optional): input directional vector. Defaults to [1].
            gb (list, optional): undesirable output directional vector. Defaults to [1].
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS.
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.x,self.y,self.b = x, y, b
        self.gy, self.gx, self.gb = gy,gx,gb
        self.fun = fun
        self.rts = rts

        self.cutactive = cutactive
        self.active = active

        # Initialize the CNLS model
        self.__model__ = ConcreteModel()

        # Initialize the sets
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))
        self.__model__.K = Set(initialize=range(len(self.y[0])))

        # Initialize the variables
        self.__model__.alpha = Var(self.__model__.I, doc='alpha')
        self.__model__.beta = Var(self.__model__.I,
                                  self.__model__.J,
                                  bounds=(0.0, None),
                                  doc='beta')
        self.__model__.gamma = Var(self.__model__.I,
                                  self.__model__.K,
                                  bounds=(0.0, None),
                                  doc='gamma')


        if type(self.b) != type(None):
            self.__model__.L = Set(initialize=range(len(self.b[0])))
            self.__model__.delta = Var(
                self.__model__.I, self.__model__.L, bounds=(0.0, None), doc='delta')


        self.__model__.epsilon = Var(self.__model__.I, doc='residual')
        self.__model__.frontier = Var(self.__model__.I,
                                      bounds=(0.0, None),
                                      doc='estimated frontier')

        # Setup the objective function and constraints
        self.__model__.objective = Objective(rule=self._CNLSDDFG1__objective_rule(),
                                             sense=minimize,
                                             doc='objective function')
        self.__model__.regression_rule = Constraint(self.__model__.I,
                                                    rule=self._CNLSDDFG1__regression_rule(),
                                                    doc='regression equation')
        self.__model__.translation_rule = Constraint(self.__model__.I,
                                                     rule=self._CNLSDDFG1__translation_property(),
                                                     doc='translation property')

        self.__model__.afriat_rule = Constraint(self.__model__.I,
                                                rule=self._CNLSDDFG1__afriat_rule(),
                                                doc='elementary Afriat approach')

        self.__model__.sweet_rule = Constraint(self.__model__.I,
                                               self.__model__.I,
                                               rule=self._CNLSDDFG1__sweet_rule(),
                                               doc='sweet spot approach')
        self.__model__.sweet_rule2 = Constraint(self.__model__.I,
                                               self.__model__.I,
                                               rule=self._CNLSDDFZG2__sweet_rule2(),
                                               doc='sweet spot-2 approach')

        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0


