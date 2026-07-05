# -*- coding: utf-8 -*-
"""Rebuilt weak CQER/CER and weak CQERDDF/CERDDF models."""

from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, Constraint
from pyomo.core.expr.numvalue import NumericValue
import numpy as np
import pandas as pd
from .constant import CET_ADDI, CET_MULT, FUN_PROD, FUN_COST, RTS_CRS, RTS_VRS, RTS_VRS1, OPT_LOCAL, OPT_DEFAULT
from .utils import tools


def _as_2d(value):
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    return arr


def _parse_data_sent(data, sent, z=None, gy=[1], gx=[1], gb=[1], baseindex=None, refindex=None):
    outputvars, inputvars, unoutputvars, zvars, gy, gx, gb = tools.assert_valid_yxbz(sent, gy, gx, gb, z)
    y, x, b, zdata, yref, xref, bref, zref, referenceflag = tools.assert_valid_yxbz2(
        baseindex, refindex, data, outputvars, inputvars, unoutputvars, zvars)
    return _as_2d(y), _as_2d(x), _as_2d(b), (None if type(zdata) == type(None) else _as_2d(zdata)), gy, gx, gb


class CQRweak:
    """Convex quantile regression with weak disposability."""

    def __init__(self, data, sent, tau, z=None, cet=CET_ADDI, fun=FUN_PROD, rts=RTS_VRS, baseindex=None, refindex=None):
        self.y, self.x, self.b, self.z, self.gy, self.gx, self.gb = _parse_data_sent(
            data, sent, z=z, gy=[1], gx=[1] * (len(sent.split('=')[0].split())), gb=[1],
            baseindex=baseindex, refindex=refindex)
        self.tau = tau
        self.cet = cet
        self.fun = fun
        self.rts = rts
        self.__model__ = ConcreteModel()
        if type(self.z) != type(None):
            self.__model__.M = Set(initialize=range(len(self.z[0])))
            self.__model__.lamda = Var(self.__model__.M, doc='z coefficient')
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))
        self.__model__.L = Set(initialize=range(len(self.b[0])))
        self.__model__.alpha = Var(self.__model__.I, doc='alpha')
        self.__model__.beta = Var(self.__model__.I, self.__model__.J, bounds=(0.0, None), doc='beta')
        self.__model__.delta = Var(self.__model__.I, self.__model__.L, bounds=(0.0, None), doc='delta')
        self.__model__.epsilon_plus = Var(self.__model__.I, bounds=(0.0, None), doc='positive error term')
        self.__model__.epsilon_minus = Var(self.__model__.I, bounds=(0.0, None), doc='negative error term')
        self.__model__.frontier = Var(self.__model__.I, bounds=(0.0, None), doc='estimated frontier')
        self.__model__.objective = Objective(rule=self.__objective_rule(), sense=minimize, doc='objective function')
        self.__model__.regression_rule = Constraint(self.__model__.I, rule=self.__regression_rule(), doc='regression equation')
        if self.cet == CET_MULT:
            raise ValueError('CQRweak currently supports CET_ADDI only.')
        self.__model__.afriat_rule = Constraint(self.__model__.I, self.__model__.I, rule=self.__afriat_rule(), doc='afriat inequality')
        self.__model__.disposability_rule = Constraint(self.__model__.I, self.__model__.I, rule=self.__disposability_rule(), doc='weak disposability')
        self.optimization_status = 0
        self.problem_status = 0

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        self.problem_status, self.optimization_status = tools.optimize_model(self.__model__, email, self.cet, solver)

    def __objective_rule(self):
        def objective_rule(model):
            return sum(self.tau * model.epsilon_plus[i] + (1 - self.tau) * model.epsilon_minus[i] for i in model.I)
        return objective_rule

    def __regression_rule(self):
        if self.rts == RTS_VRS or self.rts == RTS_VRS1:
            def regression_rule(model, i):
                expr = model.alpha[i] + sum(model.beta[i, j] * self.x[i][j] for j in model.J) + sum(model.delta[i, l] * self.b[i][l] for l in model.L)
                if type(self.z) != type(None):
                    expr = expr - sum(model.lamda[m] * self.z[i][m] for m in model.M)
                return self.y[i][0] == expr - model.epsilon_minus[i] + model.epsilon_plus[i]
            return regression_rule
        if self.rts == RTS_CRS:
            def regression_rule(model, i):
                expr = sum(model.beta[i, j] * self.x[i][j] for j in model.J) + sum(model.delta[i, l] * self.b[i][l] for l in model.L)
                if type(self.z) != type(None):
                    expr = expr - sum(model.lamda[m] * self.z[i][m] for m in model.M)
                return self.y[i][0] == expr - model.epsilon_minus[i] + model.epsilon_plus[i]
            return regression_rule
        raise ValueError('Undefined rts.')

    def __afriat_rule(self):
        op = NumericValue.__le__ if self.fun == FUN_PROD else NumericValue.__ge__
        if self.rts == RTS_VRS or self.rts == RTS_VRS1:
            def afriat_rule(model, i, h):
                if i == h:
                    return Constraint.Skip
                return op(model.alpha[i] + sum(model.beta[i, j] * self.x[i][j] for j in model.J) + sum(model.delta[i, l] * self.b[i][l] for l in model.L),
                          model.alpha[h] + sum(model.beta[h, j] * self.x[i][j] for j in model.J) + sum(model.delta[h, l] * self.b[i][l] for l in model.L))
            return afriat_rule
        if self.rts == RTS_CRS:
            def afriat_rule(model, i, h):
                if i == h:
                    return Constraint.Skip
                return op(sum(model.beta[i, j] * self.x[i][j] for j in model.J) + sum(model.delta[i, l] * self.b[i][l] for l in model.L),
                          sum(model.beta[h, j] * self.x[i][j] for j in model.J) + sum(model.delta[h, l] * self.b[i][l] for l in model.L))
            return afriat_rule
        raise ValueError('Undefined rts.')

    def __disposability_rule(self):
        if self.rts == RTS_VRS or self.rts == RTS_VRS1:
            if self.fun == FUN_PROD:
                def disposability_rule(model, i, h):
                    if i == h:
                        return Constraint.Skip
                    return model.alpha[h] + sum(model.beta[h, j] * self.x[i][j] for j in model.J) >= 0
                return disposability_rule
            if self.fun == FUN_COST:
                def disposability_rule(model, i, h):
                    if i == h:
                        return Constraint.Skip
                    return model.alpha[h] + sum(model.beta[h, j] * self.x[i][j] for j in model.J) <= 0
                return disposability_rule
        if self.rts == RTS_CRS:
            if self.fun == FUN_PROD:
                def disposability_rule(model, i, h):
                    if i == h:
                        return Constraint.Skip
                    return sum(model.beta[h, j] * self.x[i][j] for j in model.J) >= 0
                return disposability_rule
            if self.fun == FUN_COST:
                def disposability_rule(model, i, h):
                    if i == h:
                        return Constraint.Skip
                    return sum(model.beta[h, j] * self.x[i][j] for j in model.J) <= 0
                return disposability_rule
        raise ValueError('Undefined model parameters.')

    def _single_model(self):
        return self.__model__

    def get_status(self):
        return self.optimization_status

    def get_alpha(self):
        tools.assert_optimized(self.optimization_status)
        return np.asarray(list(self.__model__.alpha[:].value), dtype=float)

    def get_beta(self):
        tools.assert_optimized(self.optimization_status)
        beta = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.beta), list(self.__model__.beta[:, :].value))])
        beta = pd.DataFrame(beta, columns=['Name', 'Key', 'Value']).pivot(index='Name', columns='Key', values='Value')
        return beta.to_numpy(dtype=float)

    def get_delta(self):
        tools.assert_optimized(self.optimization_status)
        delta = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.delta), list(self.__model__.delta[:, :].value))])
        delta = pd.DataFrame(delta, columns=['Name', 'Key', 'Value']).pivot(index='Name', columns='Key', values='Value')
        return delta.to_numpy(dtype=float)

    def get_residual(self):
        tools.assert_optimized(self.optimization_status)
        ep = np.asarray(list(self.__model__.epsilon_plus[:].value), dtype=float)
        em = np.asarray(list(self.__model__.epsilon_minus[:].value), dtype=float)
        return ep - em

    def get_frontier(self):
        tools.assert_optimized(self.optimization_status)
        return np.asarray(list(self.__model__.frontier[:].value), dtype=float)


class CERweak(CQRweak):
    """Convex expectile regression with weak disposability."""

    def _CQRweak__objective_rule(self):
        def objective_rule(model):
            return sum(self.tau * model.epsilon_plus[i] ** 2 + (1 - self.tau) * model.epsilon_minus[i] ** 2 for i in model.I)
        return objective_rule


class CQRDDFweak:
    """Weak convex quantile regression with directional distance function."""

    def __init__(self, data, sent, tau, z=None, gy=[1], gx=[1], gb=[1], fun=FUN_PROD, rts=RTS_VRS, baseindex=None, refindex=None):
        self.y, self.x, self.b, self.z, self.gy, self.gx, self.gb = _parse_data_sent(
            data, sent, z=z, gy=gy, gx=gx, gb=gb, baseindex=baseindex, refindex=refindex)
        self.tau = tau
        self.fun = fun
        self.rts = rts
        self.__model__ = ConcreteModel()
        if type(self.z) != type(None):
            self.__model__.M = Set(initialize=range(len(self.z[0])))
            self.__model__.lamda = Var(self.__model__.M, doc='z coefficient')
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))
        self.__model__.K = Set(initialize=range(len(self.y[0])))
        self.__model__.L = Set(initialize=range(len(self.b[0])))
        self.__model__.alpha = Var(self.__model__.I, doc='alpha')
        self.__model__.beta = Var(self.__model__.I, self.__model__.J, bounds=(0.0, None), doc='beta')
        self.__model__.gamma = Var(self.__model__.I, self.__model__.K, bounds=(0.0, None), doc='gamma')
        self.__model__.delta = Var(self.__model__.I, self.__model__.L, bounds=(0.0, None), doc='delta')
        self.__model__.epsilon_plus = Var(self.__model__.I, bounds=(0.0, None), doc='positive error term')
        self.__model__.epsilon_minus = Var(self.__model__.I, bounds=(0.0, None), doc='negative error term')
        self.__model__.objective = Objective(rule=self.__objective_rule(), sense=minimize, doc='objective function')
        self.__model__.regression_rule = Constraint(self.__model__.I, rule=self.__regression_rule(), doc='regression equation')
        self.__model__.translation_rule = Constraint(self.__model__.I, rule=self.__translation_property(), doc='translation property')
        self.__model__.afriat_rule = Constraint(self.__model__.I, self.__model__.I, rule=self.__afriat_rule(), doc='afriat inequality')
        self.__model__.disposability_rule = Constraint(self.__model__.I, self.__model__.I, rule=self.__disposability_rule(), doc='weak disposability')
        self.optimization_status = 0
        self.problem_status = 0

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        self.problem_status, self.optimization_status = tools.optimize_model(self.__model__, email, CET_ADDI, solver)

    def __objective_rule(self):
        def objective_rule(model):
            return sum(self.tau * model.epsilon_plus[i] + (1 - self.tau) * model.epsilon_minus[i] for i in model.I)
        return objective_rule

    def __regression_rule(self):
        if self.rts == RTS_VRS or self.rts == RTS_VRS1:
            def regression_rule(model, i):
                expr = model.alpha[i] + sum(model.beta[i, j] * self.x[i][j] for j in model.J) + sum(model.delta[i, l] * self.b[i][l] for l in model.L)
                if type(self.z) != type(None):
                    expr = expr - sum(model.lamda[m] * self.z[i][m] for m in model.M)
                return sum(model.gamma[i, k] * self.y[i][k] for k in model.K) == expr - model.epsilon_minus[i] + model.epsilon_plus[i]
            return regression_rule
        if self.rts == RTS_CRS:
            def regression_rule(model, i):
                expr = sum(model.beta[i, j] * self.x[i][j] for j in model.J) + sum(model.delta[i, l] * self.b[i][l] for l in model.L)
                if type(self.z) != type(None):
                    expr = expr - sum(model.lamda[m] * self.z[i][m] for m in model.M)
                return sum(model.gamma[i, k] * self.y[i][k] for k in model.K) == expr - model.epsilon_minus[i] + model.epsilon_plus[i]
            return regression_rule
        raise ValueError('Undefined rts.')

    def __translation_property(self):
        def translation_rule(model, i):
            return sum(model.beta[i, j] * self.gx[j] for j in model.J) + sum(model.gamma[i, k] * self.gy[k] for k in model.K) + sum(model.delta[i, l] * self.gb[l] for l in model.L) == 1
        return translation_rule

    def __afriat_rule(self):
        op = NumericValue.__le__ if self.fun == FUN_PROD else NumericValue.__ge__
        if self.rts == RTS_VRS or self.rts == RTS_VRS1:
            def afriat_rule(model, i, h):
                if i == h:
                    return Constraint.Skip
                return op(model.alpha[i] + sum(model.beta[i, j] * self.x[i][j] for j in model.J) - sum(model.gamma[i, k] * self.y[i][k] for k in model.K) + sum(model.delta[i, l] * self.b[i][l] for l in model.L),
                          model.alpha[h] + sum(model.beta[h, j] * self.x[i][j] for j in model.J) - sum(model.gamma[h, k] * self.y[i][k] for k in model.K) + sum(model.delta[h, l] * self.b[i][l] for l in model.L))
            return afriat_rule
        if self.rts == RTS_CRS:
            def afriat_rule(model, i, h):
                if i == h:
                    return Constraint.Skip
                return op(sum(model.beta[i, j] * self.x[i][j] for j in model.J) - sum(model.gamma[i, k] * self.y[i][k] for k in model.K) + sum(model.delta[i, l] * self.b[i][l] for l in model.L),
                          sum(model.beta[h, j] * self.x[i][j] for j in model.J) - sum(model.gamma[h, k] * self.y[i][k] for k in model.K) + sum(model.delta[h, l] * self.b[i][l] for l in model.L))
            return afriat_rule
        raise ValueError('Undefined rts.')

    def __disposability_rule(self):
        if self.rts == RTS_VRS or self.rts == RTS_VRS1:
            if self.fun == FUN_PROD:
                def disposability_rule(model, i, h):
                    if i == h:
                        return Constraint.Skip
                    return model.alpha[h] + sum(model.beta[h, j] * self.x[i][j] for j in model.J) >= 0
                return disposability_rule
            if self.fun == FUN_COST:
                def disposability_rule(model, i, h):
                    if i == h:
                        return Constraint.Skip
                    return model.alpha[h] + sum(model.beta[h, j] * self.x[i][j] for j in model.J) <= 0
                return disposability_rule
        if self.rts == RTS_CRS:
            if self.fun == FUN_PROD:
                def disposability_rule(model, i, h):
                    if i == h:
                        return Constraint.Skip
                    return sum(model.beta[h, j] * self.x[i][j] for j in model.J) >= 0
                return disposability_rule
            if self.fun == FUN_COST:
                def disposability_rule(model, i, h):
                    if i == h:
                        return Constraint.Skip
                    return sum(model.beta[h, j] * self.x[i][j] for j in model.J) <= 0
                return disposability_rule
        raise ValueError('Undefined model parameters.')

    def _single_model(self):
        return self.__model__

    def get_status(self):
        return self.optimization_status

    def get_alpha(self):
        tools.assert_optimized(self.optimization_status)
        return np.asarray(list(self.__model__.alpha[:].value), dtype=float)

    def get_beta(self):
        tools.assert_optimized(self.optimization_status)
        beta = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.beta), list(self.__model__.beta[:, :].value))])
        beta = pd.DataFrame(beta, columns=['Name', 'Key', 'Value']).pivot(index='Name', columns='Key', values='Value')
        return beta.to_numpy(dtype=float)

    def get_gamma(self):
        tools.assert_optimized(self.optimization_status)
        gamma = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.gamma), list(self.__model__.gamma[:, :].value))])
        gamma = pd.DataFrame(gamma, columns=['Name', 'Key', 'Value']).pivot(index='Name', columns='Key', values='Value')
        return gamma.to_numpy(dtype=float)

    def get_delta(self):
        tools.assert_optimized(self.optimization_status)
        delta = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.delta), list(self.__model__.delta[:, :].value))])
        delta = pd.DataFrame(delta, columns=['Name', 'Key', 'Value']).pivot(index='Name', columns='Key', values='Value')
        return delta.to_numpy(dtype=float)

    def get_residual(self):
        tools.assert_optimized(self.optimization_status)
        ep = np.asarray(list(self.__model__.epsilon_plus[:].value), dtype=float)
        em = np.asarray(list(self.__model__.epsilon_minus[:].value), dtype=float)
        return ep - em

    def get_frontier(self):
        tools.assert_optimized(self.optimization_status)
        return np.asarray([np.sum(self.get_gamma()[i, :] * self.y[i, :]) for i in range(len(self.y))], dtype=float)


class CERDDFweak(CQRDDFweak):
    """Weak convex expectile regression with directional distance function."""

    def _CQRDDFweak__objective_rule(self):
        def objective_rule(model):
            return sum(self.tau * model.epsilon_plus[i] ** 2 + (1 - self.tau) * model.epsilon_minus[i] ** 2 for i in model.I)
        return objective_rule

from .utils import CQRweakG1, CQRweakG2, CQRweakZG1, CQRweakZG2, sweet, interpolation
from .utils import CQRDDFweakG1, CQRDDFweakG2, CQRDDFweakZG1, CQRDDFweakZG2
import time

class CQRweakG:
    """Convex quantile regression (CQR) with Genetic algorithm and weak
    """

    def __init__(self, y, x, b, tau, z=None, cet=CET_ADDI, fun=FUN_PROD, rts=RTS_VRS):
        """CQRG model

        Args:
            y (ndarray): output variable.
            x (ndarray): input variables.
            b (ndarray): undersiable variables.
            tau (float): quantile.
            z (float, optional): Contextual variable(s). Defaults to None.
            cet (String, optional): CET_ADDI (additive composite error term) or CET_MULT (multiplicative composite error term). Defaults to CET_ADDI.
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS.
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.cutactive = sweet.sweet(np.column_stack((x,b)))
        self.y, self.x, self.b, self.z = tools.assert_valid_wp_data(y, x, b, z)

        self.tau = tau
        self.cet = cet
        self.fun = fun
        self.rts = rts

        # active (added) violated concavity constraint by iterative procedure
        self.active = np.zeros((len(x), len(x)))
        # violated concavity constraint
        self.active2 = np.zeros((len(x), len(x)))

        # active (added) violated concavity constraint for weak disposbility constrains by iterative procedure
        self.activeweak = np.zeros((len(x), len(x)))
        # violated concavity constraint for weak disposbility constrains
        self.activeweak2 = np.zeros((len(x), len(x)))

        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Optimize the function by requested method"""
        # TODO(error/warning handling): Check problem status after optimization
        self.t0 = time.time()
        if type(self.z) != type(None):
            model1 = CQRweakZG1.CQRweakZG1(
                self.y, self.x, self.b, self.z, self.tau, self.cutactive, self.cet, self.fun, self.rts)
        else:
            model1 = CQRweakG1.CQRweakG1(
                self.y, self.x, self.b, self.tau, self.cutactive, self.cet, self.fun, self.rts)
        model1.optimize(email, solver)
        self.alpha = model1.get_alpha()
        self.beta = model1.get_beta()
        self.delta = model1.get_delta()
        self.__model__ = model1.__model__

        self.count = 0
        while self.__convergence_test(self.alpha, self.beta, self.delta) > 0.0001 \
                or self.__convergence_test_weak(self.alpha, self.beta, self.delta) > 0.0001:
            if type(self.z) != type(None):
                model2 = CQRweakZG2.CQRweakZG2(
                    self.y, self.x, self.b, self.z, self.tau, self.cutactive,self.active,self.activeweak,
                                self.cet, self.fun, self.rts)
            else:
                model2 = CQRweakG2.CQRweakG2(
                    self.y, self.x, self.b, self.tau, self.cutactive, self.active, self.activeweak,
                                self.cet, self.fun, self.rts)
            model2.optimize(email, solver)
            self.alpha = model2.get_alpha()
            self.beta = model2.get_beta()
            self.delta = model2.get_delta()

            # TODO: Replace print with log system
            # print("Genetic Algorithm Convergence : %8f" %
            #       (self.__convergence_test(self.alpha, self.beta)))
            self.__model__ = model2.__model__
            self.count += 1
        self.optimization_status = 1
        self.tt = time.time() - self.t0

    def __convergence_test(self, alpha, beta,delta):
        alpha = np.asarray(alpha, dtype=float).reshape(-1)
        beta = np.asarray(beta, dtype=float)
        delta = np.asarray(delta, dtype=float)
        x = np.asarray(self.x)
        b = np.asarray(self.b)
        activetmp1 = 0.0
        if self.rts == RTS_VRS and self.fun == FUN_PROD:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = alpha[i] \
                        + np.sum(beta[i, :] * x[i, :]) + np.sum(delta[i, :] * b[i, :]) - \
                        alpha[j] - np.sum(beta[j, :] * x[i, :]) - np.sum(delta[j, :] * b[i, :])

                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp


        elif self.rts == RTS_VRS and self.fun == FUN_COST:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = - alpha[i] - np.sum(beta[i, :] * x[i, :]) \
                                         - np.sum(delta[i, :] * b[i, :])+ alpha[j] \
                                         + np.sum(beta[j, :] * x[i, :]) + np.sum(delta[j, :] * b[i, :])

                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp


        elif self.rts == RTS_CRS and self.fun == FUN_PROD:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = np.sum(beta[i, :] * x[i, :]) + np.sum(delta[i, :] * b[i, :]) \
                            - np.sum(beta[j, :] * x[i, :]) - np.sum(delta[j, :] * b[i, :])

                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_CRS and self.fun == FUN_COST:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = - np.sum(beta[i, :] * x[i, :]) - np.sum(delta[i, :] * b[i, :]) \
                            + np.sum(beta[j, :] * x[i, :]) + np.sum(delta[j, :] * b[i, :])
                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp
        return activetmp


    def __convergence_test_weak(self, alpha, beta,delta):
        alpha = np.asarray(alpha, dtype=float).reshape(-1)
        beta = np.asarray(beta, dtype=float)
        delta = np.asarray(delta, dtype=float)
        x = np.asarray(self.x)
        b = np.asarray(self.b)
        activetmp1 = 0.0
        if self.rts == RTS_VRS and self.fun == FUN_PROD:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.activeweak2[i, j] = - alpha[j] - np.sum(beta[j, :] * x[i, :])

                    if self.activeweak2[i, j] > activetmp:
                        activetmp = self.activeweak2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.activeweak2[i, j] >= activetmp and activetmp > 0:
                        self.activeweak[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_VRS and self.fun == FUN_COST:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.activeweak2[i, j] = alpha[j] + np.sum(beta[j, :] * x[i, :])

                    if self.activeweak2[i, j] > activetmp:
                        activetmp = self.activeweak2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.activeweak2[i, j] >= activetmp and activetmp > 0:
                        self.activeweak[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_CRS and self.fun == FUN_PROD:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.activeweak2[i, j] = - np.sum(beta[j, :] * x[i, :])

                    if self.activeweak2[i, j] > activetmp:
                        activetmp = self.activeweak2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.activeweak2[i, j] >= activetmp and activetmp > 0:
                        self.activeweak[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_CRS and self.fun == FUN_COST:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.activeweak2[i, j] = np.sum(beta[j, :] * x[i, :])
                    if self.activeweak2[i, j] > activetmp:
                        activetmp = self.activeweak2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.activeweak2[i, j] >= activetmp and activetmp > 0:
                        self.activeweak[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp
        return activetmp


    def display_status(self):
        """Display the status of problem"""
        print(self.optimization_status)

    def display_alpha(self):
        """Display alpha value"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_various_return_to_scale(self.rts)
        self.__model__.alpha.display()

    def display_beta(self):
        """Display beta value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.beta.display()

    def display_delta(self):
        """Display beta value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.delta.display()

    def display_lamda(self):
        """Display lamda value"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_contextual_variable(self.z)
        self.__model__.lamda.display()

    # def display_residual(self):
    #     """Dispaly residual value"""
    #     tools.assert_optimized(self.optimization_status)
    #     self.__model__.epsilon.display()

    def display_positive_residual(self):
        """Dispaly positive residual value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.epsilon_plus.display()

    def display_negative_residual(self):
        """Dispaly negative residual value"""
        tools.assert_optimized(self.optimization_status)
        self.__model__.epsilon_minus.display()

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

    def get_delta(self):
        """Return delta value by array"""
        if self.optimization_status == 0:
            self.optimize()
        delta = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.delta),
                                                           list(self.__model__.delta[:, :].value))])
        delta = pd.DataFrame(delta, columns=['Name', 'Key', 'Value'])
        delta = delta.pivot(index='Name', columns='Key', values='Value')
        return delta.to_numpy()

    def get_residual(self):
        """Return residual value by array"""
        tools.assert_optimized(self.optimization_status)
        residual = np.asarray(list(self.__model__.epsilon_plus[:].value)) - np.asarray(list(self.__model__.epsilon_minus[:].value))
        return residual

    def get_positive_residual(self):
        """Return positive residual value by array"""
        tools.assert_optimized(self.optimization_status)
        residual_plus = list(self.__model__.epsilon_plus[:].value)
        return np.asarray(residual_plus)

    def get_negative_residual(self):
        """Return negative residual value by array"""
        tools.assert_optimized(self.optimization_status)
        residual_minus = list(self.__model__.epsilon_minus[:].value)
        return np.asarray(residual_minus)

    def get_lamda(self):
        """Return beta value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_contextual_variable(self.z)
        lamda = list(self.__model__.lamda[:].value)
        return np.asarray(lamda)

    def get_frontier(self):
        """Return estimated frontier value by array"""
        tools.assert_optimized(self.optimization_status)
        if self.cet == CET_MULT and type(self.z) == type(None):
            frontier = np.asarray(list(self.__model__.frontier[:].value))+1
        elif self.cet == CET_MULT and type(self.z) != type(None):
            frontier = list(np.divide(self.y, np.exp(
                - self.get_residual() - self.get_lamda()*np.asarray(self.z)[:, 0])) - 1)
        elif self.cet == CET_ADDI:
            frontier = np.asarray(self.y) + self.get_residual()
        return np.asarray(frontier)

    def get_totalconstr(self):
        """Return the number of total constraints"""
        tools.assert_optimized(self.optimization_status)
        activeconstr = np.sum(self.active) - np.trace(self.active)
        cutactiveconstr = np.sum(self.cutactive) - np.trace(self.cutactive)
        totalconstr = activeconstr + cutactiveconstr + 2 * len(self.active) + 1
        return totalconstr

    def get_runningtime(self):
        """Return the running time"""
        tools.assert_optimized(self.optimization_status)
        return self.tt

    def get_blocks(self):
        """Return the number of blocks"""
        tools.assert_optimized(self.optimization_status)
        return self.count

    def get_predict(self, x_test):
        """Return the estimated function in testing sample"""
        tools.assert_optimized(self.optimization_status)
        return interpolation.interpolation(self.get_alpha(), self.get_beta(), x_test, fun=self.fun)
    

class CERweakG(CQRweakG):
    """Convex expectile regression (CER) with Genetic algorithm
    """

    def __init__(self, y, x, b, tau, z=None, cet=CET_ADDI, fun=FUN_PROD, rts=RTS_VRS):
        """CERG model

        Args:
            y (float): output variable. 
            x (float): input variables.
            b (float): undersiable variables.
            tau (float): quantile.
            z (float, optional): Contextual variable(s). Defaults to None.
            cet (String, optional): CET_ADDI (additive composite error term) or CET_MULT (multiplicative composite error term). Defaults to CET_ADDI.
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS.
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.cutactive = sweet.sweet(x)
        self.y, self.x, self.b, self.z = tools.assert_valid_wp_data(y, x, b, z)

        self.tau = tau
        self.cet = cet
        self.fun = fun
        self.rts = rts

        # active (added) violated concavity constraint by iterative procedure
        self.active = np.zeros((len(x), len(x)))
        # violated concavity constraint
        self.active2 = np.zeros((len(x), len(x)))

        # active (added) violated concavity constraint for weak disposbility constrains by iterative procedure
        self.activeweak = np.zeros((len(x), len(x)))
        # violated concavity constraint for weak disposbility constrains
        self.activeweak2 = np.zeros((len(x), len(x)))

        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0

    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Optimize the function by requested method"""
        # TODO(error/warning handling): Check problem status after optimization
        self.t0 = time.time()
        if type(self.z) != type(None):
            model1 = CQRweakZG1.CERweakZG1(
                self.y, self.x, self.b, self.z, self.tau, self.cutactive, self.cet, self.fun, self.rts)
        else:
            model1 = CQRweakG1.CERweakG1(
                self.y, self.x, self.b, self.tau, self.cutactive, self.cet, self.fun, self.rts)
        model1.optimize(email, solver)
        self.alpha = model1.get_alpha()
        self.beta = model1.get_beta()
        self.delta = model1.get_delta()
        self.__model__ = model1.__model__

        self.count = 0
        while self.__convergence_test(self.alpha, self.beta, self.delta) > 0.0001 \
                or self.__convergence_test_weak(self.alpha, self.beta, self.delta) > 0.0001:
            if type(self.z) != type(None):
                model2 = CQRweakZG2.CERweakZG2(
                    self.y, self.x, self.b, self.z, self.tau, self.cutactive,self.active,self.activeweak,
                                self.cet, self.fun, self.rts)
            else:
                model2 = CQRweakG2.CERweakG2(
                    self.y, self.x, self.b, self.tau, self.cutactive, self.active, self.activeweak,
                                self.cet, self.fun, self.rts)
            model2.optimize(email, solver)
            self.alpha = model2.get_alpha()
            self.beta = model2.get_beta()
            self.delta = model2.get_delta()

            # TODO: Replace print with log system
            # print("Genetic Algorithm Convergence : %8f" %
            #       (self.__convergence_test(self.alpha, self.beta)))
            self.__model__ = model2.__model__
            self.count += 1
        self.optimization_status = 1
        self.tt = time.time() - self.t0

    def __convergence_test(self, alpha, beta,delta):
        alpha = np.asarray(alpha, dtype=float).reshape(-1)
        beta = np.asarray(beta, dtype=float)
        delta = np.asarray(delta, dtype=float)
        x = np.asarray(self.x)
        b = np.asarray(self.b)
        activetmp1 = 0.0
        if self.rts == RTS_VRS and self.fun == FUN_PROD:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = alpha[i] \
                        + np.sum(beta[i, :] * x[i, :]) + np.sum(delta[i, :] * b[i, :]) - \
                        alpha[j] - np.sum(beta[j, :] * x[i, :]) - np.sum(delta[j, :] * b[i, :])

                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp


        elif self.rts == RTS_VRS and self.fun == FUN_COST:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = - alpha[i] - np.sum(beta[i, :] * x[i, :]) \
                                         - np.sum(delta[i, :] * b[i, :])+ alpha[j] \
                                         + np.sum(beta[j, :] * x[i, :]) + np.sum(delta[j, :] * b[i, :])

                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp


        elif self.rts == RTS_CRS and self.fun == FUN_PROD:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = np.sum(beta[i, :] * x[i, :]) + np.sum(delta[i, :] * b[i, :]) \
                            - np.sum(beta[j, :] * x[i, :]) - np.sum(delta[j, :] * b[i, :])

                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_CRS and self.fun == FUN_COST:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = - np.sum(beta[i, :] * x[i, :]) - np.sum(delta[i, :] * b[i, :]) \
                            + np.sum(beta[j, :] * x[i, :]) + np.sum(delta[j, :] * b[i, :])
                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp
        return activetmp


    def __convergence_test_weak(self, alpha, beta,delta):
        alpha = np.asarray(alpha, dtype=float).reshape(-1)
        beta = np.asarray(beta, dtype=float)
        delta = np.asarray(delta, dtype=float)
        x = np.asarray(self.x)
        b = np.asarray(self.b)
        activetmp1 = 0.0
        if self.rts == RTS_VRS and self.fun == FUN_PROD:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.activeweak2[i, j] = - alpha[j] - np.sum(beta[j, :] * x[i, :])

                    if self.activeweak2[i, j] > activetmp:
                        activetmp = self.activeweak2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.activeweak2[i, j] >= activetmp and activetmp > 0:
                        self.activeweak[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_VRS and self.fun == FUN_COST:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.activeweak2[i, j] = alpha[j] + np.sum(beta[j, :] * x[i, :])

                    if self.activeweak2[i, j] > activetmp:
                        activetmp = self.activeweak2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.activeweak2[i, j] >= activetmp and activetmp > 0:
                        self.activeweak[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_CRS and self.fun == FUN_PROD:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.activeweak2[i, j] = - np.sum(beta[j, :] * x[i, :])

                    if self.activeweak2[i, j] > activetmp:
                        activetmp = self.activeweak2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.activeweak2[i, j] >= activetmp and activetmp > 0:
                        self.activeweak[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_CRS and self.fun == FUN_COST:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.activeweak2[i, j] = np.sum(beta[j, :] * x[i, :])
                    if self.activeweak2[i, j] > activetmp:
                        activetmp = self.activeweak2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.activeweak2[i, j] >= activetmp and activetmp > 0:
                        self.activeweak[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp
        return activetmp


# ---------------------------------------------------------------------------
# Weak CQERDDF-G / CERDDF-G models
# ---------------------------------------------------------------------------
class CQRDDFweakG:
    """Convex Nonparametric Least Square with weak disposability (weakCNLSDDF) and Genetic algorithm
    """
    def __init__(self, y, x, b, tau, z=None, gy=[1], gx=[1], gb=[1], fun=FUN_PROD, rts=RTS_VRS):
        """weakCNLSG model

        Args:
            y (ndarray): output variable.
            x (ndarray): input variables.
            b (ndarray): undersiable variables.
            tau (float): quantile.
            z (ndarray, optional): Contextual variable(s). Defaults to None.
            gy (list, optional): output directional vector. Defaults to [1].
            gx (list, optional): input directional vector. Defaults to [1].
            gb (list, optional): undesirable output directional vector. Defaults to [1].
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS.
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.cutactive = sweet.sweet(np.column_stack((x,b)))
        self.y, self.x, self.b, self.z,self.gy, self.gx, self.gb = \
            tools.assert_valid_direciontal_data_with_z(y,x,b,z,gy,gx,gb)

        self.tau = tau

        self.fun = fun
        self.rts = rts

        # active (added) violated concavity constraint by iterative procedure
        self.active = np.zeros((len(x), len(x)))
        # violated concavity constraint
        self.active2 = np.zeros((len(x), len(x)))

        # active (added) violated concavity constraint for weak disposbility constrains by iterative procedure
        self.activeweak = np.zeros((len(x), len(x)))
        # violated concavity constraint for weak disposbility constrains
        self.activeweak2 = np.zeros((len(x), len(x)))

        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0



    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Optimize the function by requested method"""
        # TODO(error/warning handling): Check problem status after optimization
        self.t0 = time.time()
        if type(self.z) != type(None):
            model1 = CQRDDFweakZG1.CQRDDFweakZG1(
                self.y, self.x, self.b, self.z, self.tau,self.cutactive, self.gy, self.gx, self.gb, self.fun, self.rts)
        else:
            model1 = CQRDDFweakG1.CQRDDFweakG1(
                self.y, self.x, self.b,self.tau, self.cutactive, self.gy, self.gx, self.gb, self.fun, self.rts)
        model1.optimize(email, solver)
        self.alpha = model1.get_alpha()
        self.beta = model1.get_beta()
        self.gamma = model1.get_gamma()
        self.delta = model1.get_delta()
        self.__model__ = model1.__model__

        self.count = 0
        while self.__convergence_test(self.alpha, self.beta, self.gamma, self.delta) > 0.0001 \
                or self.__convergence_test_weak(self.alpha, self.beta, self.gamma, self.delta) > 0.0001:
            if type(self.z) != type(None):
                model2 = CQRDDFweakZG2.CQRDDFweakZG2(
                    self.y, self.x, self.b, self.z, self.tau, self.cutactive, self.active,self.activeweak,
                    self.gy, self.gx, self.gb,  self.fun, self.rts)
            else:
                model2 = CQRDDFweakG2.CQRDDFweakG2(
                    self.y, self.x, self.b, self.tau, self.cutactive, self.active,self.activeweak,
                    self.gy, self.gx, self.gb, self.fun, self.rts)
            model2.optimize(email, solver)
            self.alpha = model2.get_alpha()
            self.beta = model2.get_beta()
            self.gamma = model2.get_gamma()
            self.delta = model2.get_delta()
            # TODO: Replace print with log system
            # print("Genetic Algorithm Convergence : %8f" %
            #       (self.__convergence_test(self.alpha, self.beta)))
            self.__model__ = model2.__model__
            self.count += 1
        self.optimization_status = 1
        self.tt = time.time() - self.t0

    def __convergence_test(self, alpha, beta, gamma, delta):
        alpha = np.asarray(alpha, dtype=float).reshape(-1)
        beta = np.asarray(beta, dtype=float)
        gamma = np.asarray(gamma, dtype=float)
        delta = np.asarray(delta, dtype=float)
        x = np.asarray(self.x)
        y = np.asarray(self.y)
        b = np.asarray(self.b)
        activetmp1 = 0.0
        if self.rts == RTS_VRS and self.fun == FUN_PROD:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = alpha[i] + np.sum(beta[i, :] * x[i, :]) \
                        + np.sum(delta[i, :] * b[i, :]) - np.sum(gamma[i, :] * y[i, :]) \
                                         - alpha[j] - np.sum(beta[j, :] * x[i, :]) \
                        - np.sum(delta[j, :] * b[i, :]) + np.sum(gamma[j, :] * y[i, :])

                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp


        elif self.rts == RTS_VRS and self.fun == FUN_COST:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = -alpha[i] - np.sum(beta[i, :] * x[i, :]) \
                        - np.sum(delta[i, :] * b[i, :]) + np.sum(gamma[i, :] * y[i, :]) \
                                         + alpha[j] + np.sum(beta[j, :] * x[i, :]) \
                        + np.sum(delta[j, :] * b[i, :]) - np.sum(gamma[j, :] * y[i, :])

                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp


        elif self.rts == RTS_CRS and self.fun == FUN_PROD:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = np.sum(beta[i, :] * x[i, :]) \
                        + np.sum(delta[i, :] * b[i, :]) - np.sum(gamma[i, :] * y[i, :]) \
                                          - np.sum(beta[j, :] * x[i, :]) \
                        - np.sum(delta[j, :] * b[i, :]) + np.sum(gamma[j, :] * y[i, :])

                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_CRS and self.fun == FUN_COST:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = - np.sum(beta[i, :] * x[i, :]) \
                        - np.sum(delta[i, :] * b[i, :]) + np.sum(gamma[i, :] * y[i, :]) \
                                          + np.sum(beta[j, :] * x[i, :]) \
                        + np.sum(delta[j, :] * b[i, :]) - np.sum(gamma[j, :] * y[i, :])
                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp
        return activetmp



    def __convergence_test_weak(self, alpha, beta, gamma,delta):
        alpha = np.asarray(alpha, dtype=float).reshape(-1)
        beta = np.asarray(beta, dtype=float)
        gamma = np.asarray(gamma, dtype=float)
        delta = np.asarray(delta, dtype=float)
        x = np.asarray(self.x)
        y = np.asarray(self.y)
        b = np.asarray(self.b)
        activetmp1 = 0.0
        if self.rts == RTS_VRS and self.fun == FUN_PROD:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.activeweak2[i, j] = - alpha[j] - np.sum(beta[j, :] * x[i, :])

                    if self.activeweak2[i, j] > activetmp:
                        activetmp = self.activeweak2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.activeweak2[i, j] >= activetmp and activetmp > 0:
                        self.activeweak[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_VRS and self.fun == FUN_COST:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.activeweak2[i, j] = alpha[j] + np.sum(beta[j, :] * x[i, :])

                    if self.activeweak2[i, j] > activetmp:
                        activetmp = self.activeweak2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.activeweak2[i, j] >= activetmp and activetmp > 0:
                        self.activeweak[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_CRS and self.fun == FUN_PROD:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.activeweak2[i, j] = - np.sum(beta[j, :] * x[i, :])

                    if self.activeweak2[i, j] > activetmp:
                        activetmp = self.activeweak2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.activeweak2[i, j] >= activetmp and activetmp > 0:
                        self.activeweak[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_CRS and self.fun == FUN_COST:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.activeweak2[i, j] = np.sum(beta[j, :] * x[i, :])
                    if self.activeweak2[i, j] > activetmp:
                        activetmp = self.activeweak2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.activeweak2[i, j] >= activetmp and activetmp > 0:
                        self.activeweak[i, j] = 1
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

    # def display_residual(self):
    #     """Dispaly residual value"""
    #     tools.assert_optimized(self.optimization_status)
    #     self.__model__.epsilon.display()

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
        residual = np.asarray(list(self.__model__.epsilon_plus[:].value)) - np.asarray(list(self.__model__.epsilon_minus[:].value))
        return residual

    def get_lamda(self):
        """Return beta value by array"""
        tools.assert_optimized(self.optimization_status)
        tools.assert_contextual_variable(self.z)
        lamda = list(self.__model__.lamda[:].value)
        return np.asarray(lamda)

    def get_frontier(self):
        """Return estimated frontier value by array"""
        raise ValueError("DDF hsa no frontier.")

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


class CERDDFweakG(CQRDDFweakG):
    """Convex Nonparametric Least Square with weak disposability (weakCNLSDDF) and Genetic algorithm
    """
    def __init__(self, y, x, b, tau, z=None, gy=[1], gx=[1], gb=[1], fun=FUN_PROD, rts=RTS_VRS):
        """weakCNLSG model

        Args:
            y (ndarray): output variable.
            x (ndarray): input variables.
            b (ndarray): undersiable variables.
            tau (float): quantile.
            z (ndarray, optional): Contextual variable(s). Defaults to None.
            gy (list, optional): output directional vector. Defaults to [1].
            gx (list, optional): input directional vector. Defaults to [1].
            gb (list, optional): undesirable output directional vector. Defaults to [1].
            fun (String, optional): FUN_PROD (production frontier) or FUN_COST (cost frontier). Defaults to FUN_PROD.
            rts (String, optional): RTS_VRS (variable returns to scale) or RTS_CRS (constant returns to scale). Defaults to RTS_VRS.
        """
        # TODO(error/warning handling): Check the configuration of the model exist
        self.cutactive = sweet.sweet(np.column_stack((x,b)))
        self.y, self.x, self.b, self.z,self.gy, self.gx, self.gb = \
            tools.assert_valid_direciontal_data_with_z(y,x,b,z,gy,gx,gb)

        self.tau = tau

        self.fun = fun
        self.rts = rts

        # active (added) violated concavity constraint by iterative procedure
        self.active = np.zeros((len(x), len(x)))
        # violated concavity constraint
        self.active2 = np.zeros((len(x), len(x)))

        # active (added) violated concavity constraint for weak disposbility constrains by iterative procedure
        self.activeweak = np.zeros((len(x), len(x)))
        # violated concavity constraint for weak disposbility constrains
        self.activeweak2 = np.zeros((len(x), len(x)))

        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0



    def optimize(self, email=OPT_LOCAL, solver=OPT_DEFAULT):
        """Optimize the function by requested method"""
        # TODO(error/warning handling): Check problem status after optimization
        self.t0 = time.time()
        if type(self.z) != type(None):
            model1 = CQRDDFweakZG1.CERDDFweakZG1(
                self.y, self.x, self.b, self.z,self.tau, self.cutactive, self.gy, self.gx, self.gb, self.fun, self.rts)
        else:
            model1 = CQRDDFweakG1.CERDDFweakG1(
                self.y, self.x, self.b, self.tau, self.cutactive, self.gy, self.gx, self.gb, self.fun, self.rts)
        model1.optimize(email, solver)
        self.alpha = model1.get_alpha()
        self.beta = model1.get_beta()
        self.gamma = model1.get_gamma()
        self.delta = model1.get_delta()
        self.__model__ = model1.__model__

        self.count = 0
        while self.__convergence_test(self.alpha, self.beta, self.gamma, self.delta) > 0.0001 \
                or self.__convergence_test_weak(self.alpha, self.beta, self.gamma, self.delta) > 0.0001:
            if type(self.z) != type(None):
                model2 = CQRDDFweakZG2.CERDDFweakZG2(
                    self.y, self.x, self.b, self.z, self.tau, self.cutactive, self.active,self.activeweak,
                    self.gy, self.gx, self.gb,  self.fun, self.rts)
            else:
                model2 = CQRDDFweakG2.CERDDFweakG2(
                    self.y, self.x, self.b, self.tau, self.cutactive, self.active,self.activeweak,
                    self.gy, self.gx, self.gb, self.fun, self.rts)
            model2.optimize(email, solver)
            self.alpha = model2.get_alpha()
            self.beta = model2.get_beta()
            self.gamma = model2.get_gamma()
            self.delta = model2.get_delta()
            # TODO: Replace print with log system
            # print("Genetic Algorithm Convergence : %8f" %
            #       (self.__convergence_test(self.alpha, self.beta)))
            self.__model__ = model2.__model__
            self.count += 1
        self.optimization_status = 1
        self.tt = time.time() - self.t0

    def __convergence_test(self, alpha, beta, gamma, delta):
        alpha = np.asarray(alpha, dtype=float).reshape(-1)
        beta = np.asarray(beta, dtype=float)
        gamma = np.asarray(gamma, dtype=float)
        delta = np.asarray(delta, dtype=float)
        x = np.asarray(self.x)
        y = np.asarray(self.y)
        b = np.asarray(self.b)
        activetmp1 = 0.0
        if self.rts == RTS_VRS and self.fun == FUN_PROD:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = alpha[i] + np.sum(beta[i, :] * x[i, :]) \
                        + np.sum(delta[i, :] * b[i, :]) - np.sum(gamma[i, :] * y[i, :]) \
                                         - alpha[j] - np.sum(beta[j, :] * x[i, :]) \
                        - np.sum(delta[j, :] * b[i, :]) + np.sum(gamma[j, :] * y[i, :])

                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp


        elif self.rts == RTS_VRS and self.fun == FUN_COST:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = -alpha[i] - np.sum(beta[i, :] * x[i, :]) \
                        - np.sum(delta[i, :] * b[i, :]) + np.sum(gamma[i, :] * y[i, :]) \
                                         + alpha[j] + np.sum(beta[j, :] * x[i, :]) \
                        + np.sum(delta[j, :] * b[i, :]) - np.sum(gamma[j, :] * y[i, :])

                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp


        elif self.rts == RTS_CRS and self.fun == FUN_PROD:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = np.sum(beta[i, :] * x[i, :]) \
                        + np.sum(delta[i, :] * b[i, :]) - np.sum(gamma[i, :] * y[i, :]) \
                                          - np.sum(beta[j, :] * x[i, :]) \
                        - np.sum(delta[j, :] * b[i, :]) + np.sum(gamma[j, :] * y[i, :])

                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_CRS and self.fun == FUN_COST:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.active2[i, j] = - np.sum(beta[i, :] * x[i, :]) \
                        - np.sum(delta[i, :] * b[i, :]) + np.sum(gamma[i, :] * y[i, :]) \
                                          + np.sum(beta[j, :] * x[i, :]) \
                        + np.sum(delta[j, :] * b[i, :]) - np.sum(gamma[j, :] * y[i, :])
                    if self.active2[i, j] > activetmp:
                        activetmp = self.active2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.active2[i, j] >= activetmp and activetmp > 0:
                        self.active[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp
        return activetmp



    def __convergence_test_weak(self, alpha, beta, gamma,delta):
        alpha = np.asarray(alpha, dtype=float).reshape(-1)
        beta = np.asarray(beta, dtype=float)
        gamma = np.asarray(gamma, dtype=float)
        delta = np.asarray(delta, dtype=float)
        x = np.asarray(self.x)
        y = np.asarray(self.y)
        b = np.asarray(self.b)
        activetmp1 = 0.0
        if self.rts == RTS_VRS and self.fun == FUN_PROD:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.activeweak2[i, j] = - alpha[j] - np.sum(beta[j, :] * x[i, :])

                    if self.activeweak2[i, j] > activetmp:
                        activetmp = self.activeweak2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.activeweak2[i, j] >= activetmp and activetmp > 0:
                        self.activeweak[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_VRS and self.fun == FUN_COST:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.activeweak2[i, j] = alpha[j] + np.sum(beta[j, :] * x[i, :])

                    if self.activeweak2[i, j] > activetmp:
                        activetmp = self.activeweak2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.activeweak2[i, j] >= activetmp and activetmp > 0:
                        self.activeweak[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_CRS and self.fun == FUN_PROD:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.activeweak2[i, j] = - np.sum(beta[j, :] * x[i, :])

                    if self.activeweak2[i, j] > activetmp:
                        activetmp = self.activeweak2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.activeweak2[i, j] >= activetmp and activetmp > 0:
                        self.activeweak[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp

        elif self.rts == RTS_CRS and self.fun == FUN_COST:
        # go into the loop
            for i in range(len(x)):
                activetmp = 0.0
                # go into the sub-loop and find the violated concavity constraints
                for j in range(len(x)):
                    self.activeweak2[i, j] = np.sum(beta[j, :] * x[i, :])
                    if self.activeweak2[i, j] > activetmp:
                        activetmp = self.activeweak2[i, j]

            # find the maximal violated constraint in sub-loop and added into the active matrix
                for j in range(len(x)):
                    if self.activeweak2[i, j] >= activetmp and activetmp > 0:
                        self.activeweak[i, j] = 1
                if activetmp > activetmp1:
                    activetmp1 = activetmp
        return activetmp

# ---------------------------------------------------------------------------
# Deterministic efficiency helpers for weak CQER/CER models
# ---------------------------------------------------------------------------
def _cqerw_as_2d(value):
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 0:
        return arr.reshape(1, 1)
    if arr.ndim == 1:
        return arr.reshape(-1, 1)
    return arr


def _cqerw_as_1d(value):
    arr = np.asarray(value, dtype=float)
    if arr.ndim == 2 and arr.shape[1] == 1:
        return arr[:, 0]
    return arr.reshape(-1)


def _cqerw_to_numpy(value):
    if hasattr(value, "to_numpy"):
        return np.asarray(value.to_numpy(), dtype=float)
    return np.asarray(value, dtype=float)


def _cqerw_alpha(self):
    n = len(_cqerw_as_2d(self.y))
    try:
        alpha = _cqerw_to_numpy(self.get_alpha()).reshape(-1)
        if alpha.size == n:
            return alpha
    except Exception:
        pass
    return np.zeros(n, dtype=float)


def _cqerw_lamda_effect(self):
    if getattr(self, "z", None) is None:
        return None
    try:
        z = _cqerw_as_2d(self.z)
        lamda = _cqerw_to_numpy(self.get_lamda()).reshape(-1)
        return z @ lamda
    except Exception:
        return None


def _cqerw_linear_frontier(self):
    tools.assert_optimized(self.optimization_status)
    x = _cqerw_as_2d(self.x)
    b = _cqerw_as_2d(self.b)
    beta = _cqerw_as_2d(_cqerw_to_numpy(self.get_beta()))
    delta = _cqerw_as_2d(_cqerw_to_numpy(self.get_delta()))
    value = _cqerw_alpha(self) + np.sum(beta * x, axis=1) + np.sum(delta * b, axis=1)
    z_effect = _cqerw_lamda_effect(self)
    if z_effect is not None:
        if getattr(self, "cet", CET_ADDI) == CET_MULT:
            value = value * np.exp(-z_effect)
        else:
            value = value - z_effect
    return np.asarray(value, dtype=float)


def _cqerw_get_frontier(self):
    """Return the fitted deterministic weak CQR/CER frontier at observed DMUs."""
    return _cqerw_linear_frontier(self)


def _cqerw_ratio(numerator, denominator):
    numerator = np.asarray(numerator, dtype=float)
    denominator = np.asarray(denominator, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.divide(
            numerator,
            denominator,
            out=np.full_like(numerator, np.nan, dtype=float),
            where=np.isfinite(numerator) & np.isfinite(denominator) & (denominator > 1e-12),
        )
    return np.clip(ratio, 0.0, 1.0)


def _cqerw_get_technical_efficiency(self):
    """Return deterministic weak CQR/CER frontier efficiency.

    This is not a StoNED/JLMS efficiency and does not use RED_MOM, RED_QLE,
    or RED_KDE.  It is a deterministic ratio against the fitted weak frontier.
    """
    y = _cqerw_as_2d(self.y)[:, 0]
    frontier = _cqerw_as_1d(self.get_frontier())
    if self.fun == FUN_PROD:
        te = _cqerw_ratio(y, frontier)
    elif self.fun == FUN_COST:
        te = _cqerw_ratio(frontier, y)
    else:
        raise ValueError("Undefined model parameters.")
    self.TE = te
    return te


def _cqerw_get_technical_efficiency_components(self):
    """Return deterministic weak CQR/CER efficiency components."""
    y = _cqerw_as_2d(self.y)[:, 0]
    frontier = _cqerw_as_1d(self.get_frontier())
    te = self.get_technical_efficiency()
    return pd.DataFrame({
        "y": y,
        "frontier": frontier,
        "residual": _cqerw_as_1d(self.get_residual()),
        "TE": te,
    })


def _cqerw_ddf_raw_distance(self):
    tools.assert_optimized(self.optimization_status)
    y = _cqerw_as_2d(self.y)
    gamma = _cqerw_as_2d(_cqerw_to_numpy(self.get_gamma()))
    return _cqerw_linear_frontier(self) - np.sum(gamma * y, axis=1)


def _cqerw_ddf_get_directional_distance(self):
    """Return nonnegative deterministic DDF distance to the fitted weak frontier."""
    raw = _cqerw_ddf_raw_distance(self)
    distance = np.where(np.isfinite(raw), np.maximum(raw, 0.0), np.nan)
    self.directional_distance = distance
    return distance


def _cqerw_ddf_components(self):
    distance = _cqerw_as_1d(self.get_directional_distance())
    x = _cqerw_as_2d(self.x)
    y = _cqerw_as_2d(self.y)
    b = _cqerw_as_2d(self.b)
    gx = np.asarray(getattr(self, "gx", []), dtype=float).reshape(-1)
    gy = np.asarray(getattr(self, "gy", []), dtype=float).reshape(-1)
    gb = np.asarray(getattr(self, "gb", []), dtype=float).reshape(-1)
    columns = {
        "raw_directional_distance": _cqerw_ddf_raw_distance(self),
        "directional_distance": distance,
    }
    scores = []
    for j, g in enumerate(gx[:x.shape[1]]):
        if g > 0:
            projected = x[:, j] - distance * g
            te = _cqerw_ratio(projected, x[:, j])
            columns[f"inefficiency_x_{j}"] = distance * g
            columns[f"TE_x_{j}"] = te
            scores.append(te)
    for k, g in enumerate(gy[:y.shape[1]]):
        if g > 0:
            projected = y[:, k] + distance * g
            te = _cqerw_ratio(y[:, k], projected)
            columns[f"inefficiency_y_{k}"] = distance * g
            columns[f"TE_y_{k}"] = te
            scores.append(te)
    for l, g in enumerate(gb[:b.shape[1]]):
        if g > 0:
            projected = b[:, l] - distance * g
            te = _cqerw_ratio(projected, b[:, l])
            columns[f"inefficiency_b_{l}"] = distance * g
            columns[f"TE_b_{l}"] = te
            scores.append(te)
    if scores:
        columns["TE"] = np.nanmin(np.column_stack(scores), axis=1)
    else:
        columns["TE"] = np.full_like(distance, np.nan, dtype=float)
    return pd.DataFrame(columns)


def _cqerw_ddf_get_technical_efficiency_components(self):
    """Return deterministic weak DDF efficiency components."""
    comp = _cqerw_ddf_components(self)
    self.TE = comp["TE"].to_numpy(dtype=float)
    return comp


def _cqerw_ddf_get_technical_efficiency(self):
    """Return binding deterministic weak DDF efficiency."""
    return self.get_technical_efficiency_components()["TE"].to_numpy(dtype=float)


for _cls in [CQRweak, CERweak, CQRweakG, CERweakG]:
    _cls.get_frontier = _cqerw_get_frontier
    _cls.get_technical_efficiency = _cqerw_get_technical_efficiency
    _cls.get_technical_efficiency_components = _cqerw_get_technical_efficiency_components

for _cls in [CQRDDFweak, CERDDFweak, CQRDDFweakG, CERDDFweakG]:
    _cls.get_directional_distance = _cqerw_ddf_get_directional_distance
    _cls.get_technical_efficiency = _cqerw_ddf_get_technical_efficiency
    _cls.get_technical_efficiency_components = _cqerw_ddf_get_technical_efficiency_components

