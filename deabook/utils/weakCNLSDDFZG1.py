from ._g_common import cg_pairs, replace_constraint
from ..CNLSSDFDDFweak import CNLSDDFweak
from pyomo.environ import Constraint

class weakCNLSDDFZG1(CNLSDDFweak):
    def __init__(self, data, sent, z, cutactive, gy, gx, gb, fun, rts, baseindex=None, refindex=None):
        super().__init__(data, sent, z=z, gy=gy, gx=gx, gb=gb, fun=fun, rts=rts)
        pairs = cg_pairs(cutactive)
        replace_constraint(self, 'afriat_rule', pairs, self._CNLSDDFweak__afriat_rule(), 'active afriat inequality')
        if hasattr(self.__model__, 'red_factor_rule'):
            self.__model__.del_component(getattr(self.__model__, 'red_factor_rule'))
            def red_rule(model, i, h):
                return model.alpha[h] + sum(model.delta[h, j] * self.x[i][j] for j in model.J) >= 0
            replace_constraint(self, 'red_factor_rule', pairs, red_rule, 'active weak disposability')

    # Compatibility accessors for the outer weakCNLSDDFG loop.
    # CNLSDDFweak names input coefficients delta and bad-output coefficients kappa,
    # while weakCNLSDDFG expects beta/gamma/delta.
    def get_beta(self):
        return self.get_delta()

    def get_delta(self):
        return self.get_kappa()

class weakCNLSDDFZG2(weakCNLSDDFZG1):
    def __init__(self, data, sent, z, cutactive, active, activeweak, gy, gx, gb, fun, rts, baseindex=None, refindex=None):
        super().__init__(data, sent, z, cutactive, gy, gx, gb, fun, rts, baseindex=baseindex, refindex=refindex)
        replace_constraint(self, 'afriat_rule', cg_pairs(cutactive, active), self._CNLSDDFweak__afriat_rule(), 'active afriat inequality')
        if hasattr(self.__model__, 'red_factor_rule'):
            def red_rule(model, i, h):
                return model.alpha[h] + sum(model.delta[h, j] * self.x[i][j] for j in model.J) >= 0
            replace_constraint(self, 'red_factor_rule', cg_pairs(cutactive, activeweak), red_rule, 'active weak disposability')
