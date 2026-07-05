from .CQERG1 import CQRG1 as _CQRG1, CERG1 as _CERG1
from ._g_common import cg_pairs, replace_constraint

class CQRG2(_CQRG1):
    def __init__(self, y, x, tau, cutactive, active, cet, fun, rts, z=None):
        super().__init__(y, x, tau, cutactive, cet, fun, rts, z=z)
        replace_constraint(self, 'afriat_rule', cg_pairs(cutactive, active), self._CQR__afriat_rule(), 'active afriat inequality')

class CERG2(_CERG1):
    def __init__(self, y, x, tau, cutactive, active, cet, fun, rts, z=None):
        super().__init__(y, x, tau, cutactive, cet, fun, rts, z=z)
        replace_constraint(self, 'afriat_rule', cg_pairs(cutactive, active), self._CQR__afriat_rule(), 'active afriat inequality')
