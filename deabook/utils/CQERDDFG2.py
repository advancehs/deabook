from .CQERDDFG1 import CQRDDFG1 as _CQRDDFG1, CERDDFG1 as _CERDDFG1
from ._g_common import cg_pairs, replace_constraint

class CQRDDFG2(_CQRDDFG1):
    def __init__(self, y, x, b, tau, cutactive, active, gy, gx, gb, fun, rts, z=None):
        super().__init__(y, x, b, tau, cutactive, gy, gx, gb, fun, rts, z=z)
        replace_constraint(self, 'afriat_rule', cg_pairs(cutactive, active), self._CQRDDF__afriat_rule(), 'active afriat inequality')

class CERDDFG2(_CERDDFG1):
    def __init__(self, y, x, b, tau, cutactive, active, gy, gx, gb, fun, rts, z=None):
        super().__init__(y, x, b, tau, cutactive, gy, gx, gb, fun, rts, z=z)
        replace_constraint(self, 'afriat_rule', cg_pairs(cutactive, active), self._CQRDDF__afriat_rule(), 'active afriat inequality')
