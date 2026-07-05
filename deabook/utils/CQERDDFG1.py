from ._g_common import cg_pairs, replace_constraint
from ..CQERDFDDF import CQRDDF, CERDDF

class CQRDDFG1(CQRDDF):
    def __init__(self, y, x, b, tau, cutactive, gy, gx, gb, fun, rts, z=None):
        super().__init__(y, x, tau, b=b, z=z, gy=gy, gx=gx, gb=gb, fun=fun, rts=rts)
        replace_constraint(self, 'afriat_rule', cg_pairs(cutactive), self._CQRDDF__afriat_rule(), 'active afriat inequality')

class CERDDFG1(CERDDF):
    def __init__(self, y, x, b, tau, cutactive, gy, gx, gb, fun, rts, z=None):
        super().__init__(y, x, tau, b=b, z=z, gy=gy, gx=gx, gb=gb, fun=fun, rts=rts)
        replace_constraint(self, 'afriat_rule', cg_pairs(cutactive), self._CQRDDF__afriat_rule(), 'active afriat inequality')
