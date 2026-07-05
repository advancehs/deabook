from ._g_common import cg_pairs, replace_constraint
from ..CQERDFDDF import CQR, CER

class CQRG1(CQR):
    def __init__(self, y, x, tau, cutactive, cet, fun, rts, z=None):
        super().__init__(y, x, tau, z=z, cet=cet, fun=fun, rts=rts)
        replace_constraint(self, 'afriat_rule', cg_pairs(cutactive), self._CQR__afriat_rule(), 'active afriat inequality')

class CERG1(CER):
    def __init__(self, y, x, tau, cutactive, cet, fun, rts, z=None):
        super().__init__(y, x, tau, z=z, cet=cet, fun=fun, rts=rts)
        replace_constraint(self, 'afriat_rule', cg_pairs(cutactive), self._CQR__afriat_rule(), 'active afriat inequality')
