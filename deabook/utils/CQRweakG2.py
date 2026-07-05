from .CQRweakG1 import CQRweakG1 as _CQRweakG1, CERweakG1 as _CERweakG1
from ._g_common import cg_pairs, replace_constraint

class CQRweakG2(_CQRweakG1):
    def __init__(self, y, x, b, tau, cutactive, active, activeweak, cet, fun, rts, z=None):
        super().__init__(y, x, b, tau, cutactive, cet, fun, rts, z=z)
        replace_constraint(self, 'afriat_rule', cg_pairs(cutactive, active), self._CQRweak__afriat_rule(), 'active afriat inequality')
        replace_constraint(self, 'disposability_rule', cg_pairs(cutactive, activeweak), self._CQRweak__disposability_rule(), 'active weak disposability')

class CERweakG2(_CERweakG1):
    def __init__(self, y, x, b, tau, cutactive, active, activeweak, cet, fun, rts, z=None):
        super().__init__(y, x, b, tau, cutactive, cet, fun, rts, z=z)
        replace_constraint(self, 'afriat_rule', cg_pairs(cutactive, active), self._CQRweak__afriat_rule(), 'active afriat inequality')
        replace_constraint(self, 'disposability_rule', cg_pairs(cutactive, activeweak), self._CQRweak__disposability_rule(), 'active weak disposability')
