from .CQRDDFweakG1 import CQRDDFweakG1 as _CQRDDFweakG1, CERDDFweakG1 as _CERDDFweakG1
from ._g_common import cg_pairs, replace_constraint

class CQRDDFweakG2(_CQRDDFweakG1):
    def __init__(self, y, x, b, tau, cutactive, active, activeweak, gy, gx, gb, fun, rts, z=None):
        super().__init__(y, x, b, tau, cutactive, gy, gx, gb, fun, rts, z=z)
        replace_constraint(self, 'afriat_rule', cg_pairs(cutactive, active), self._CQRDDFweak__afriat_rule(), 'active afriat inequality')
        replace_constraint(self, 'disposability_rule', cg_pairs(cutactive, activeweak), self._CQRDDFweak__disposability_rule(), 'active weak disposability')

class CERDDFweakG2(_CERDDFweakG1):
    def __init__(self, y, x, b, tau, cutactive, active, activeweak, gy, gx, gb, fun, rts, z=None):
        super().__init__(y, x, b, tau, cutactive, gy, gx, gb, fun, rts, z=z)
        replace_constraint(self, 'afriat_rule', cg_pairs(cutactive, active), self._CQRDDFweak__afriat_rule(), 'active afriat inequality')
        replace_constraint(self, 'disposability_rule', cg_pairs(cutactive, activeweak), self._CQRDDFweak__disposability_rule(), 'active weak disposability')
