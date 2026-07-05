from ._g_common import cg_pairs, replace_constraint, make_wp_dataframe
from ..CQERDFDDFweak import CQRDDFweak, CERDDFweak

class CQRDDFweakG1(CQRDDFweak):
    def __init__(self, y, x, b, tau, cutactive, gy, gx, gb, fun, rts, z=None):
        data, sent, z_name = make_wp_dataframe(y, x, b, z)
        super().__init__(data, sent, tau, z=z_name, gy=gy, gx=gx, gb=gb, fun=fun, rts=rts)
        pairs = cg_pairs(cutactive)
        replace_constraint(self, 'afriat_rule', pairs, self._CQRDDFweak__afriat_rule(), 'active afriat inequality')
        replace_constraint(self, 'disposability_rule', pairs, self._CQRDDFweak__disposability_rule(), 'active weak disposability')

class CERDDFweakG1(CERDDFweak):
    def __init__(self, y, x, b, tau, cutactive, gy, gx, gb, fun, rts, z=None):
        data, sent, z_name = make_wp_dataframe(y, x, b, z)
        super().__init__(data, sent, tau, z=z_name, gy=gy, gx=gx, gb=gb, fun=fun, rts=rts)
        pairs = cg_pairs(cutactive)
        replace_constraint(self, 'afriat_rule', pairs, self._CQRDDFweak__afriat_rule(), 'active afriat inequality')
        replace_constraint(self, 'disposability_rule', pairs, self._CQRDDFweak__disposability_rule(), 'active weak disposability')
