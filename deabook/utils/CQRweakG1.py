from ._g_common import cg_pairs, replace_constraint, make_wp_dataframe
from ..CQERDFDDFweak import CQRweak, CERweak

class CQRweakG1(CQRweak):
    def __init__(self, y, x, b, tau, cutactive, cet, fun, rts, z=None):
        data, sent, z_name = make_wp_dataframe(y, x, b, z)
        super().__init__(data, sent, tau, z=z_name, cet=cet, fun=fun, rts=rts)
        pairs = cg_pairs(cutactive)
        replace_constraint(self, 'afriat_rule', pairs, self._CQRweak__afriat_rule(), 'active afriat inequality')
        replace_constraint(self, 'disposability_rule', pairs, self._CQRweak__disposability_rule(), 'active weak disposability')

class CERweakG1(CERweak):
    def __init__(self, y, x, b, tau, cutactive, cet, fun, rts, z=None):
        data, sent, z_name = make_wp_dataframe(y, x, b, z)
        super().__init__(data, sent, tau, z=z_name, cet=cet, fun=fun, rts=rts)
        pairs = cg_pairs(cutactive)
        replace_constraint(self, 'afriat_rule', pairs, self._CQRweak__afriat_rule(), 'active afriat inequality')
        replace_constraint(self, 'disposability_rule', pairs, self._CQRweak__disposability_rule(), 'active weak disposability')
