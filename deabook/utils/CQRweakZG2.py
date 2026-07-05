from .CQRweakG2 import CQRweakG2, CERweakG2

class CQRweakZG2(CQRweakG2):
    def __init__(self, y, x, b, z, tau, cutactive, active, activeweak, cet, fun, rts):
        super().__init__(y, x, b, tau, cutactive, active, activeweak, cet, fun, rts, z=z)

class CERweakZG2(CERweakG2):
    def __init__(self, y, x, b, z, tau, cutactive, active, activeweak, cet, fun, rts):
        super().__init__(y, x, b, tau, cutactive, active, activeweak, cet, fun, rts, z=z)
