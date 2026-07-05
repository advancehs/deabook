from .CQRweakG1 import CQRweakG1, CERweakG1

class CQRweakZG1(CQRweakG1):
    def __init__(self, y, x, b, z, tau, cutactive, cet, fun, rts):
        super().__init__(y, x, b, tau, cutactive, cet, fun, rts, z=z)

class CERweakZG1(CERweakG1):
    def __init__(self, y, x, b, z, tau, cutactive, cet, fun, rts):
        super().__init__(y, x, b, tau, cutactive, cet, fun, rts, z=z)
