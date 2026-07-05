from .CQRDDFweakG2 import CQRDDFweakG2, CERDDFweakG2

class CQRDDFweakZG2(CQRDDFweakG2):
    def __init__(self, y, x, b, z, tau, cutactive, active, activeweak, gy, gx, gb, fun, rts):
        super().__init__(y, x, b, tau, cutactive, active, activeweak, gy, gx, gb, fun, rts, z=z)

class CERDDFweakZG2(CERDDFweakG2):
    def __init__(self, y, x, b, z, tau, cutactive, active, activeweak, gy, gx, gb, fun, rts):
        super().__init__(y, x, b, tau, cutactive, active, activeweak, gy, gx, gb, fun, rts, z=z)
