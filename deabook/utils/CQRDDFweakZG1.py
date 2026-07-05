from .CQRDDFweakG1 import CQRDDFweakG1, CERDDFweakG1

class CQRDDFweakZG1(CQRDDFweakG1):
    def __init__(self, y, x, b, z, tau, cutactive, gy, gx, gb, fun, rts):
        super().__init__(y, x, b, tau, cutactive, gy, gx, gb, fun, rts, z=z)

class CERDDFweakZG1(CERDDFweakG1):
    def __init__(self, y, x, b, z, tau, cutactive, gy, gx, gb, fun, rts):
        super().__init__(y, x, b, tau, cutactive, gy, gx, gb, fun, rts, z=z)
