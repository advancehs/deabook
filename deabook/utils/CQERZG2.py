from .CQERG2 import CQRG2, CERG2

class CQRZG2(CQRG2):
    def __init__(self, y, x, z, tau, cutactive, active, cet, fun, rts):
        super().__init__(y, x, tau, cutactive, active, cet, fun, rts, z=z)

class CERZG2(CERG2):
    def __init__(self, y, x, z, tau, cutactive, active, cet, fun, rts):
        super().__init__(y, x, tau, cutactive, active, cet, fun, rts, z=z)
