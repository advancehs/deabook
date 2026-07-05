from .CQERG1 import CQRG1, CERG1

class CQRZG1(CQRG1):
    def __init__(self, y, x, z, tau, cutactive, cet, fun, rts):
        super().__init__(y, x, tau, cutactive, cet, fun, rts, z=z)

class CERZG1(CERG1):
    def __init__(self, y, x, z, tau, cutactive, cet, fun, rts):
        super().__init__(y, x, tau, cutactive, cet, fun, rts, z=z)
