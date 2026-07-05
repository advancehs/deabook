from .CQERDDFG2 import CQRDDFG2, CERDDFG2

class CQRDDFZG2(CQRDDFG2):
    def __init__(self, y, x, b, z, tau, cutactive, active, gy, gx, gb, fun, rts):
        super().__init__(y, x, b, tau, cutactive, active, gy, gx, gb, fun, rts, z=z)

class CERDDFZG2(CERDDFG2):
    def __init__(self, y, x, b, z, tau, cutactive, active, gy, gx, gb, fun, rts):
        super().__init__(y, x, b, tau, cutactive, active, gy, gx, gb, fun, rts, z=z)
