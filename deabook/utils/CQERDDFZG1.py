from .CQERDDFG1 import CQRDDFG1, CERDDFG1

class CQRDDFZG1(CQRDDFG1):
    def __init__(self, y, x, b, z, tau, cutactive, gy, gx, gb, fun, rts):
        super().__init__(y, x, b, tau, cutactive, gy, gx, gb, fun, rts, z=z)

class CERDDFZG1(CERDDFG1):
    def __init__(self, y, x, b, z, tau, cutactive, gy, gx, gb, fun, rts):
        super().__init__(y, x, b, tau, cutactive, gy, gx, gb, fun, rts, z=z)
