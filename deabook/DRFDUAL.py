"""Directional-response-function dual shadow-price model."""

from .constant import LEFT, OPT_DEFAULT, RTS_VRS
from .dual import DirectionalResponseDualBase


class DRFDUAL(DirectionalResponseDualBase):
    """Dual of the directional response function."""

    def __init__(
        self,
        data,
        year,
        sent="inputvar=outputvar:unoutputvar",
        fenmu="unoutputvar",
        fenzi="inputvar",
        side=LEFT,
        rts=RTS_VRS,
        baseindex=None,
        refindex=None,
    ):
        super().__init__(
            data=data,
            year=year,
            sent=sent,
            fenmu=fenmu,
            fenzi=fenzi,
            side=side,
            rts=rts,
            baseindex=baseindex,
            refindex=refindex,
        )

    def optimize(self, solver=OPT_DEFAULT):
        return super().optimize(solver=solver)
