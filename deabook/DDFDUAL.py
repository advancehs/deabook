"""Directional-distance-function dual shadow-price model."""

from .constant import OPT_DEFAULT, RTS_VRS
from .dual import DirectionalDistanceDualBase


class DDFDUAL(DirectionalDistanceDualBase):
    """Dual of the directional distance function.

    The public constructor and output columns are kept compatible with the
    legacy implementation, while the shared dual-model lifecycle lives in
    :class:`deabook.dual.DirectionalDistanceDualBase`.
    """

    def __init__(
        self,
        data,
        sent="inputvar=outputvar:unoutputvar",
        gy=[1],
        gx=[1],
        gb=[1],
        rts=RTS_VRS,
        baseindex=None,
        refindex=None,
    ):
        super().__init__(
            data=data,
            sent=sent,
            gy=gy,
            gx=gx,
            gb=gb,
            rts=rts,
            baseindex=baseindex,
            refindex=refindex,
        )

    def optimize(self, solver=OPT_DEFAULT):
        return super().optimize(solver=solver)
