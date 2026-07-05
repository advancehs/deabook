"""Teaching-friendly DEA entry points built on the object-oriented core."""

from __future__ import annotations

from collections.abc import Sequence

import pandas as pd

from .DEA import DEA2
from .constant import OPT_DEFAULT, RTS_CRS, RTS_VRS1
from .core import parse_formula, prepare_dea_data


_RTS_ALIASES = {
    "crs": RTS_CRS,
    "vrs": RTS_VRS1,
    RTS_CRS: RTS_CRS,
    RTS_VRS1: RTS_VRS1,
}


class DEAModel:
    """Small object-oriented wrapper for the book's introductory DEA examples.

    Parameters
    ----------
    formula : str
        Teaching-style formula, with outputs on the left and inputs on the
        right, for example ``"Y = K L"``.
    dataframe : pandas.DataFrame
        DEA data. Each row is one DMU.
    orient : {"oo", "io"}, default "oo"
        Output orientation or input orientation.
    rts : {"crs", "vrs", RTS_CRS, RTS_VRS1}, default "crs"
        Returns-to-scale assumption.
    """

    def __init__(self, formula, dataframe, orient="oo", rts="crs"):
        if orient not in {"oo", "io"}:
            raise ValueError("orient must be 'oo' or 'io'.")
        if rts not in _RTS_ALIASES:
            raise ValueError("rts must be 'crs', 'vrs', RTS_CRS or RTS_VRS1.")
        self.formula = formula
        self.dataframe = dataframe
        self.orient = orient
        self.rts = _RTS_ALIASES[rts]
        self.spec = parse_formula(formula)
        self.sent = f"{' '.join(self.spec.x)}={' '.join(self.spec.y)}"
        self.model = None
        self.results = None

    def solve(self, eval_query=None, ref_query=None, solver=OPT_DEFAULT):
        """Solve the DEA model.

        Parameters
        ----------
        eval_query : str or None
            Optional pandas query selecting evaluated DMUs.
        ref_query : str or None
            Optional pandas query selecting reference DMUs.
        solver : str or None
            Pyomo solver name.

        Returns
        -------
        pandas.DataFrame
            Result table with ``theta`` and ``te`` columns.
        """
        gx = [1] * len(self.spec.x) if self.orient == "io" else [0] * len(self.spec.x)
        gy = [1] * len(self.spec.y) if self.orient == "oo" else [0] * len(self.spec.y)
        self.model = DEA2(
            self.dataframe,
            sent=self.sent,
            gx=gx,
            gy=gy,
            rts=self.rts,
            baseindex=eval_query,
            refindex=ref_query,
        )
        result = self.model.optimize(solver=solver)
        result = result.rename(columns={"rho": "theta"})
        self.results = result
        return result


def validate_dea_inputs(dataframe, input_cols, output_cols, undesirable_cols=None, orient="oo", rts="crs"):
    """Validate a DataFrame-based DEA input specification.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        Source DEA data.
    input_cols, output_cols : sequence of str
        Input and desirable-output variable names.
    undesirable_cols : sequence of str or None
        Optional undesirable-output names.
    orient : {"oo", "io"}
        Orientation used by introductory radial DEA examples.
    rts : {"crs", "vrs", RTS_CRS, RTS_VRS1}
        Returns-to-scale assumption.

    Returns
    -------
    bool
        ``True`` when validation succeeds.
    """
    if not isinstance(dataframe, pd.DataFrame):
        raise TypeError("dataframe must be a pandas.DataFrame.")
    if not isinstance(input_cols, Sequence) or isinstance(input_cols, str):
        raise TypeError("input_cols must be a sequence of column names.")
    if not isinstance(output_cols, Sequence) or isinstance(output_cols, str):
        raise TypeError("output_cols must be a sequence of column names.")
    bad_cols = [] if undesirable_cols is None else list(undesirable_cols)
    if orient not in {"oo", "io"}:
        raise ValueError("orient must be 'oo' or 'io'.")
    if rts not in _RTS_ALIASES:
        raise ValueError("rts must be 'crs', 'vrs', RTS_CRS or RTS_VRS1.")
    formula = f"{' '.join(input_cols)}={' '.join(output_cols)}"
    if bad_cols:
        formula += f":{' '.join(bad_cols)}"
    prepare_dea_data(dataframe, formula, input_side="left")
    if dataframe.shape[0] < 2:
        raise ValueError("DEA requires at least two DMUs.")
    return True


def simple_dea_solve(dataframe, input_cols, output_cols, orient="oo", rts="crs", solver=OPT_DEFAULT):
    """Solve a simple radial DEA model from column lists.

    Parameters
    ----------
    dataframe : pandas.DataFrame
        Source DEA data.
    input_cols, output_cols : sequence of str
        Input and output columns.
    orient : {"oo", "io"}, default "oo"
        Output or input orientation.
    rts : {"crs", "vrs"}, default "crs"
        Returns-to-scale assumption.
    solver : str or None
        Pyomo solver name.

    Returns
    -------
    pandas.DataFrame
        Result table with ``theta`` and ``te``.
    """
    validate_dea_inputs(dataframe, input_cols, output_cols, orient=orient, rts=rts)
    formula = f"{' '.join(output_cols)}={' '.join(input_cols)}"
    return DEAModel(formula, dataframe, orient=orient, rts=rts).solve(solver=solver)
