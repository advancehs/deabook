# -*- coding: utf-8 -*-
"""Shared engine for the Malmquist interfaces exposed on :mod:`deabook.DEA`.

This is a **leaf module**: it imports nothing from ``.DEA`` so that ``DEA.py``
can import it without creating a cycle (``MQDEA.py`` already imports
``DEA.py``).  All DEA solving is injected as a callable.

Reference-set semantics follow Chapter 3 §3.4 (panel-data reference sets):
``CONTEMPORARY`` (当期), ``WINDOW``, ``SEQUENTIAL`` (时序) and ``GLOBAL``
(全局), with ``CONTEMPORARY`` the default.
"""

from __future__ import annotations

from typing import Callable

import numpy as np
import pandas as pd

from .constant import (
    CONTEMPORARY,
    GLOBAL,
    OPT_DEFAULT,
    OPT_LOCAL,
    SEQUENTIAL,
    WINDOW,
)

# canonical value so downstream dispatch sees a single global token.
_TECH_CANON = {
    GLOBAL: GLOBAL,
    CONTEMPORARY: CONTEMPORARY,
    SEQUENTIAL: SEQUENTIAL,
    WINDOW: WINDOW,
}


def normalize_tech(tech: str) -> str:
    """Validate ``tech`` and return its canonical value."""
    canonical = _TECH_CANON.get(tech)
    if canonical is None:
        allowed = sorted(set(_TECH_CANON.values()))
        raise ValueError(f"tech must be one of {allowed}; got {tech!r}.")
    return canonical


def normalize_window(window: int) -> int:
    """Return ``window`` as a full width (int, >= 1)."""
    width = int(window)
    if width < 1:
        raise ValueError(f"window must be a positive integer (full width); got {window!r}.")
    return width


def reference_years(years: list, i: int, tech: str, window: int) -> list:
    """Reference years for the ``i``-th period of ``years`` (Chapter 3 §3.4).

    Positional, so it works for non-integer or non-contiguous year labels.

    - ``CONTEMPORARY`` -> ``[years[i]]``                       (§3.4-3)
    - ``WINDOW``       -> centred window, ``half = window//2``  (§3.4-4)
    - ``SEQUENTIAL``   -> ``years[: i+1]``                     (§3.4-5)
    - ``GLOBAL``       -> all years                            (§3.4-6)
    """
    count = len(years)
    if tech == CONTEMPORARY:
        return [years[i]]
    if tech == GLOBAL:
        return list(years)
    if tech == SEQUENTIAL:
        return list(years[: i + 1])
    if tech == WINDOW:
        half = normalize_window(window) // 2
        return list(years[max(0, i - half): min(count - 1, i + half) + 1])
    raise ValueError(f"Unsupported tech {tech!r}.")


def year_selector(year_col: str, values: list) -> str:
    """Robust legacy selector ``Year=[2016,2017]`` parsed by ``core.select_rows``.

    Uses ``repr`` so string year labels are quoted correctly.
    """
    return f"{year_col}=[{','.join(repr(v) for v in values)}]"


def panel_years(data: pd.DataFrame, id: str, year: str) -> list:
    """Sorted unique years as *Python* scalars; raise on missing/empty panel."""
    for column in (id, year):
        if column not in data.columns:
            raise ValueError(f"Column '{column}' not found in data.")
    if data.duplicated(subset=[id, year]).any():
        raise ValueError(f"Duplicate ({id!r}, {year!r}) pairs in the panel.")
    return sorted(pd.Series(data[year]).drop_duplicates().tolist())


def build_distance_frame(
    solve_distance: Callable[[str, str, str, str], pd.DataFrame],
    work: pd.DataFrame,
    years: list,
    tech: str,
    window: int,
    year_col: str,
    *,
    id_col: str,
    email: str,
    solver: str,
    bundle: str = "E",
) -> pd.DataFrame:
    """Compute ``D11``/``D12``/``D21`` (one efficiency bundle) for every period.

    ``solve_distance(eval_sel, ref_sel, email, solver)`` must return a DataFrame
    that contains a column named ``bundle``, indexed by the rows of ``work`` for
    ``eval_sel``.  ``D12`` is ``NaN`` for the first period, ``D21`` for the last.

    Under ``GLOBAL`` only ``D11`` is solved and reused for ``D12``/``D21``, so the
    generic Malmquist formulas collapse to ``mqpi = D11/D11_prev`` and
    ``MTECHCH = 1`` without a separate branch.
    """
    tech = normalize_tech(tech)
    window = normalize_window(window)
    last = len(years) - 1
    parts = {"D11": [], "D12": [], "D21": []}
    for i, period in enumerate(years):
        eval_sel = year_selector(year_col, [period])

        def solve_for(refs: list) -> pd.Series:
            frame = solve_distance(eval_sel, year_selector(year_col, refs), email, solver)
            return frame[bundle]

        parts["D11"].append(solve_for(reference_years(years, i, tech, window)))
        if tech == GLOBAL:
            parts["D12"].append(parts["D11"][-1])
            parts["D21"].append(parts["D11"][-1])
        else:
            if i > 0:
                parts["D12"].append(solve_for(reference_years(years, i - 1, tech, window)))
            if i < last:
                parts["D21"].append(solve_for(reference_years(years, i + 1, tech, window)))

    out = pd.DataFrame(index=work.index)
    out[id_col] = work[id_col].values
    for component, series_list in parts.items():
        if series_list:
            out[component] = pd.concat(series_list).reindex(work.index)
        else:
            out[component] = np.nan
    return out


def _shift(frame: pd.DataFrame, id: str, column: str) -> pd.Series:
    """Previous-period value of ``column`` within each ``id`` group."""
    return frame.groupby(id, sort=False)[column].shift(1)


def _sqrt_clip(series: pd.Series) -> pd.Series:
    return np.sqrt(series.clip(lower=0))


def two_factor(frame: pd.DataFrame, id: str, suffix: str = "") -> pd.DataFrame:
    """Two-factor Malmquist decomposition (eq-4-4/4-5, 4-9/4-10, 4-22/4-24)."""
    d11, d12, d21 = frame["D11"], frame["D12"], frame["D21"]
    d11_prev = _shift(frame, id, "D11")
    d21_prev = _shift(frame, id, "D21")
    out = pd.DataFrame(index=frame.index)
    out[f"MQ{suffix}"] = _sqrt_clip((d12 / d11_prev) * (d11 / d21_prev))
    out[f"MEFFCH{suffix}"] = d11 / d11_prev
    out[f"MTECHCH{suffix}"] = _sqrt_clip((d12 / d11) * (d11_prev / d21_prev))
    return out


def fgnz_decompose(vrs: pd.DataFrame, crs: pd.DataFrame, id: str) -> pd.DataFrame:
    """FGNZ (1994) three-factor decomposition (eq-4-11..4-14)."""
    v11, v11_prev = vrs["D11"], _shift(vrs, id, "D11")
    c11, c11_prev = crs["D11"], _shift(crs, id, "D11")
    c12, c21 = crs["D12"], crs["D21"]
    c21_prev = _shift(crs, id, "D21")

    peffch = v11 / v11_prev
    seffch = (c11 / v11) / (c11_prev / v11_prev)
    techch = _sqrt_clip((c12 / c11) * (c11_prev / c21_prev))

    out = pd.DataFrame(index=vrs.index)
    out["PEFFCH"] = peffch
    out["SEFFCH"] = seffch
    out["TECHCH"] = techch
    out["MQ"] = peffch * seffch * techch
    out["MEFFCH"] = peffch * seffch
    out["MTECHCH"] = techch
    return out


def rd_decompose(vrs: pd.DataFrame, crs: pd.DataFrame, id: str) -> pd.DataFrame:
    """Ray & Desli (1997) decomposition (eq-4-15..4-17).

    ``EFFCH`` equals FGNZ's ``PEFFCH``; ``TECHCH_RD = M_o / EFFCH`` (=
    ``SEFFCH * TECHCH``), so the Malmquist index equals FGNZ's ``MQ`` exactly and
    only the split differs.
    """
    fgnz = fgnz_decompose(vrs, crs, id)
    out = pd.DataFrame(index=vrs.index)
    out["EFFCH"] = fgnz["PEFFCH"]
    out["TECHCH_RD"] = fgnz["MQ"] / fgnz["PEFFCH"]
    out["MQ"] = fgnz["MQ"]
    return out


def global_decompose(vrs: pd.DataFrame, crs: pd.DataFrame, id: str) -> pd.DataFrame:
    """Pastor & Lovell (2005) global-frontier three-factor decomposition.

    Same three-factor split as FGNZ, plus ``mqpi = D11^c / D11^c_prev``.
    """
    out = fgnz_decompose(vrs, crs, id)
    out["mqpi"] = crs["D11"] / _shift(crs, id, "D11")
    return out


def pl_global_decompose(
    cont_vrs: pd.DataFrame,
    cont_crs: pd.DataFrame,
    glob_crs: pd.DataFrame,
    id: str,
) -> pd.DataFrame:
    """Pastor & Lovell (2005) three-factor decomposition of the *global* index.

    Decomposes the global Malmquist index ``mqpi = D^G(x^t)/D^G(x^{t-1})`` into

    - ``PEFFCH``  : contemporary VRS efficiency change
                   ``D^t,v(x^t) / D^{t-1,v}(x^{t-1})``;
    - ``SEFFCH``  : contemporary scale-efficiency change
                   ``(D^t,c/D^t,v) / (D^{t-1,c}/D^{t-1,v})``;
    - ``TECHCH``  : best-practice-gap change
                   ``BPG_t / BPG_{t-1}``, ``BPG_t = D^G(x^t)/D^t(x^t)``;

    so that ``mqpi = PEFFCH * SEFFCH * TECHCH`` exactly, with
    ``MEFFCH = PEFFCH * SEFFCH`` and ``MTECHCH = TECHCH`` matching the two-factor
    Pastor-Lovell decomposition of :func:`two_factor` under ``GLOBAL``.
    """
    v11, v11_prev = cont_vrs["D11"], _shift(cont_vrs, id, "D11")
    c11, c11_prev = cont_crs["D11"], _shift(cont_crs, id, "D11")
    g11, g11_prev = glob_crs["D11"], _shift(glob_crs, id, "D11")

    peffch = v11 / v11_prev
    seffch = (c11 / v11) / (c11_prev / v11_prev)
    bpg_frame = pd.DataFrame({id: cont_crs[id], "BPG": g11 / c11}, index=cont_crs.index)
    techch = bpg_frame["BPG"] / _shift(bpg_frame, id, "BPG")

    out = pd.DataFrame(index=cont_vrs.index)
    out["PEFFCH"] = peffch
    out["SEFFCH"] = seffch
    out["TECHCH"] = techch
    out["MQ"] = peffch * seffch * techch
    out["MEFFCH"] = peffch * seffch
    out["MTECHCH"] = techch
    out["mqpi"] = g11 / g11_prev
    return out


class MalmquistModel:
    """Lazy Malmquist model: built by a ``DEA.mq*`` factory, solved by
    :meth:`optimize`.

    Mirrors the Chapter-3 lifecycle (``DEA.ccr(...).optimize(...)``) instead of
    the legacy ``MQDEA(...)`` eager ``__init__``.
    """

    def __init__(
        self,
        data: pd.DataFrame,
        id: str,
        year: str,
        *,
        tech: str,
        window: int,
        build: Callable[[pd.DataFrame, list, str, int, str, str, str], pd.DataFrame],
        kind: str,
        orient: str | None = None,
        verbose: bool = False,
    ) -> None:
        self.data = data
        self.id = id
        self.year = year
        self.tech = normalize_tech(tech)
        self.window = normalize_window(window)
        self._build = build
        self.kind = kind
        self.orient = orient
        self.verbose = verbose
        self._result: pd.DataFrame | None = None
        self._years = panel_years(data, id, year)

    @property
    def years(self) -> list:
        return list(self._years)

    def __repr__(self) -> str:
        return (
            f"MalmquistModel(kind={self.kind!r}, tech={self.tech!r}, "
            f"window={self.window}, periods={len(self._years)})"
        )

    def optimize(self, email: str = OPT_LOCAL, solver: str = OPT_DEFAULT) -> pd.DataFrame:
        """Solve (cached) and return a DataFrame indexed like the input ``data``."""
        if self._result is not None:
            return self._result

        work = self.data.sort_values([self.id, self.year])
        result = self._build(work, self._years, self.tech, self.window, self.year, email, solver)
        if self.verbose:
            print(f"[{self.kind}] tech={self.tech} solved {len(self._years)} periods.")

        out = self.data.copy()
        for column in result.columns:
            out[column] = result[column].reindex(self.data.index)
        self._result = out
        return out
