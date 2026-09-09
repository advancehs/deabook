# -*- coding: utf-8 -*-
"""Panel-data wrapper for the *static* (single-period) DEA efficiency models.

This module exposes a single class, :class:`Panel`, that turns **any** of the
book's static (single-period) efficiency estimators into a *panel* estimator.
It is deliberately method-agnostic: the efficiency method is *injected* as a
callable, so the same ``Panel`` drives every static interface in the package --
Chapter 3's radial / directional / slack models (:mod:`deabook.DEA`), the
weak-disposability family (:mod:`deabook.DEAweak`: ``ddf``, ``sbmu``,
``nddf``, ``hyper``, ``ushephard`` ...), and the material-balance models
(:mod:`deabook.DEAmb`: ``MB`` / ``MBx`` / ``MBxy`` / ``MB2``).

Design
------
A ``Panel`` declares the DMU identifier (``id``), the time variable
(``year``), the DEA ``sent`` declaration, the efficiency ``method`` to apply
(e.g. ``DEA.ccr`` or ``DEAweak.ddf`` or ``DEAmb.MB``), and the reference
``tech`` (Chapter 3 §3.4).  :meth:`solve` then loops over the observed
periods and, for each period ``t``:

* evaluates the DMUs in period ``t`` against the reference set dictated by
  ``tech`` (``CONTEMPORARY`` -> only ``t``; ``SEQUENTIAL`` -> ``<= t``;
  ``WINDOW`` -> a centred window around ``t``; ``GLOBAL`` -> every period);
* calls the injected ``method`` (``method(data, sent, evalindex=...,
  refindex=..., **kwargs).optimize(...)``);
* attaches the ``id``/``year`` columns (looked up by the *original* data index,
  never by a group-by on a possibly non-unique column) and stacks the period
  results into one long panel ``DataFrame``.

The reference-set primitives -- ``reference_years``, ``panel_years``,
``year_selector`` -- are **reused verbatim** from :mod:`deabook.mqcore`, so the
static panel frontier and the Chapter-4 Malmquist productivity frontier share
exactly the same technology definitions.  (Productivity / Malmquist indices
themselves remain out of scope here; this is the static-only companion the
book's Chapter 3 discusses.)

Robustness
----------
The static methods in the book do not share one uniform Python signature.
``Panel`` therefore adapts to each injected method instead of assuming a fixed
contract:

* **construction** -- only the keyword arguments the method actually accepts
  (``rts``, ``orient``, ``gx``/``gy``/``gb``, ``sx``/``sy``, ``level`` ...) are
  forwarded; ``evalindex``/``refindex`` are passed only if the method accepts
  them (all book static methods do);
* **solving** -- ``optimize`` is called with whatever of ``solver`` / ``email``
  / ``dmu`` the method supports, and any tuple return such as
  ``DEAmb.MB``'s ``(DataFrame, info)`` is unwrapped to its ``DataFrame``;
* **pivoting** -- :meth:`pivot` infers the efficiency column from a fallback
  list (``te`` -> ``rho`` -> ``obj`` -> ``objective_value`` -> ``beta_*``) so it
  works whether the method reports a ready-made efficiency score, a raw ``rho``,
  or a material-balance ``obj`` / ``beta`` value.
"""

from __future__ import annotations

import inspect
from typing import Callable

import pandas as pd

from .constant import (
    CET_ADDI,
    CONTEMPORARY,
    FUN_PROD,
    OPT_DEFAULT,
    OPT_LOCAL,
    RED_QLE,
    RTS_VRS1,
)
from .mqcore import (
    normalize_tech,
    normalize_window,
    panel_years,
    reference_years,
    year_selector,
)

# Columns :meth:`pivot` tries, in order, when no explicit ``value`` is given.
# Different static methods expose efficiency under different names: ready-made
# radial ``te``; directional / slack ``rho``; material-balance ``obj``; and the
# generic ``objective_value`` (e.g. NDDF's weighted ``beta`` sum).  ``beta_*``
# columns are caught separately as a last resort.
_EFFICIENCY_FALLBACK = ("te", "rho", "obj", "objective_value")


class Panel:
    """Loop any static DEA efficiency method over the periods of a panel.

    Parameters
    ----------
    data : pandas.DataFrame
        Panel data, one row per DMU-period observation.
    id : str
        Column name of the DMU identifier.
    year : str
        Column name of the time/period variable.
    sent : str
        ``"inputs=outputs"`` (or, for weak / material-balance models, the
        appropriate ``":unoutput"`` / ``"+/-"`` extension) DEA declaration.
    method : callable
        A static efficiency factory that accepts
        ``(data, sent, *, evalindex=None, refindex=None, **kwargs)`` and returns
        an optimizable model (``.optimize(...)``).  This is the *injected*
        efficiency method.  Examples from the book:

        * ``DEA.ccr`` / ``DEA.bcc`` / ``DEA.hyperbolic`` / ``DEA.ite`` /
          ``DEA.ddf`` / ``DEA.sbm`` / ``DEA.sddf``;
        * ``DEAweak.ddf`` / ``DEAweak.sbmu`` / ``DEAweak.sddf`` /
          ``DEAweak.nddf`` / ``DEAweak.hyper`` / ``DEAweak.ushephard``;
        * ``DEAmb.MB`` / ``DEAmb.MBx`` / ``DEAmb.MBxy`` / ``DEAmb.MB2``.
    tech : str, default ``CONTEMPORARY``
        Reference technology (Chapter 3 §3.4): ``CONTEMPORARY``, ``SEQUENTIAL``,
        ``WINDOW`` or ``GLOBAL``.
    window : int, default 3
        Full window width for ``tech=WINDOW`` (a centred window of this width
        around each period).
    solver : str, optional
        Solver name passed through to the method's ``optimize``.  Defaults to the
        package default (``OPT_DEFAULT`` -> the historical default solver).
    **method_kwargs
        Extra keyword arguments forwarded to ``method`` (e.g. ``orient`` for
        ``ccr``/``bcc``, ``rts`` for most methods, ``gx``/``gy``/``gb`` for the
        directional / weak models, ``sx``/``sy``/``level`` for the material-
        balance models).  Arguments the chosen ``method`` does not accept are
        silently ignored.
    """

    def __init__(
        self,
        data: pd.DataFrame,
        id: str,
        year: str,
        sent: str,
        method: Callable,
        *,
        tech: str = CONTEMPORARY,
        window: int = 3,
        solver: str | None = OPT_DEFAULT,
        **method_kwargs,
    ) -> None:
        if not isinstance(data, pd.DataFrame):
            raise TypeError("data must be a pandas.DataFrame (one row per DMU-period).")
        if not callable(method):
            raise TypeError("method must be a callable DEA factory, e.g. deabook.DEA.ccr.")

        self.data = data
        self.id = id
        self.year = year
        self.sent = sent
        self.method = method
        self.tech = normalize_tech(tech)
        self.window = normalize_window(window)
        self.solver = OPT_DEFAULT if solver is None else solver
        self.method_kwargs = dict(method_kwargs)

        # Reuse the canonical panel validation: sorts unique years, raises on
        # duplicate (id, year) pairs.
        self._years: list = panel_years(data, id, year)
        self._result: pd.DataFrame | None = None

    # -- introspection -------------------------------------------------------
    @property
    def years(self) -> list:
        """Sorted unique periods found in the panel."""
        return list(self._years)

    def __repr__(self) -> str:
        name = getattr(self.method, "__name__", repr(self.method))
        return (
            f"Panel(method={name!r}, tech={self.tech!r}, window={self.window}, "
            f"id={self.id!r}, year={self.year!r}, periods={len(self._years)})"
        )

    # -- internals -----------------------------------------------------------
    @staticmethod
    def _accepted(callable_: Callable, kwargs: dict) -> dict:
        """Keep only the keyword arguments ``callable_`` actually accepts."""
        try:
            params = inspect.signature(callable_).parameters
        except (TypeError, ValueError):
            return dict(kwargs)
        allowed = set(params)
        return {key: value for key, value in kwargs.items() if key in allowed}

    def _build_model(self, eval_sel: str, ref_sel: str):
        """Construct the injected method for one (eval, ref) period pair."""
        params = inspect.signature(self.method).parameters
        call_kwargs = dict(self._accepted(self.method, self.method_kwargs))
        if "evalindex" in params:
            call_kwargs["evalindex"] = eval_sel
        if "refindex" in params:
            call_kwargs["refindex"] = ref_sel
        if "evalindex" not in params and "refindex" not in params:
            raise ValueError(
                f"The injected method {getattr(self.method, '__name__', self.method)!r} "
                "does not accept evalindex/refindex, so it cannot be sliced into "
                "per-period panels. Use a method that supports sample/reference "
                "selection (all deabook static efficiency methods do)."
            )
        return self.method(self.data, self.sent, **call_kwargs)

    @staticmethod
    def _call_optimize(model, email: str, solver) -> pd.DataFrame:
        """Call ``model.optimize(...)`` tolerant of differing signatures/returns."""
        opt = getattr(model, "optimize", None)
        if not callable(opt):
            raise TypeError(
                "The injected method did not return an optimizable model "
                "(no optimize() callable)."
            )
        params = inspect.signature(opt).parameters
        solve_kwargs = {}
        if "solver" in params:
            solve_kwargs["solver"] = solver
        if "email" in params:
            solve_kwargs["email"] = email
        if "dmu" in params:
            # DEAmb.MB's optimize(dmu=...) triggers a verbose info print; passing
            # None suppresses it while still returning the (DataFrame, info) pair.
            solve_kwargs["dmu"] = None
        result = opt(**solve_kwargs)

        # Unwrap tuple/list returns such as DEAmb.MB's (DataFrame, info).
        if isinstance(result, (tuple, list)):
            for item in result:
                if isinstance(item, pd.DataFrame):
                    result = item
                    break
            else:
                result = result[0] if result else None

        if not isinstance(result, pd.DataFrame):
            raise TypeError(
                f"optimize() must return a DataFrame (or a tuple containing one); "
                f"got {type(result)}."
            )
        return result

    def _solve_period(self, period: object, email: str, solver) -> pd.DataFrame:
        """Evaluate one period and return its efficiency rows with id/year."""
        eval_years = [period]
        ref_years = reference_years(self._years, self._years.index(period), self.tech, self.window)
        eval_sel = year_selector(self.year, eval_years)
        ref_sel = year_selector(self.year, ref_years)

        model = self._build_model(eval_sel, ref_sel)
        result = self._call_optimize(model, email, solver)

        # Attach id/year by the ORIGINAL row index (robust to non-unique values).
        out = result.copy()
        out[self.id] = self.data.loc[out.index, self.id].values
        out[self.year] = self.data.loc[out.index, self.year].values
        return out

    # -- public API ----------------------------------------------------------
    def solve(self, email: str = OPT_LOCAL, solver: str | None = None) -> pd.DataFrame:
        """Solve every period and return a long panel ``DataFrame``.

        The returned frame has one row per DMU-period and contains the original
        ``id`` and ``year`` columns plus the efficiency columns produced by the
        injected ``method`` (e.g. ``optimization_status``, ``rho``,
        ``objective_value``, ``te``/``tei``/``teo``, ``obj``/``beta_*`` ...).
        Cached: a second call returns the stored result without re-solving.
        """
        if self._result is not None:
            return self._result

        solver = self.solver if solver is None else solver
        frames = [self._solve_period(period, email, solver) for period in self._years]
        result = pd.concat(frames, axis=0, ignore_index=False)

        # Tidy ordering: group by DMU, then by period.
        result = result.sort_values([self.id, self.year], kind="stable").reset_index(drop=True)

        self._result = result
        return result

    def pivot(self, value: str | None = None) -> pd.DataFrame:
        """Return a wide table: rows = ``id``, columns = ``year``, cells = efficiency.

        Parameters
        ----------
        value : str, optional
            Column to pivot into the cells.  When omitted, the efficiency column
            is inferred from a fallback list (``te`` -> ``rho`` -> ``obj`` ->
            ``objective_value`` -> first ``beta_*``), so the call works whether
            the method reports a ready-made efficiency score, a raw ``rho``, or a
            material-balance ``obj``/``beta`` value.

        Raises
        ------
        RuntimeError
            If :meth:`solve` has not been called yet.
        KeyError
            If no usable efficiency column can be found (pass ``value=`` then).
        """
        if self._result is None:
            raise RuntimeError("Call solve() before pivot().")
        if value is None:
            for column in _EFFICIENCY_FALLBACK:
                if column in self._result.columns:
                    value = column
                    break
        if value is None:
            beta_cols = [c for c in self._result.columns if c.startswith("beta_")]
            if beta_cols:
                value = beta_cols[0]
        if value is None or value not in self._result.columns:
            raise KeyError(
                "Could not infer an efficiency column to pivot; pass value= "
                f"explicitly. Available columns: {list(self._result.columns)}"
            )
        return self._result.pivot(index=self.id, columns=self.year, values=value)


class StoNEDPanel:
    """Panel-compatible adapter for the DDF-StoNED estimator (CNLSDDFweak + StoNED).

    StoNED is a *two-step* estimator: a CNLS frontier is first fit on a (possibly
    windowed) subset of the data, after which the residuals are decomposed into
    noise + inefficiency through :class:`deabook.StoNED`.  ``Panel`` therefore
    cannot inject the raw ``CNLSDDFweak`` factory directly -- it neither accepts
    ``evalindex``/``refindex`` nor returns an efficiency frame in one
    ``.optimize()`` call.  ``StoNEDPanel`` bridges that gap:

    * it honours the same ``(data, sent, *, evalindex, refindex, **kwargs)``
      contract ``Panel`` uses to inject a static method;
    * it fits the CNLS frontier on the *reference* subset (the window / pooled
      sample selected by ``refindex``);
    * it recovers the technical-efficiency score for every fitted observation via
      :meth:`deabook.StoNED.get_technical_efficiency_components`;
    * it returns **only** the *evaluation* subset (the current period selected by
      ``evalindex``), re-indexed to the original data index -- exactly the shape
      ``Panel`` needs to attach ``id``/``year`` and stack the periods.

    Because the CNLS fit already conditions on the windowed reference set, the
    ``tech`` (``CONTEMPORARY``/``SEQUENTIAL``/``WINDOW``/``GLOBAL``) supplied to
    ``Panel`` flows through ``year_selector`` into ``refindex``/``evalindex`` and
    produces the same rolling / sequential / global frontier concept as for the
    DEA models.  Only the **additive** DDF variant (``CNLSDDFweak``, solved with
    ``mosek``) is supported here; the multiplicative Shephard variant
    (``CNLSSDweak``) needs a nonlinear solver (``knitro``) and is out of scope.
    """

    def __init__(self, data, sent, *, evalindex=None, refindex=None,
                 rts=RTS_VRS1, gy=None, gx=None, gb=None, fun=FUN_PROD,
                 cet=CET_ADDI, red_method=RED_QLE, solver=None):
        self.data = data
        self.sent = sent
        self.evalindex = evalindex
        self.refindex = refindex
        self.rts = rts
        self.gy = [1] if gy is None else list(gy)
        self.gx = [1] if gx is None else list(gx)
        self.gb = [1] if gb is None else list(gb)
        self.fun = fun
        self.cet = cet
        self.red_method = red_method
        self._solver = solver

    def optimize(self, email=OPT_LOCAL, solver=None):
        """Fit the windowed CNLS frontier and return the evaluation-period TE frame."""
        from .CNLSSDFDDFweak import CNLSDDFweak
        from .StoNED import StoNED
        from .core import select_rows

        solver = solver or self._solver or OPT_DEFAULT
        ref_data = select_rows(self.data, self.refindex, "refindex")
        eval_data = select_rows(self.data, self.evalindex, "evalindex")

        model = CNLSDDFweak(
            ref_data, sent=self.sent, gy=self.gy, gx=self.gx, gb=self.gb,
            fun=self.fun, rts=self.rts,
        )
        model.optimize(solver=solver)

        rd = StoNED(model)
        res = rd.get_technical_efficiency_components(self.red_method)
        # The components frame comes back with a default RangeIndex; re-attach the
        # original (windowed) row indices, keep only the evaluation period, and
        # expose a lowercase `te` column so Panel.pivot can find it.
        res.index = ref_data.index
        res = res.loc[eval_data.index].rename(columns={"TE": "te"})
        return res
