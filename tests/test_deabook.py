from deabook import PANEL
from deabook.constant import CONTEMPORARY, GLOBAL, SEQUENTIAL, WINDOW
#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for the ``deabook`` modules actually used by the book (chapters 3–10).

These tests verify *dependency completeness* of the surviving module set: every
public API the manuscript calls must import, construct, and (where an LP solver
is available locally) solve without a missing-class / missing-function error.

Coverage by chapter
-------------------
* Chapter 3–5, 8  -> ``deabook.DEA``      (ccr, bcc, hyperbolic, ite, ddf, sbm,
                                            mq, mqddf, mqfgnz, mqrd, gmqdecom,
                                            ep, epddf)
* Chapter 6–8     -> ``deabook.DEAweak``  (ushephard, hyper, ddf, nddf, dep,
                                            depddf, lddf, sbmu, and the
                                            DDFweak / NDDFweak class family
                                            migrated back from the backup on
                                            2026-09-05)
* Chapter 9       -> ``deabook.DEAmb``    (MB)
* Chapter 10      -> ``deabook.CNLSSDFDDFweak`` (CNLSSDweak, CNLSDDFweak)
                     ``deabook.CQERDFDDFweak``  (CQRweak, CQRDDFweak -- needs
                     the optional ``amplpy`` solver)

Environment notes
------------------
* An LP solver (``highs``/``cbc``) is required for the DEA / DEAweak optimising
  tests. They are skipped if none is installed.
* ``CNLSSDweak`` / ``CNLSDDFweak`` construct without an external solver, but
  their ``.optimize()`` needs ``mosek``/``knitro``; that part is skipped when
  those solvers are absent.
* ``CQRweak`` / ``CQRDDFweak`` import ``amplpy`` at module load time, so the
  whole test is skipped when ``amplpy`` is not installed. This is an *optional
  external dependency*, not a missing package module, and is therefore reported
  as a skip rather than a failure.
"""

import importlib.util
import numpy as np
import pandas as pd
import pytest
import pyomo.environ  # noqa: F401  (registers solver plugins for detection)


# ---------------------------------------------------------------------------
# Solver detection
# ---------------------------------------------------------------------------
def _solver_available(name):
    try:
        return bool(pyomo.environ.SolverFactory(name).available(exception_flag=False))
    except Exception:
        return False


def _local_solver():
    for name in ("highs", "appsi_highs", "cbc"):
        if _solver_available(name):
            return name
    return None


SOLVER = _local_solver()
HAS_CNLS_SOLVER = _solver_available("mosek") or _solver_available("knitro")
HAS_AMPLPY = importlib.util.find_spec("amplpy") is not None


# ---------------------------------------------------------------------------
# Synthetic but *feasible* data (random data can make the LP infeasible for
# some DMUs, which is a data artefact, not a dependency problem).
# ---------------------------------------------------------------------------
def _feasible_df(n=12, seed=2):
    rng = np.random.default_rng(seed)
    K = rng.uniform(10, 50, n)
    L = rng.uniform(20, 60, n)
    E = rng.uniform(5, 15, n)
    Y = 2 * K + 3 * L + 5 * E + rng.uniform(0, 20, n)      # monotone in inputs
    CO2 = 0.5 * E + rng.uniform(0, 2, n)                   # pollution tied to E
    firms = [1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4, 4]
    years = [2010, 2011, 2012] * 4
    return pd.DataFrame(
        {"K": K, "L": L, "E": E, "Y": Y, "CO2": CO2,
         "firm": firms[:n], "year": years[:n]}
    )


# ===========================================================================
# Chapter 3–5, 8 : deabook.DEA
# ===========================================================================
from deabook.DEA import (  # noqa: E402
    ccr, bcc, hyperbolic, ite, ddf, sbm, sddf,
    mq, mqddf, mqfgnz, mqrd, gmqdecom, ep, epddf,
)
from deabook.constant import RTS_VRS1, RTS_VRS2, CET_MULT, CET_ADDI, FUN_PROD


@pytest.mark.skipif(SOLVER is None, reason="No local LP solver (highs/cbc) available")
class TestDEA:
    def test_ccr(self):
        df = _feasible_df()
        res = ccr(df, sent="K L E=Y", orient="oo").optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)
        assert res["rho"].notna().all()

    def test_bcc(self):
        df = _feasible_df()
        res = bcc(df, sent="K L E=Y", orient="oo").optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)
        assert res["rho"].notna().all()

    def test_hyperbolic(self):
        df = _feasible_df()
        res = hyperbolic(df, sent="K L E=Y").optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)

    def test_ite(self):
        df = _feasible_df()
        res = ite(df, sent="K L E=Y", rts=RTS_CRS).optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)

    def test_ddf(self):
        df = _feasible_df()
        res = ddf(df, sent="K L E=Y", gx=[0, 0, 0], gy=[1]).optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)

    def test_sbm(self):
        df = _feasible_df()
        res = sbm(df, sent="K L E=Y", gx=[1, 1, 1], gy=[1]).optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)
        assert res["rho"].notna().all()

    def test_sddf(self):
        # Färe-Grosskopf (2010) SD model (no undesirable outputs).
        df = _feasible_df()
        m = sddf(df, sent="K L E=Y", rts=RTS_VRS1)
        res = m.optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)
        assert "objective_value" in res.columns and "te" in res.columns
        assert res["objective_value"].notna().all()
        assert any(c.startswith("slack_x_") for c in res.columns)
        assert any(c.startswith("slack_y_") for c in res.columns)
        # A dominated DMU must be flagged inefficient (objective > 0).
        dominated = pd.DataFrame({
            "K": [10.0, 20.0, 30.0], "L": [20.0, 40.0, 60.0],
            "E": [5.0, 10.0, 8.0], "Y": [100.0, 200.0, 150.0],
        })
        res_d = sddf(dominated, sent="K L E=Y", rts=RTS_VRS1).optimize(solver=SOLVER)
        assert res_d["objective_value"].iloc[2] > 0
        # Undesirable outputs must be rejected by DEA.sddf.
        try:
            sddf(df, sent="K L E=Y:CO2")
            raise AssertionError("DEA.sddf should reject undesirable outputs")
        except ValueError:
            pass

    def test_mq(self):
        df = _feasible_df()
        res = mq(df, id="firm", year="year", sent="K L E=Y", orient="oo",
                tech=CONTEMPORARY).optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)

    def test_mqddf(self):
        df = _feasible_df()
        res = mqddf(df, id="firm", year="year", sent="K L E=Y", gy=[1], gx=[0, 0, 0],
                    tech=CONTEMPORARY).optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)

    def test_mqfgnz(self):
        df = _feasible_df()
        res = mqfgnz(df, id="firm", year="year", sent="K L E=Y", orient="oo",
                     tech=CONTEMPORARY).optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)

    def test_mqrd(self):
        df = _feasible_df()
        res = mqrd(df, id="firm", year="year", sent="K L E=Y", orient="oo",
                  tech=CONTEMPORARY).optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)

    def test_gmqdecom(self):
        df = _feasible_df()
        res = gmqdecom(df, id="firm", year="year", sent="K L E=Y", orient="oo",
                       tech=GLOBAL).optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)

    def test_ep(self):
        df = _feasible_df()
        res = ep(df, id="firm", year="year", sent="K L E=Y", rts=RTS_CRS,
                 tech=CONTEMPORARY).optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)

    def test_epddf(self):
        df = _feasible_df()
        res = epddf(df, id="firm", year="year", sent="K L E=Y", gy=[0], gx=[0, 0, 1],
                    tech=CONTEMPORARY).optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)


# ===========================================================================
# Chapter 6–8 : deabook.DEAweak  (+ migrated DDFweak / NDDFweak family)
# ===========================================================================
from deabook.DEAweak import (  # noqa: E402
    DEAweak, DEAweakRho,
    DDFweak, DDFweakRho, NDDFweak, NDDFweakRho,
    ushephard, hyper, ddf as weak_ddf, nddf,
    dep, depddf, lddf, sbmu, sddf as weak_sddf,
)


@pytest.mark.skipif(SOLVER is None, reason="No local LP solver (highs/cbc) available")
class TestDEAweak:
    def test_ushephard(self):
        df = _feasible_df()
        res = ushephard(df, sent="K L E=Y:CO2", rts=RTS_CRS).optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)

    def test_hyper(self):
        df = _feasible_df()
        res = hyper(df, sent="K L E=Y:CO2", rts=RTS_CRS).optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)

    def test_ddf(self):
        df = _feasible_df()
        res = weak_ddf(df, sent="K L E=Y:CO2", rts=RTS_CRS).optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)

    def test_nddf(self):
        df = _feasible_df()
        res = nddf(df, sent="K L E=Y:CO2", rts=RTS_CRS).optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)

    def test_dep(self):
        df = _feasible_df()
        res = dep(df, id="firm", year="year", sent="K L E=Y:CO2", rts=RTS_CRS,
                  tech=CONTEMPORARY).optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)

    def test_depddf(self):
        df = _feasible_df()
        res = depddf(df, id="firm", year="year", sent="K L E=Y:CO2", rts=RTS_CRS,
                     tech=CONTEMPORARY).optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)

    def test_lddf(self):
        df = _feasible_df()
        res = lddf(df, id="firm", year="year", sent="K L E=Y:CO2", orient="oo",
                   rts=RTS_CRS, tech=CONTEMPORARY).optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)

    def test_sbmu(self):
        df = _feasible_df()
        m = sbmu(df, sent="K L E=Y:CO2", rts=RTS_VRS1)
        res = m.optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)

    def test_sddf(self):
        # Pham-Zelenyuk (2018) SD model with undesirable outputs.
        df = _feasible_df()
        m = weak_sddf(df, sent="K L E=Y:CO2", rts=RTS_VRS1)
        res = m.optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(df)
        assert "objective_value" in res.columns and "te" in res.columns
        assert any(c.startswith("slack_b_") for c in res.columns)
        # A dominated DMU must be flagged inefficient (objective > 0).
        dominated = pd.DataFrame({
            "K": [10.0, 20.0, 30.0], "L": [20.0, 40.0, 60.0],
            "E": [5.0, 10.0, 8.0], "Y": [100.0, 200.0, 150.0],
            "CO2": [2.0, 4.0, 10.0],
        })
        res_d = weak_sddf(dominated, sent="K L E=Y:CO2", rts=RTS_VRS1).optimize(solver=SOLVER)
        assert res_d["objective_value"].iloc[2] > 0
        # A model without an undesirable output must be rejected.
        try:
            weak_sddf(df, sent="K L E=Y")
            raise AssertionError("DEAweak.sddf should require undesirable outputs")
        except ValueError:
            pass

    def test_migrated_classes_present(self):
        """The DDFweak/NDDFweak family was migrated back from the backup on
        2026-09-05; its presence is what repairs DEAweak.ddf/nddf/dep/depddf/lddf.
        """
        for cls in (DDFweak, DDFweakRho, NDDFweak, NDDFweakRho):
            assert isinstance(cls, type), f"{cls} is not a class"
        df = _feasible_df()
        # Exercise DDFweak / NDDFweak construction directly (proves the
        # assert_DEAweak_with_indices / assert_NDDFweak_with_indices helpers
        # already in tools.py are sufficient -- no further migration needed).
        # sent "K L E=Y:CO2" -> 3 inputs, 1 output, 1 undesirable output.
        m1 = DDFweak(df, sent="K L E=Y:CO2", gy=[1], gx=[0, 0, 0], gb=[1])
        m2 = NDDFweak(df, sent="K L E=Y:CO2", gy=[1], gx=[0, 0, 0], gb=[1])
        assert m1 is not None and m2 is not None


# ===========================================================================
# Chapter 9 : deabook.DEAmb.MB
# ===========================================================================
from deabook.DEAmb import MB  # noqa: E402


@pytest.mark.skipif(SOLVER is None, reason="No local LP solver (highs/cbc) available")
class TestDEAmbMB:
    def test_mb_construct_and_optimize(self):
        rng = np.random.default_rng(3)
        mb = pd.DataFrame({
            "K": rng.uniform(10, 50, 8),
            "L": rng.uniform(20, 60, 8),
            "F": rng.uniform(5, 15, 8),     # polluted input
            "E": rng.uniform(40, 100, 8),   # desirable output (this chapter only)
            "CO2": rng.uniform(10, 40, 8),
        })
        # sent: "K L+F=E:CO2"  -> non-polluted inputs {K,L}, polluted {F},
        #                                  output {E}, undesirable {CO2}
        model = MB(mb, sent="K L+F=E:CO2", sx=[[0, 0, 1.0]], sy=[[1.0]], level=3)
        res, _info = model.optimize(solver=SOLVER)
        assert isinstance(res, pd.DataFrame) and len(res) == len(mb)


# ===========================================================================
# Chapter 10 : deabook.CNLSSDFDDFweak  (StoNED weak-disposability)
# ===========================================================================
from deabook.CNLSSDFDDFweak import CNLSSDweak, CNLSDDFweak  # noqa: E402


class TestCNLSWeak:
    def test_construct(self):
        df = _feasible_df()
        m1 = CNLSSDweak(df, sent="K L E=Y:CO2", gy=[0], gx=[0, 0, 0], gb=[1],
                        cet=CET_MULT, fun=FUN_PROD, rts=RTS_VRS1)
        m2 = CNLSDDFweak(df, sent="K L E=Y:CO2", gy=[0], gx=[0, 0, 0], gb=[1],
                         fun=FUN_PROD, rts=RTS_VRS1)
        assert m1 is not None and m2 is not None

    @pytest.mark.skipif(not HAS_CNLS_SOLVER,
                        reason="CNLS requires external solver mosek/knitro")
    def test_optimize(self):
        df = _feasible_df()
        m = CNLSDDFweak(df, sent="K L E=Y:CO2", gy=[0], gx=[0, 0, 0], gb=[1],
                        fun=FUN_PROD, rts=RTS_VRS1)
        # Availability is not enough to guarantee a MOSEK license on CI.
        selected = "mosek" if _solver_available("mosek") else "knitro"
        if selected == "mosek":
            import mosek
            try:
                with mosek.Env() as env:
                    env.checkoutlicense(mosek.feature.pton)
            except mosek.Error as exc:
                pytest.skip(f"MOSEK license unavailable: {exc.errno}")
        m.optimize(solver=selected)
        assert m.optimization_status == 1
        assert np.isfinite(m.get_residual()).all()



# ===========================================================================
# Chapter 10 : deabook.CQERDFDDFweak  (convex QR with weak disposability)
# ---------------------------------------------------------------------------
# CQRweak / CQRDDFweak import ``amplpy`` at module load time, so the whole
# group is skipped when that optional solver package is not installed. This is
# an *external* dependency, not a missing package module.
# ===========================================================================
@pytest.mark.skipif(  # pragma: no cover - depends on optional solver
    not HAS_AMPLPY,
    reason="CQERDFDDFweak requires the optional 'amplpy' solver package",
)
class TestCQRWeak:
    def test_cqrweak(self):
        from deabook.CQERDFDDFweak import CQRweak
        df = _feasible_df()
        m = CQRweak(df, sent="K L E=Y:CO2", tau=0.5, cet=CET_ADDI,
                    fun=FUN_PROD, rts=RTS_VRS)
        m.optimize()
        te = m.get_technical_efficiency()
        assert isinstance(te, (pd.DataFrame, pd.Series))

    def test_cqrddfweak(self):
        from deabook.CQERDFDDFweak import CQRDDFweak
        df = _feasible_df()
        m = CQRDDFweak(df, sent="K L E=Y:CO2", tau=0.5, gy=[0], gx=[0, 0, 0], gb=[1],
                       fun=FUN_PROD, rts=RTS_VRS)
        m.optimize()
        te = m.get_technical_efficiency()
        assert isinstance(te, (pd.DataFrame, pd.Series))


# ===========================================================================
# Chapter 3 (static panel wrapper): PANEL.tech
# ---------------------------------------------------------------------------
# Wraps any static efficiency factory (ccr/bcc/hyperbolic/ite/ddf/sbm) into a
# panel estimator that loops over periods under a reference technology (§3.4).
# ===========================================================================
from deabook import PANEL  # noqa: E402
from deabook import DEA, DEAweak, DEAmb
from deabook.constant import RTS_CRS, RTS_VRS1


@pytest.mark.skipif(SOLVER is None, reason="No local LP solver (highs/cbc) available")
class TestPanel:
    def test_ccr_all_techs(self):
        df = _feasible_df()
        for tech in (CONTEMPORARY, GLOBAL, SEQUENTIAL, WINDOW):
            res = PANEL.tech(df, "firm", "year", "K L E=Y", ccr,
                        orient="oo", tech=tech, solver=SOLVER).solve()
            assert isinstance(res, pd.DataFrame) and len(res) == len(df)
            assert res["firm"].notna().all() and res["year"].notna().all()
            assert res["te"].notna().all()

    def test_pivot_fallback_to_rho(self):
        df = _feasible_df()
        p = PANEL.tech(df, "firm", "year", "K L E=Y", ddf,
                  gy=[1], gx=[0, 0, 0], rts=RTS_CRS,
                  tech=CONTEMPORARY, solver=SOLVER)
        res = p.solve()
        wide = p.pivot()                                  # falls back to "rho"
        assert wide.shape == (4, 3)                       # 4 firms x 3 years
        assert (res["rho"].notna()).all()

    def test_global_frontier_is_larger(self):
        # GLOBAL reference set >= CONTEMPORARY reference set, so the global
        # efficiency score is weakly <= the contemporary score for every DMU.
        df = _feasible_df()
        g = PANEL.tech(df, "firm", "year", "K L E=Y", ccr,
                  orient="oo", tech=GLOBAL, solver=SOLVER).solve().set_index(["firm", "year"])
        c = PANEL.tech(df, "firm", "year", "K L E=Y", ccr,
                  orient="oo", tech=CONTEMPORARY, solver=SOLVER).solve().set_index(["firm", "year"])
        idx = c.index
        assert (g["te"].reindex(idx).values <= c["te"].reindex(idx).values + 1e-9).all()

    def test_panel_with_sddf(self):
        # DEA.sddf (Färe-Grosskopf 2010) and DEAweak.sddf (Pham-Zelenyuk 2018)
        # must plug into Panel like any other static factory.
        df = _feasible_df()
        p1 = PANEL.tech(df, "firm", "year", "K L E=Y", DEA.sddf,
                   rts=RTS_VRS1, tech=CONTEMPORARY, solver=SOLVER)
        r1 = p1.solve()
        assert isinstance(r1, pd.DataFrame) and len(r1) == len(df)
        assert "te" in r1.columns
        wide1 = p1.pivot()
        assert wide1.shape == (4, 3)

        r2 = PANEL.tech(df, "firm", "year", "K L E=Y:CO2", DEAweak.sddf,
                   rts=RTS_VRS1, tech=CONTEMPORARY, solver=SOLVER).solve()
        assert isinstance(r2, pd.DataFrame) and len(r2) == len(df)
        assert any(c.startswith("slack_b_") for c in r2.columns)

    def test_rts_ignored_by_ccr(self):
        # Passing rts= to ccr must not raise (kwarg filtering).
        df = _feasible_df()
        res = PANEL.tech(df, "firm", "year", "K L E=Y", ccr,
                    orient="oo", rts=RTS_VRS1, tech=CONTEMPORARY, solver=SOLVER).solve()
        assert len(res) == len(df)

    @staticmethod
    def _panel_with_bad(n=24, seed=11):
        """Panel with a polluting input F and undesirable output CO2."""
        rng = np.random.default_rng(seed)
        K = rng.uniform(10, 50, n)
        L = rng.uniform(20, 60, n)
        F = rng.uniform(5, 15, n)
        E = rng.uniform(40, 100, n)
        Y = 2 * K + 3 * L + F + 5 * E + rng.uniform(0, 20, n)
        CO2 = 0.3 * F + rng.uniform(0, 2, n)
        firms = np.repeat(range(1, 7), 4)
        years = np.tile([2016, 2017, 2018, 2019], 6)
        return pd.DataFrame(
            {"K": K, "L": L, "F": F, "E": E, "Y": Y, "CO2": CO2,
             "firm": firms, "year": years}
        )

    def test_weak_ddf_via_panel(self):
        # DEAweak.ddf (undesirable-output DDF) drives Panel and reports rho.
        df = self._panel_with_bad()
        p = PANEL.tech(df, "firm", "year", "K L F E=Y:CO2", DEAweak.ddf,
                  rts=RTS_CRS, tech=CONTEMPORARY, solver=SOLVER)
        res = p.solve()
        assert len(res) == len(df)
        assert res["firm"].notna().all() and res["year"].notna().all()
        assert res["rho"].notna().all()
        assert p.pivot(value="rho").shape == (6, 4)

    def test_amb_mb_via_panel(self):
        # DEAmb.MB returns a (DataFrame, info) tuple; Panel unwraps it to the df.
        df = self._panel_with_bad()
        p = PANEL.tech(df, "firm", "year", "K L+F=E:CO2", DEAmb.MB,
                  sx=[[0, 0, 1.0]], sy=[[1.0]], level=3,
                  tech=CONTEMPORARY, solver=SOLVER)
        res = p.solve()
        assert len(res) == len(df)
        assert res["obj"].notna().all()
        assert p.pivot(value="obj").shape == (6, 4)
