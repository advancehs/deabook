"""Shared OOP helpers for dynamic productivity-index modules."""

from __future__ import annotations

import numpy as np
import pandas as pd

from .constant import CONTEMPORARY, TOTAL
from .utils import tools


class DynamicProductivityMixin:
    """Small shared layer for Malmquist/Luenberger wrappers.

    The existing public classes still contain the model-specific schedules and
    formulas, but construction, orientation naming, year validation and the
    ``optimize`` method live here so future dynamic DEA variants share one
    lifecycle.
    """

    def _init_dynamic_productivity(self, data, id, year, tech, email, solver):
        if year not in data.columns:
            raise ValueError(f"Year column '{year}' not found in data.")
        if id not in data.columns:
            raise ValueError(f"ID column '{id}' not found in data.")
        self.id = id
        self.year_name = year
        self.tlt = pd.Series(data[year]).drop_duplicates().sort_values().reset_index(drop=True)
        if self.tlt.empty:
            raise ValueError("The panel contains no time periods.")
        self.tech = tech
        self.email = email
        self.solver = solver
        self.datazz = data.copy()

    def _init_dynamic_dea(self, data, id, year, sent, gy, gx, rts, tech, email, solver):
        """Initialize shared state for dynamic DEA/DDF productivity classes."""
        self._init_dynamic_productivity(data, id, year, tech, email, solver)
        self.gy, self.gx, self.inputvars, self.outputvars = tools.assert_MQDEA(data, sent, gy, gx)
        self.xcol = list(self.inputvars)
        self.ycol = list(self.outputvars)
        self.rts = rts
        self._set_dynamic_orientation(self.gx, self.gy)
        self._dispatch_dynamic_technology(data, sent, id, year)

    def _init_dynamic_weak_dea(
        self,
        data,
        id,
        year,
        sent,
        gy,
        gx,
        gb,
        rts,
        tech,
        email,
        solver,
        dynamic=None,
    ):
        """Initialize shared state for dynamic weak-disposability classes."""
        self._init_dynamic_productivity(data, id, year, tech, email, solver)
        self.gy, self.gx, self.gb, self.inputvars, self.outputvars, self.unoutputvars = tools.assert_MQDEAweak(
            data, sent, gy, gx, gb
        )
        self.xcol = list(self.inputvars)
        self.ycol = list(self.outputvars)
        self.bcol = list(self.unoutputvars)
        self.rts = rts
        if dynamic is not None:
            self.dynamic = dynamic
        self._set_dynamic_orientation(self.gx, self.gy, self.gb)
        self._dispatch_dynamic_technology(data, sent, id, year)

    def _set_dynamic_orientation(self, gx, gy, gb=None):
        sum_gx = sum(gx)
        sum_gy = sum(gy)
        sum_gb = sum(gb or [])
        self.input_oriented = sum_gx >= 1 and sum_gy == 0 and sum_gb == 0
        self.output_oriented = sum_gy >= 1 and sum_gx == 0 and sum_gb == 0
        self.undesirable_oriented = sum_gb >= 1 and sum_gx == 0 and sum_gy == 0
        self.unoutput_oriented = self.undesirable_oriented
        self.hyper_oriented = sum_gx >= 1 and sum_gy >= 1 and sum_gb == 0
        self.hyper_orientedyx = self.hyper_oriented
        self.hyper_orientedyb = sum_gb >= 1 and sum_gy >= 1 and sum_gx == 0
        self.hyper_orientedxb = sum_gx >= 1 and sum_gb >= 1 and sum_gy == 0
        self.hyper_orientedyxb = sum_gx >= 1 and sum_gy >= 1 and sum_gb >= 1

    def _dispatch_dynamic_technology(self, data, sent, id, year):
        if self.tech == TOTAL:
            print("Calculating D11 (Total frontier) for all periods...")
            self.get_total(data, sent, id, year)
            print("TOTAL tech calculation finished.")
            return
        if self.tech == CONTEMPORARY:
            print("Calculating CONTEMPORARY tech components (D11, D12, D21)...")
            self.get_contemp(data, sent, id, year)
            return
        raise ValueError(f"Unsupported technology type '{self.tech}'. Must be '{TOTAL}' or '{CONTEMPORARY}'.")

    def optimize(self):
        if not hasattr(self, "datazz"):
            raise RuntimeError("Dynamic productivity calculation failed during initialization.")
        return self.datazz

    def _year_selector(self, year_name, value):
        return f"{year_name}=[{value}]"

    def _all_year_selector(self, year_name):
        return f"{year_name}=[{','.join(map(str, self.tlt.tolist()))}]"

    def _join_component(self, frame):
        self.datazz = self.datazz.join(frame, how="left")

    def _ratio_to_previous(self, source_col, result_col, group_col):
        if source_col not in self.datazz.columns or self.datazz[source_col].isnull().all():
            self.datazz[result_col] = np.nan
            return
        self.datazz[source_col] = pd.to_numeric(self.datazz[source_col], errors="coerce")
        previous_col = f"{source_col}_prev"
        self.datazz[previous_col] = self.datazz.groupby(group_col)[source_col].transform(lambda x: x.shift(1))
        self.datazz[result_col] = self.datazz[source_col] / self.datazz[previous_col]
        self.datazz.drop(columns=[previous_col], inplace=True)
