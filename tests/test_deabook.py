#!/usr/bin/env python

"""Tests for `deabook` package and its OOP DEA core."""

import unittest

import numpy as np
import pandas as pd

import deabook as deabook_pkg
from deabook import deabook
from deabook.DEA import DEA2, DDF2
from deabook.DEAweak import DEAweak2, DDFweak2
from deabook.DDFDUAL import DDFDUAL as StandaloneDDFDUAL
from deabook.DRFDUAL import DRFDUAL
from deabook.constant import LEFT, RTS_CRS, RTS_VRS, RTS_VRS1
from deabook.core import parse_formula, parse_sent, prepare_dea_data


SOLVER = "glpk"


class TestDeabook(unittest.TestCase):
    """Smoke and regression tests for documented public APIs."""

    def setUp(self):
        self.data = pd.DataFrame(
            {
                "K": [1.0, 2.0, 3.0, 4.0],
                "L": [1.0, 1.0, 2.0, 3.0],
                "Y": [1.0, 2.0, 2.5, 3.0],
                "CO2": [1.0, 1.5, 2.0, 2.5],
                "Year": [2018, 2018, 2019, 2019],
            }
        )

    def test_import_deabook_module(self):
        self.assertTrue(hasattr(deabook, "DEAModel"))
        self.assertTrue(hasattr(deabook_pkg, "DEAModel"))
        self.assertTrue(hasattr(deabook_pkg, "simple_dea_solve"))
        self.assertTrue(hasattr(deabook_pkg, "validate_dea_inputs"))

    def test_parse_formula_styles(self):
        formula = parse_formula("Y:CO2 = K L")
        self.assertEqual(formula.x, ["K", "L"])
        self.assertEqual(formula.y, ["Y"])
        self.assertEqual(formula.b, ["CO2"])
        sent = parse_sent("K L=Y:CO2")
        self.assertEqual(sent.x, ["K", "L"])
        self.assertEqual(sent.y, ["Y"])
        self.assertEqual(sent.b, ["CO2"])

    def test_prepare_dea_data_keeps_indices_and_filters(self):
        prepared = prepare_dea_data(
            self.data,
            "K L=Y",
            eval_query="Year=[2019]",
            ref_query="Year==2018",
        )
        self.assertEqual(prepared.evaluated_index, [2, 3])
        self.assertEqual(prepared.reference_index, [0, 1])
        self.assertEqual(prepared.x.shape, (2, 2))
        self.assertEqual(prepared.yref.shape, (2, 1))

    def test_prepare_dea_data_rejects_missing_and_nonpositive_values(self):
        with self.assertRaises(KeyError):
            prepare_dea_data(self.data, "K Missing=Y")
        bad = self.data.copy()
        bad.loc[0, "K"] = 0.0
        with self.assertRaises(ValueError):
            prepare_dea_data(bad, "K L=Y")
        bad = self.data.copy()
        bad.loc[0, "CO2"] = -1.0
        with self.assertRaises(ValueError):
            prepare_dea_data(bad, "K L=Y:CO2")

    def test_teaching_dea_model_output_and_input_orientation(self):
        output_result = deabook.DEAModel("Y = K L", self.data, orient="oo", rts="crs").solve(
            solver=SOLVER
        )
        self.assertIn("theta", output_result.columns)
        self.assertIn("te", output_result.columns)
        self.assertTrue(((output_result["te"] > 0) & (output_result["te"] <= 1)).all())

        input_result = deabook.simple_dea_solve(
            self.data,
            ["K", "L"],
            ["Y"],
            orient="io",
            rts="crs",
            solver=SOLVER,
        )
        self.assertTrue(np.allclose(input_result["theta"], input_result["te"]))

    def test_dea2_documented_columns(self):
        input_model = DEA2(self.data, "K L=Y", gy=[0], gx=[1, 0], rts=RTS_CRS)
        input_result = input_model.optimize(solver=SOLVER)
        self.assertIn("te", input_result.columns)
        self.assertTrue(((input_result["te"] > 0) & (input_result["te"] <= 1)).all())

        output_model = DEA2(self.data, "K L=Y", gy=[1], gx=[0, 0], rts=RTS_CRS)
        output_result = output_model.optimize(solver=SOLVER)
        self.assertIn("te", output_result.columns)
        self.assertTrue(((output_result["te"] > 0) & (output_result["te"] <= 1)).all())

        hyper_model = DEA2(self.data, "K L=Y", gy=[1], gx=[1, 0], rts=RTS_VRS1)
        hyper_result = hyper_model.optimize(solver=SOLVER)
        self.assertIn("tei", hyper_result.columns)
        self.assertIn("teo", hyper_result.columns)

    def test_ddf2_documented_columns(self):
        input_model = DDF2(self.data, "K L=Y", gy=[0], gx=[1, 0], rts=RTS_CRS)
        input_result = input_model.optimize(solver=SOLVER)
        self.assertIn("tei", input_result.columns)
        self.assertIn("objective_value", input_result.columns)

        output_model = DDF2(self.data, "K L=Y", gy=[1], gx=[0, 0], rts=RTS_VRS1)
        output_result = output_model.optimize(solver=SOLVER)
        self.assertIn("teo", output_result.columns)
        self.assertIn("objective_value", output_result.columns)

    def test_weak_dea_models_use_shared_oop_result_contract(self):
        bad_model = DEAweak2(
            self.data,
            "K L=Y:CO2",
            gy=[0],
            gx=[0, 0],
            gb=[1],
            rts=RTS_CRS,
        )
        bad_result = bad_model.optimize(solver=SOLVER)
        self.assertIn("te", bad_result.columns)
        self.assertIn("rho", bad_result.columns)
        self.assertIsInstance(bad_model.get_lamda(), pd.DataFrame)

        ddf_model = DDFweak2(
            self.data,
            "K L=Y:CO2",
            gy=[0],
            gx=[0, 0],
            gb=[1],
            rts=RTS_CRS,
        )
        ddf_result = ddf_model.optimize(solver=SOLVER)
        self.assertIn("teuo", ddf_result.columns)
        self.assertIn("objective_value", ddf_result.columns)

    def test_standalone_dual_models_share_shadow_price_contract(self):
        dual = StandaloneDDFDUAL(
            self.data,
            sent="K L=Y:CO2",
            gy=[1],
            gx=[1, 0],
            gb=[1],
            rts=RTS_CRS,
        )
        result = dual.optimize(solver=SOLVER)
        self.assertIn("obj", result.columns)
        self.assertIn("Input0's shadow price", result.columns)
        self.assertIn("Output0's shadow price", result.columns)
        self.assertIn("Undesirable Output0's shadow price", result.columns)

        drf = DRFDUAL(
            self.data,
            year=self.data["Year"],
            sent="K L=Y:CO2",
            fenmu="CO2",
            fenzi="K",
            side=LEFT,
            rts=RTS_CRS,
        )
        drf_result = drf.optimize(solver=SOLVER)
        self.assertEqual(list(result.columns), list(drf_result.columns))

    def test_legacy_rts_alias_is_available(self):
        self.assertEqual(RTS_VRS, RTS_VRS1)

    def test_direction_length_error_is_clear(self):
        with self.assertRaises(ValueError):
            DEA2(self.data, "K L=Y", gy=[0], gx=[1], rts=RTS_CRS)


if __name__ == "__main__":
    unittest.main()
