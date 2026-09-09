"""Regression checks for book-case indexing, outputs, and solver contracts."""
import importlib.util
from pathlib import Path
import shutil

import numpy as np
import pandas as pd
import pytest
from pyomo.environ import ConcreteModel, Constraint, NonNegativeReals, Objective, Var, SolverFactory

from deabook import data_path
from deabook import core
from deabook.example.book_cases import run
from deabook.example.book_cases.loaders import DATA_FILES, read_raw


def panel():
    return pd.DataFrame({'id': ['a', 'a'], 'year': [2018, 2019],
                         'x': [1., 1.], 'y': [1., 2.]})


@pytest.mark.parametrize('labels', [[0, 1], [1, 0], [10, 20], ['first', 'last']])
def test_dynamic_hand_calculation_independent_of_labels(labels):
    result = run.dynamic(panel().set_axis(labels), ['x'], ['y']).iloc[0]
    # Single CRS producer doubles output with unchanged input.
    assert result.MQ == pytest.approx(2.)
    assert result.EC == pytest.approx(1.)
    assert result.TC == pytest.approx(2.)


def test_dynamic_row_order_invariance():
    d = panel()
    expected = run.dynamic(d, ['x'], ['y'])
    pd.testing.assert_frame_equal(run.dynamic(d.iloc[::-1], ['x'], ['y']), expected)


def test_dynamic_duplicates_rejected():
    d = panel()
    with pytest.raises(ValueError, match='Duplicate id-year'):
        run.dynamic(pd.concat([d, d.iloc[:1]]), ['x'], ['y'])


def test_dynamic_unmatched_retained_as_missing():
    d = panel()
    d.loc[1, 'id'] = 'new'
    result = run.dynamic(d, ['x'], ['y'])
    assert len(result) == 1 and result.iloc[0]['id'] == 'new'
    assert result[['D00', 'D10', 'D01', 'D11', 'MQ', 'EC', 'TC']].isna().all().all()


def test_dynamic_energy_hand_calculation():
    d = pd.DataFrame({'id': ['a', 'a'], 'year': [2018, 2019],
                      'energy': [2., 1.], 'output': [1., 1.]})
    result = run.dynamic(d, ['energy'], ['output'], mode='energy').iloc[0]
    assert result.MQ == pytest.approx(2.)
    assert result.EC == pytest.approx(1.)
    assert result.TC == pytest.approx(2.)


def test_dynamic_caches_repeated_distances(monkeypatch):
    d = pd.DataFrame({'id': ['a'] * 3, 'year': [2017, 2018, 2019],
                      'x': [1., 1., 1.], 'y': [1., 2., 4.]})
    original = run.radial
    calls = []
    def counted(*args, **kwargs):
        calls.append(1)
        return original(*args, **kwargs)
    monkeypatch.setattr(run, 'radial', counted)
    result = run.dynamic(d, ['x'], ['y'])
    np.testing.assert_allclose(result.MQ, [2., 2.])
    assert len(calls) == 7  # Eight requests, one shared contemporary distance.


def test_installed_import_does_not_write(tmp_path, monkeypatch):
    source = Path(run.__file__)
    installed = tmp_path/'site-packages/deabook/example/book_cases/run.py'
    installed.parent.mkdir(parents=True)
    shutil.copy2(source, installed)
    monkeypatch.chdir(tmp_path)
    spec = importlib.util.spec_from_file_location('installed_case_run', installed)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert not (installed.parent/'results').exists()
    assert not (tmp_path/'book-case-output').exists()
    assert module.OUT == tmp_path/'book-case-output/results'
    module.configure_output(tmp_path/'custom')
    assert not (tmp_path/'custom').exists()


def test_custom_output_creates_parents_on_save(tmp_path, monkeypatch):
    monkeypatch.setattr(run, 'OUT', tmp_path/'nested/results')
    monkeypatch.setattr(run, 'IMG', tmp_path/'nested/images')
    run.FUN[1]()
    assert (run.OUT/'ch01_full.csv').is_file()
    assert (run.IMG/'case_ch01.png').is_file()
    assert (run.IMG/'case_ch01.pdf').is_file()


@pytest.mark.parametrize('name', DATA_FILES)
def test_official_data_locator(name):
    assert Path(data_path(DATA_FILES[name])).is_file()
    assert not read_raw(name).empty


def infeasible_model():
    m = ConcreteModel()
    m.x = Var(domain=NonNegativeReals)
    m.c = Constraint(expr=m.x <= -1)
    m.obj = Objective(expr=m.x)
    return m


@pytest.mark.skipif(not SolverFactory('highs').available(exception_flag=False), reason='HiGHS required')
def test_strict_solver_failure_and_recovery(monkeypatch):
    monkeypatch.setattr(core, '_SOLVER_OK', {})
    config = core.SolverConfig(solver='highs')
    with pytest.raises(RuntimeError, match='infeasible'):
        core.solve_pyomo_model(infeasible_model(), config, strict=True)
    outcome = core.solve_pyomo_model(infeasible_model(), config, strict=False)
    assert not outcome.ok and outcome.termination_condition == 'infeasible'
    # A failed model must not disable HiGHS for the next feasible model.
    m = infeasible_model()
    m.c.set_value(m.x >= 1)
    assert core.solve_pyomo_model(m, config, strict=True).ok
    assert m.x.value == pytest.approx(1.)
