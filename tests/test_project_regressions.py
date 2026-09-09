import contextlib
import io
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from pyomo.environ import ConcreteModel, Var, Constraint, Objective, NonNegativeReals
from deabook import DEA
from deabook.constant import GLOBAL, CONTEMPORARY, WINDOW, OPT_LOCAL, CET_ADDI
from deabook.mqcore import reference_years, normalize_window
from deabook.utils import tools


@pytest.mark.parametrize('orient', ['oo', 'io'])
def test_cross_period_radial_productivity(orient):
    d = pd.DataFrame({'id':['a','a'],'year':[2018,2019], 'x':[1.,1.], 'y':[1.,2.]})
    r = DEA.mq(d, 'id', 'year', 'x=y', orient=orient).optimize(solver='highs')
    assert r.loc[1, 'MQ'] == pytest.approx(2.)
    assert r.loc[1, 'MEFFCH'] == pytest.approx(1.)
    assert r.loc[1, 'MTECHCH'] == pytest.approx(2.)


def test_duplicate_index_rejected_without_losing_rows():
    with pytest.raises(ValueError, match='index must be unique'):
        DEA.ccr(pd.DataFrame({'x':[1.,2.],'y':[1.,1.]},index=[0,0]), 'x=y')


@pytest.mark.parametrize('tech', [GLOBAL, CONTEMPORARY])
def test_missing_year_not_bridged(tech):
    d=pd.DataFrame({'id':['a','b','b','a','b'],'year':[2016,2016,2017,2018,2018],
                    'x':[1.]*5,'y':[1.,1.,2.,4.,4.]})
    r=DEA.mq(d,'id','year','x=y',tech=tech).optimize(solver='highs')
    col='mqpi' if tech==GLOBAL else 'MQ'
    assert pd.isna(r.loc[3,col])
    assert pd.isna(r.loc[3,'MEFFCH']) and pd.isna(r.loc[3,'MTECHCH'])
    assert r.loc[4,col] == pytest.approx(2.)


@pytest.mark.parametrize('width', [1,2,3,4,5])
def test_window_exact_width(width):
    years=list(range(11))
    assert len(reference_years(years,5,WINDOW,width)) == width
    if width==2: assert reference_years(years,5,WINDOW,width)==[5,6]


@pytest.mark.parametrize('width', [0,-1,2.5,'3',True])
def test_window_invalid_values(width):
    with pytest.raises(ValueError):normalize_window(width)


def test_legacy_solver_does_not_report_infeasible_as_success():
    m=ConcreteModel();m.x=Var(domain=NonNegativeReals)
    m.c=Constraint(expr=m.x<=-1);m.obj=Objective(expr=m.x)
    with pytest.raises(RuntimeError,match='infeasible'):
        tools.optimize_model(m,OPT_LOCAL,CET_ADDI,'highs')
    m.c.set_value(m.x>=1)
    status,ok=tools.optimize_model(m,OPT_LOCAL,CET_ADDI,'highs')
    assert ok==1 and m.x.value==pytest.approx(1.)


def test_quickstart_notebook_executes():
    root=Path(__file__).resolve().parents[1]
    notebook=json.loads((root/'docs/examples/intro.ipynb').read_text(encoding='utf-8'))
    scope={}
    with contextlib.redirect_stdout(io.StringIO()):
        for cell in notebook['cells']:
            if cell['cell_type']=='code':exec(''.join(cell['source']),scope)
    assert scope['wide'].shape == (6,2)
    assert np.isfinite(scope['wide'].to_numpy()).all()

@pytest.mark.parametrize('gx,gy', [([0],[1]),([1],[0])])
def test_ddf_cross_period_agrees_with_radial(gx,gy):
    d=pd.DataFrame({'id':['a','a'],'year':[2018,2019],'x':[1.,1.],'y':[1.,2.]})
    r=DEA.mqddf(d,'id','year','x=y',gx=gx,gy=gy).optimize(solver='highs')
    assert r.loc[1,'MQ']==pytest.approx(2.)


def test_rd_uses_vrs_cross_distances_and_retains_scale_factor():
    from deabook.mqcore import rd_decompose
    crs=pd.DataFrame({'id':['a','a'],'D11':[.5,.75],'D12':[np.nan,.9],'D21':[.4,np.nan]})
    vrs=pd.DataFrame({'id':['a','a'],'D11':[.8,.9],'D12':[np.nan,.95],'D21':[.7,np.nan]})
    out=rd_decompose(vrs,crs,'id').iloc[1]
    assert out.TECHCH_RD==pytest.approx(np.sqrt((.95/.9)*(.8/.7)))
    assert out.MQ==pytest.approx(out.EFFCH*out.TECHCH_RD*out.SCALECH_RD)
    assert out.SCALECH_RD != pytest.approx(1.)
