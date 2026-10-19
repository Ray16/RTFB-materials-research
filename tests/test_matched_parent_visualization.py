"""The plotted 12-potential benchmark excludes spacings and preserves step meaning."""
import importlib.util
from pathlib import Path
import pytest


def module():
    p=Path(__file__).resolve().parents[1]/'scripts/plotting/validation/plot_matched_parent_validation.py'
    spec=importlib.util.spec_from_file_location('parent_visualization',p)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def test_first_and_second_step_depend_on_the_chemical_ladder():
    m=module()
    assert m.step(dict(family='phenothiazine (p-type)',event='ox1->neu'))==1
    assert m.step(dict(family="viologen (4,4'-bipyridinium)",event='ox1->neu'))==2
    assert m.step(dict(family='aromatic imide',event='neu->red1'))==1
    assert m.step(dict(family='aromatic imide',event='red1->red2'))==2


def test_exact_parent_snapshot_reproduces_the_reported_benchmark():
    m=module();root=Path(__file__).resolve().parents[1]
    rows,stats,_=m.load_snapshot(root/'results/reviews/2026-10-19/full_matched_parent_final')
    assert len(rows)==12 and len({r['id'] for r in rows})==8
    assert sum(r['step']==1 for r in rows)==8 and sum(r['step']==2 for r in rows)==4
    assert stats['MAE']==pytest.approx(.2090416666666667)
    assert stats['RMSE']==pytest.approx(.24571566392885905)
    outlier=next(r for r in rows if r['label']=='DMeOEPT' and r['step']==2)
    assert outlier['exp_V']==.655 and outlier['calc_V']==1.202
    assert outlier['error_V']==pytest.approx(.547)
    assert all(r['kind']=='E' and r['tier']=='A' for r in rows)
