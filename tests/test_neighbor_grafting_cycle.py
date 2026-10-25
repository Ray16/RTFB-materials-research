import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('neighbor_cycle_test',ROOT/'scripts/validation/analyze_neighbor_grafting.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_neighbor_penalty_cancels_shared_base_and_both_attachment_paths_close():
    values=dict(empty=-100.,site2=-110.,site3=-111.,pair=-122.)
    a,b=[module.cycle(values,constant) for constant in (-50.,75.)]
    assert a['deltaDeltaG_kJ_mol']==pytest.approx(-module.EV_KJ)
    assert b['deltaDeltaG_kJ_mol']==pytest.approx(a['deltaDeltaG_kJ_mol'])
    assert a['direct_double_graft_kJ_mol']==pytest.approx(a['graft_2_on_empty_kJ_mol']+a['graft_3_given_2_kJ_mol'])
    assert a['direct_double_graft_kJ_mol']==pytest.approx(a['graft_3_on_empty_kJ_mol']+a['graft_2_given_3_kJ_mol'])
    assert b['graft_3_given_2_kJ_mol']-a['graft_3_given_2_kJ_mol']==pytest.approx(125*module.EV_KJ)


@pytest.mark.parametrize('bad',[{'empty':0.,'site2':1.,'site3':2.},
    {'empty':0.,'site2':1.,'site3':2.,'pair':float('nan')},
    {'empty':0.,'site2':1.,'site3':2.,'pair':float('inf')}])
def test_missing_or_nonfinite_native_state_never_becomes_a_neighbor_result(bad):
    with pytest.raises(ValueError):module.cycle(bad,0.)


def test_favorable_absolute_grafting_does_not_require_favorable_nonadditivity():
    result=module.cycle(dict(empty=0.,site2=-1.,site3=-2.,pair=-2.5),-10.)
    assert result['graft_2_given_3_kJ_mol']<0
    assert result['deltaDeltaG_kJ_mol']>0
