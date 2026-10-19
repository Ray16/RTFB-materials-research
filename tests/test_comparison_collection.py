"""The maintained comparison must retain pending candidates without promoting values."""
import importlib.util
from pathlib import Path

import pytest


def module():
    path=Path(__file__).resolve().parents[1]/'scripts/plotting/candidates/update_comparison_collection.py'
    spec=importlib.util.spec_from_file_location('candidate_collection_figures',path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result)
    return result


def row(gid, **kwargs):
    return dict(id=gid,family='imide',E1_vs_Fc_V=-1.,E2_vs_Fc_V=-2.,
                capacity_nominal_mAh_g=150.,accessible_capacity_mAh_g=150.,**kwargs)


def test_pending_and_new_candidates_remain_visible_without_unverified_energies():
    m=module();old=row('old');new=row('new')
    result=m.comparison_rows(dict(original39=[old]),[old,new],[old])
    assert {r['id'] for r in result}=={'old','new'}
    pending=next(r for r in result if r['id']=='new')
    assert pending['cohort']=='Additional' and pending['E1_vs_Fc_V'] is None
    assert pending['accessible_capacity_mAh_g'] is None
    assert pending['capacity_nominal_mAh_g']==150.


def test_model_hold_keeps_calculated_potentials_but_withholds_window_capacity():
    m=module();held=row('held',model_primary_admissible=False)
    result=m.comparison_rows(dict(original39=[held]),[held],[held])[0]
    assert result['E1_vs_Fc_V']==-1. and result['E2_vs_Fc_V']==-2.
    assert result['accessible_capacity_mAh_g'] is None
    assert result['capacity_nominal_mAh_g']==150.


def test_duplicate_cohort_ids_cannot_double_count_a_candidate():
    m=module();old=row('old')
    with pytest.raises(ValueError,match='unique'):
        m.comparison_rows(dict(original39=[old]),[old,old],[old])


def test_a_new_candidate_changes_the_collection_without_a_hardcoded_roster():
    m=module();old=row('old');new=row('registered_after_build')
    result=m.comparison_rows(dict(original39=[old]),[old,new],[old,new])
    assert len(result)==2 and all(r['E1_vs_Fc_V']==-1. for r in result)
