import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

PATH=Path(__file__).parents[1]/'scripts/pipeline/smooth_cds_control.py'
spec=importlib.util.spec_from_file_location('smooth_CDS_control_test',PATH)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_energy_only_request_does_not_evaluate_gradient():
    calls=[];stats=dict(energy_calls=0,gradient_calls=0)
    item=module.LazyCDS(None,lambda _:calls.append('energy') or .2,
        lambda _:calls.append('gradient') or np.zeros((2,3)),stats)
    assert item[0]==item[0]==.2
    assert calls==['energy']
    assert np.array_equal(item[1],np.zeros((2,3)))
    assert calls==['energy','gradient']
    assert stats==dict(energy_calls=1,gradient_calls=1)


def test_nonfinite_smooth_energy_is_rejected():
    item=module.LazyCDS(None,lambda _:float('nan'),lambda _:None,dict(energy_calls=0,gradient_calls=0))
    with pytest.raises(ValueError,match='Nonfinite'):item[0]


@pytest.fixture
def study(tmp_path):
    roles=['contact','empty','capped','site2','site3','precursor','local','chloride','base','protonated','pair']
    row=dict(status='experimental_smooth_CDS_component_study_complete',source_sha256={},
        states={r:dict(finest_grid_change_kJ_mol=.01) for r in roles})
    row['states']['contact'].update(sampled_gradient=dict(tested_components=4,maximum_error_Ha_Bohr=1e-9),
        CPU_legacy_contact_reference_error_Ha=0.,rigid_rotation_spread_kJ_mol=.01,
        contact_scan=dict(centered_curvature_Ha_Bohr2={'a':.002,'b':.00200001,'c':.00200001}))
    (tmp_path/'independent_review.json').write_text(json.dumps(dict(status='smooth_CDS_fixed_geometry_component_review_passed',source_sha256={})))
    def save(): (tmp_path/'summary.json').write_text(json.dumps(row))
    save()
    return tmp_path,row,save


def test_valid_component_evidence_keeps_qualification_scope(study):
    root,_,_=study
    assert str(root/'summary.json') in module.component_evidence(root)
    assert module.MODEL['primary_candidate_admission'] is False


@pytest.mark.parametrize('problem',['nan','missing','large_gradient','unstable_curvature'])
def test_bad_component_evidence_cannot_authorize_quantum_variant(study,problem):
    root,row,save=study
    if problem=='missing':del row['states']['site3']
    elif problem=='nan':row['states']['contact']['sampled_gradient']['maximum_error_Ha_Bohr']=float('nan')
    elif problem=='large_gradient':row['states']['contact']['sampled_gradient']['maximum_error_Ha_Bohr']=.01
    else:row['states']['contact']['contact_scan']['centered_curvature_Ha_Bohr2']['c']=1.
    save()
    with pytest.raises(ValueError):module.component_evidence(root)
