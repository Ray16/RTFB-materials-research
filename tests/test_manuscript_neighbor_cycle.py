from copy import deepcopy
import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('neighbor_paper_claims',ROOT/'manuscript/check_claims.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


@pytest.fixture
def data():
    folder=ROOT/'manuscript/evidence/2026-10-25_neighbor_cycle'
    names=['neighbor_cycle.json','pmdi_first_site_grafting.json','control_addendum.json',
        'control_archive_relocation.json','control_full_replay.json','control_example_manifest.json']
    return {name:json.loads((folder/name).read_text()) for name in names}


def test_complete_finite_cycle_and_extracted_control_replay_reconstruct(data):
    result=module.audit_neighbor_cycle(data.__getitem__)
    assert result['matched_states']==4
    assert result['deltaDeltaG_kJ_mol']==pytest.approx(-1.86161975)
    assert not result['conformer_convergence_inferred']
    assert module.audit_control_addendum(data.__getitem__)['displaced_evaluations']==384


@pytest.mark.parametrize('change',['missing_state','wrong_sites','G','conditional','convergence','nan'])
def test_changed_neighbor_claims_fail(data,change):
    data=deepcopy(data);row=data['neighbor_cycle.json']
    if change=='missing_state':row['states'].pop('site3')
    elif change=='wrong_sites':row['states']['pair']['task']['occupied_sites']=[2,4]
    elif change=='G':row['states']['pair']['G_eV']+=.1
    elif change=='conditional':row['by_cutoff']['50']['models']['free_ions']['graft_3_given_2_kJ_mol']+=1
    elif change=='convergence':row['conformer_convergence_established']=True
    else:row['states']['pair']['E_eV']=float('nan')
    with pytest.raises(ValueError):module.audit_neighbor_cycle(data.__getitem__)


@pytest.mark.parametrize('change',['counts','original_read','GPU','nan'])
def test_control_bundle_coverage_and_isolation_cannot_be_overstated(data,change):
    data=deepcopy(data)
    if change=='counts':data['control_archive_relocation.json']['displaced_evaluations']=392
    elif change=='original_read':data['control_archive_relocation.json']['isolated_replays']['full']['original_project_open_attempts']=['original record']
    elif change=='GPU':data['control_full_replay.json']['CUDA_or_xTB_loaded']=True
    else:data['control_full_replay.json']['states'][0]['maximum_force_error_Ha_Bohr']=float('nan')
    with pytest.raises(ValueError):module.audit_control_addendum(data.__getitem__)
