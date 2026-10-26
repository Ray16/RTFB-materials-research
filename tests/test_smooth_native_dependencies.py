"""Reject model/resource/parent substitutions before native fan-out."""
from copy import deepcopy
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).parents[1]
spec = importlib.util.spec_from_file_location('native_dependencies_test', ROOT/'scripts/pipeline/prepare_smooth_native_batch.py')
module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)


@pytest.fixture
def state():
    task = dict(smiles='C[C@H](Cl)F', atoms=8, active_basis_functions=128, charge=0, mult=1)
    point = dict(tasks=[task], experimental_method_id='smooth', experimental_model={'surface':'smooth'},
        quantum_settings={'basis':'svp'}, runtime_profile='runtime.json', maximum_host_RSS_GiB=80,
        maximum_peak_GPU_MiB=31744)
    review = dict(status='alternate_native_state_independently_reconstructed', experimental_method_id='smooth',
        primary_candidate_admission=False, conformer_convergence_established=False,
        peak_host_RSS_bytes=10*2**30, observed_peak_GPU_MiB=20000, G_eV=-100.,
        G_by_cutoff_eV={'25':-100.1,'50':-100.,'100':-99.9}, minimum_frequency_cm=5.,
        maximum_sampled_column_error=1e-5, frequency_reconstruction_error_cm=1e-8,
        nuclear_reconstruction_errors_eV={'25':0.,'50':0.,'100':0.})
    return deepcopy(point), point, review


def test_native_resources_may_transfer_between_configurations_without_transferring_G(state):
    plan, point, review = state
    plan['tasks'][0]['smiles']='C[C@@H](Cl)F'
    before=deepcopy(plan)
    module.resource_compatible(plan,point,review)
    assert plan==before and 'G_eV' not in plan


@pytest.mark.parametrize('change',['constitution','basis','charge','method','surface','runtime','settings',
                                 'promotion','unreviewed','RAM','GPU','nan_RAM','different_bound'])
def test_native_resource_context_mismatches_fail_closed(state,change):
    plan,point,review=state
    if change=='constitution':plan['tasks'][0]['smiles']='FCCCl'
    elif change=='basis':plan['tasks'][0]['active_basis_functions']+=1
    elif change=='charge':plan['tasks'][0]['charge']=-1
    elif change=='method':plan['experimental_method_id']='legacy'
    elif change=='surface':plan['experimental_model']['surface']='legacy'
    elif change=='runtime':plan['runtime_profile']='other.json'
    elif change=='settings':plan['quantum_settings']['basis']='svpd'
    elif change=='promotion':review['primary_candidate_admission']=True
    elif change=='unreviewed':review['status']='finished'
    elif change=='RAM':review['peak_host_RSS_bytes']=81*2**30
    elif change=='GPU':review['observed_peak_GPU_MiB']=32000
    elif change=='nan_RAM':review['peak_host_RSS_bytes']=float('nan')
    else:point['maximum_host_RSS_GiB']=256
    with pytest.raises(ValueError):module.resource_compatible(plan,point,review)


@pytest.mark.parametrize('change',['nan_G','infinite_G','cutoff','closure','mode','derivative','negative_error','nuclear'])
def test_native_review_requires_finite_reconstructed_G_and_physics(state,change):
    _,_,review=state
    if change=='nan_G':review['G_by_cutoff_eV']['50']=float('nan')
    elif change=='infinite_G':review['G_eV']=float('inf')
    elif change=='cutoff':review['G_by_cutoff_eV'].pop('25')
    elif change=='closure':review['G_eV']+=1
    elif change=='mode':review['minimum_frequency_cm']=-6
    elif change=='derivative':review['maximum_sampled_column_error']=.02
    elif change=='negative_error':review['frequency_reconstruction_error_cm']=-1.
    else:review['nuclear_reconstruction_errors_eV']['50']=1e-4
    with pytest.raises(ValueError):module.validate_review(review,'smooth')


def test_pending_or_failed_electronic_parent_never_prepares_native(tmp_path,monkeypatch):
    geometry=tmp_path/'geometry';geometry.mkdir();(geometry/'plan.json').write_text('{}')
    case=tmp_path/'energy';case.mkdir()
    entry=dict(campaign=str(tmp_path/'native'),geometry_campaign=str(geometry),
        geometry_plan_sha256=module.api.digest(geometry/'plan.json'),existing_checked_reference=False,
        electronic_options=[dict(campaign=str(case),review=str(case/'review.json'))])
    monkeypatch.setattr(module.api,'load_module',lambda *args:SimpleNamespace(owner_alive=lambda owner:False))
    assert module.prepare_ready(entry,{})=='upstream_pending'
    (case/'experimental_energy').mkdir()
    (case/'experimental_energy/status.json').write_text(json.dumps(dict(status='failed',finished=1,owner={})))
    assert module.prepare_ready(entry,{})=='upstream_failed_preserved'
    assert not (tmp_path/'native').exists()


def test_unreviewed_success_or_wrong_geometry_cannot_enter_native(tmp_path,monkeypatch):
    geometry=tmp_path/'geometry';geometry.mkdir();(geometry/'plan.json').write_text('{}')
    case=tmp_path/'energy';(case/'experimental_energy').mkdir(parents=True)
    (case/'experimental_energy/status.json').write_text(json.dumps(dict(
        status='experimental_hybrid_energy_finished_independent_review_pending',finished=1,owner={})))
    review=case/'review.json'
    entry=dict(campaign=str(tmp_path/'native'),geometry_campaign=str(geometry),
        geometry_plan_sha256=module.api.digest(geometry/'plan.json'),existing_checked_reference=False,
        electronic_options=[dict(campaign=str(case),review=str(review))])
    monkeypatch.setattr(module.api,'load_module',lambda *args:SimpleNamespace(owner_alive=lambda owner:False))
    assert module.prepare_ready(entry,{})=='upstream_pending'
    review.write_text('{}')
    monkeypatch.setattr(module.api,'read_plan',lambda *args:dict(source_geometry_campaign='another_geometry'))
    with pytest.raises(ValueError,match='different geometry'):module.prepare_ready(entry,{})
    assert not (tmp_path/'native').exists()


def test_native_placement_template_builds_without_GPU_imports():
    spec=importlib.util.spec_from_file_location('native_placement_test',ROOT/'scripts/pipeline/coordinate_smooth_native_batch.py')
    coordinator=importlib.util.module_from_spec(spec);spec.loader.exec_module(coordinator)
    placement=coordinator.native_placement({'library_path':'isolated'})
    assert callable(placement.launch_local)
