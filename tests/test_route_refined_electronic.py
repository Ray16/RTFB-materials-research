"""Fresh geometry scoring cannot inherit an old energy or relax state/resources."""
import json
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
import shlex

import pytest

from redox.core.protocol import ACTIVE_SP
from redox.core.hashing import sha256_file
from redox.workflows import route_refined_electronic as api
from redox.workflows.screening_energy_panel import output_budget
from redox.workflows.molecular_native_namespace import load
from redox.core.common import ROOT


def plan(key='phenothiazine_neu_mode_minus'):
    q,m=api.STATES[key];budget=output_budget(976,m,False)
    return dict(role='screening_energy_control',energy_control=True,first_route_refined_electronic=True,
        route_electronic_profile=api.PROFILE,candidate_alias='phenothiazine',source_id=key,
        source_state=dict(id=key,smiles=api.SMILES,charge=q,mult=m),electronic_protocol=deepcopy(ACTIVE_SP),
        reference_energy_Ha=None,stored_energy_comparison_performed=False,atoms=39,nao=976,has_iron=False,
        maximum_host_RSS_GiB=40,minimum_free_GiB=32,minimum_GPU_memory_MiB=31744,
        maximum_peak_GPU_MiB=31744,maximum_wall_seconds=7200,minimum_live_storage_GiB=32,
        comparison_absolute_tolerance_Ha=1e-7,native_G_qualified=False,chemical_free_energy_assembled=False,
        redox_potential_qualified=False,material_admission=False,maximum_durable_output_bytes=budget,
        minimum_initial_storage_GiB=32+budget/2**30)


@pytest.mark.parametrize('key',list(api.STATES))
def test_frozen_known_dimensions_match_shared_resource_definition_and_exact_state(key):
    p=plan(key);assert api.validate_plan(p)['reference_energy_Ha'] is None
    assert p['maximum_durable_output_bytes']==api.OUTPUT_BYTES[p['source_state']['mult']]


@pytest.mark.parametrize('field,value',[('reference_energy_Ha',-1.),('stored_energy_comparison_performed',True),
    ('minimum_GPU_memory_MiB',15360),('electronic_protocol',{}),('nao',397),('native_G_qualified',True),
    ('maximum_durable_output_bytes',1),('minimum_initial_storage_GiB',32),('route_electronic_profile','other')])
def test_historical_reference_state_method_and_resource_substitution_are_refused(field,value):
    p=plan();p[field]=value
    with pytest.raises(ValueError):api.validate_plan(p)


def test_oxidized_state_cannot_silently_become_a_singlet():
    p=plan('phenothiazine_ox_stationarity');p['source_state']['mult']=1
    with pytest.raises(ValueError):api.validate_plan(p)


def test_fresh_input_requires_exact_accepted_geometry_and_route(tmp_path):
    p=plan();(tmp_path/'input.xyz').write_text('unchanged checked geometry');p['geometry_sha256']=sha256_file(tmp_path/'input.xyz')
    review=dict(status='route_primary_geometry_independently_checked',key=p['source_id'],smiles=api.SMILES,
        total_charge=0,multiplicity=1,geometry_sha256=p['geometry_sha256'],native_G_qualified=False,Hessian_class_qualified=False)
    (tmp_path/'geometry_review.json').write_text(json.dumps(review));(tmp_path/'route_gate.json').write_text(json.dumps(dict(proceed_to_full_redox=True,exact_route_binding_verified=True)))
    api.validate_inputs(p,tmp_path)
    (tmp_path/'input.xyz').write_text('old or different geometry')
    with pytest.raises(ValueError):api.validate_inputs(p,tmp_path)


@pytest.mark.parametrize('launcher',['redox.workflows.route_refined_electronic','redox.workflows.counterion_thermochemistry'])
def test_only_explicit_maintained_scientific_adapters_reach_the_remote_launcher(tmp_path,monkeypatch,launcher):
    module=load(ROOT/'scripts/fleet/run_screening_energy_panel.py','refined_adapter_test');calls=[]
    def run(argv,**kwargs):calls.append(argv);return SimpleNamespace(stdout='{"status":"recorded"}')
    monkeypatch.setattr(module.subprocess,'run',run)
    module.launch(dict(namespace=str(tmp_path/'namespace.json'),launcher_module=launcher),'host',0,'/python',tmp_path/'claim',1)
    assert shlex.split(calls[0][-1])[-1]==launcher


def test_scientific_child_stays_under_the_original_budget_observer_until_completion(tmp_path,monkeypatch):
    from redox.workflows import molecular_native_namespace as namespace_api,screening_control
    case=tmp_path/'case';case.mkdir();(case/'plan.json').write_text('{}')
    spec=dict(case=str(case),sources=str(tmp_path/'sources'));child=object();record=dict(status='submitted');calls=[]
    monkeypatch.setattr(namespace_api,'verify',lambda p:spec)
    monkeypatch.setattr(api,'validate_inputs',lambda *args:None)
    monkeypatch.setattr(screening_control,'launch',lambda *args:(child,record))
    def observe(process,actual_spec,launch_record,output):
        assert process is child and actual_spec is spec and launch_record is record
        calls.append(output);return dict(status='closed',returncode=0)
    monkeypatch.setattr(screening_control,'observe_budget',observe)
    result=api.launch(tmp_path/'namespace.json',0,tmp_path/'placement')
    assert result==dict(status='closed',returncode=0) and calls==[tmp_path/'placement']


def test_near_degenerate_guess_cannot_replace_the_actual_selected_spin():
    reviewer=load(ROOT/'scripts/validation/review_screening_energy.py','selected_refined_spin_regression')
    near=dict(energy_Ha=-100.+4.4e-10,spin_contamination=.0177544,checkpoint='near')
    exact=dict(energy_Ha=-100.,spin_contamination=.0177476,checkpoint='selected')
    result=dict(e_smd_Ha=-100.,smd_spin_contam=.0177)
    assert reviewer.selected_spin_record([near,exact],result)['checkpoint']=='selected'
    result['smd_spin_contam']=.0178
    with pytest.raises(ValueError):reviewer.selected_spin_record([near,exact],result)
