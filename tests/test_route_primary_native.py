"""The first-native pilot cannot silently change the molecule, state or coverage."""
import json
import pytest
import numpy as np

from redox.core.common import ROOT
from redox.workflows.route_primary_native import validate_plan,REFINED_PROFILE,CDS_PROFILE,CDS_POLICY
from redox.workflows.molecular_native_namespace import load


@pytest.fixture
def plan():
    return json.loads((ROOT/'validation/phenothiazine_primary_2026-11-02/neutral_native_plan.json').read_text())


def test_declared_first_native_profile_requires_sulfur_columns_and_keeps_G_pending(plan):
    assert validate_plan(plan)['native_nao']==397
    assert plan['extra_columns']==[[15,0],[15,1],[15,2]] and not plan['native_G_qualified']


@pytest.mark.parametrize('field,value',[('charge',2),('mult',3),('native_nao',450),('atoms',40),
    ('extra_columns',[[15,2]]),('geometry_repair_permitted',True),('native_G_qualified',True),
    ('smiles','CCO'),('maximum_peak_GPU_MiB',32768),('native_memory_cap_GiB',16)])
def test_changed_state_resource_or_missing_sulfur_evidence_requires_a_new_pilot(plan,field,value):
    plan[field]=value
    with pytest.raises(ValueError):validate_plan(plan)


def test_extra_columns_use_recorded_center_instead_of_another_Bohr_constant(tmp_path):
    api=load(ROOT/'scripts/validation/run_route_primary_native.py','tested_native_center')
    folder=tmp_path/'SCF_evaluations/0000/0000';folder.mkdir(parents=True)
    exact=np.array([[1.,2.,3.],[-4.,5.,-6.]])/.52917721067
    p=folder/'settings.json';p.write_text(json.dumps(dict(coordinates_Bohr=exact.tolist(),charge=0,mult=1)))
    _,actual=api.recorded_center(tmp_path,0,1,2)
    assert np.array_equal(actual,exact)
    assert not np.array_equal(actual,np.array([[1.,2.,3.],[-4.,5.,-6.]])/.52917721092)
    with pytest.raises(ValueError):api.recorded_center(tmp_path,1,2,2)


@pytest.mark.parametrize('key,charge,mult,state',[
    ('phenothiazine_neu_mode_minus',0,1,'neu'),('phenothiazine_neu_mode_plus',0,1,'neu'),
    ('phenothiazine_ox_stationarity',1,2,'ox')])
def test_refined_native_inputs_keep_full_own_hessian_and_sulfur_contract(plan,key,charge,mult,state):
    plan.update(route_native_profile=REFINED_PROFILE,source_id=key,charge=charge,mult=mult,native_state=state,
        refined_route_native_component=True,electronic_origin_status='route_refined_primary_electronic_independently_checked')
    assert validate_plan(plan)['extra_columns']==[[15,0],[15,1],[15,2]]
    plan['electronic_origin_status']='benchmark_screening_energy_control_independently_checked'
    with pytest.raises(ValueError):validate_plan(plan)


def test_cds_accuracy_control_is_explicit_and_does_not_relax_physical_gates(plan):
    from copy import deepcopy
    plan.update(route_native_profile=CDS_PROFILE,source_id='phenothiazine_neu_mode_minus',
        refined_route_native_component=True,electronic_origin_status='route_refined_primary_electronic_independently_checked',
        CDS_Hessian_numerical_control=True,native_numerical_policy=deepcopy(CDS_POLICY))
    assert validate_plan(plan)['minimum_GPU_memory_MiB']==31744
    assert not plan['native_G_qualified'] and plan['extra_columns']==[[15,0],[15,1],[15,2]]
    plan['native_numerical_policy']['energy_function_changed']=True
    with pytest.raises(ValueError):validate_plan(plan)


def test_numerical_override_cannot_be_hidden_in_a_default_profile(plan):
    plan['CDS_Hessian_numerical_control']=True
    with pytest.raises(ValueError):validate_plan(plan)
