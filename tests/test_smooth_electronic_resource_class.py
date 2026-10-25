"""Transfer execution bounds without transferring molecular/model evidence."""
from copy import deepcopy
import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('smooth_electronic_class_test',ROOT/'scripts/pipeline/prepare_smooth_electronic_batch.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


@pytest.fixture
def state():
    task=dict(smiles='C[C@H](Cl)F',atoms=8,active_basis_functions=128,charge=0,mult=1)
    point=dict(tasks=[task],experimental_method_id='smooth-model',experimental_model={'surface':'smooth'},
        quantum_settings={'geometry_SCF_conv_tol':1e-12},runtime_profile='runtime.json',
        maximum_host_RSS_GiB=80,maximum_peak_GPU_MiB=31744)
    review=dict(status='alternate_electronic_geometry_and_resource_checked_native_G_pending',
        portable_method_id='smooth-model',native_Hessian_or_G_qualified=False,primary_candidate_admission=False,
        peak_host_RSS_bytes=75*2**30,observed_peak_GPU_MiB=20000)
    return deepcopy(point),point,review


def test_different_stereo_configuration_can_share_only_resource_class(state):
    plan,point,review=state;plan['tasks'][0]['smiles']='C[C@@H](Cl)F'
    module.resource_compatible(plan,point,review)


@pytest.mark.parametrize('change',['constitution','basis','charge','method','surface','runtime',
                                  'settings','G_promotion','RAM','nan'])
def test_incompatible_resource_or_science_context_is_rejected(state,change):
    plan,point,review=state
    if change=='constitution':plan['tasks'][0]['smiles']='FCCCl'
    elif change=='basis':plan['tasks'][0]['active_basis_functions']+=1
    elif change=='charge':plan['tasks'][0]['charge']=-1
    elif change=='method':plan['experimental_method_id']='legacy-model'
    elif change=='surface':plan['experimental_model']={'surface':'legacy'}
    elif change=='runtime':plan['runtime_profile']='another-runtime.json'
    elif change=='settings':plan['quantum_settings']['geometry_SCF_conv_tol']=1e-9
    elif change=='G_promotion':review['native_Hessian_or_G_qualified']=True
    elif change=='RAM':review['peak_host_RSS_bytes']=81*2**30
    else:review['peak_host_RSS_bytes']=float('nan')
    with pytest.raises(ValueError):module.resource_compatible(plan,point,review)
