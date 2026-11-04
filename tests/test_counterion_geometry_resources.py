"""Geometry-only T4 admission keeps physics, scope, and the Hessian gate fixed."""
from copy import deepcopy

import pytest

from redox.workflows import counterion_geometry_resources as profile
from redox.workflows import resources


def plan():
    return dict(role='primary_counterion_geometry_control', geometry_resource_profile=profile.PROFILE,
        geometry_memory_pilot=True, atoms=42, nao=412, charge=1, mult=1,
        parent_formal_core_charge=2, counterion_charge=-1, smiles=profile.SMILES,
        maximum_host_RSS_GiB=32, minimum_free_GiB=32, minimum_GPU_memory_MiB=15360,
        maximum_peak_GPU_MiB=12288, maximum_wall_seconds=14400,
        maximum_durable_output_bytes=2*2**30, Hessian_computed=False,
        native_G_qualified=False, material_admission=False, method=deepcopy(profile.METHOD),
        optimizer_convergence=deepcopy(profile.CONVERGENCE), source_id='N0_below',
        seed=dict(source='verified negative native mode', sign=-1, maximum_atom_displacement_A=.1),
        fixed_gradient_review='/durable/fixed_gradient_review.json')


@pytest.mark.parametrize('field,value', [('role', 'primary_counterion_thermochemistry'),
    ('atoms', 43), ('nao', 1003), ('charge', 2), ('mult', 2), ('Hessian_computed', True),
    ('native_G_qualified', True), ('maximum_peak_GPU_MiB', 15360), ('material_admission', True),
    ('smiles', profile.SMILES.replace('[Cl-]', '[I-]')), ('fixed_gradient_review', None),
    ('source_id', 'unprepared_new_contact')])
def test_scope_and_resource_changes_need_separate_qualification(field, value):
    data=plan();data[field]=value
    with pytest.raises(ValueError): profile.validate_plan(data)


def test_strict_stationarity_and_primary_smd_cannot_be_relaxed():
    for field in ('method', 'optimizer_convergence', 'seed'):
        data=plan();data[field]={}
        with pytest.raises(ValueError): profile.validate_plan(data)


def test_a_free_T4_is_eligible_only_for_the_declared_profile(monkeypatch):
    monkeypatch.setattr(resources,'inspect_gpu',lambda i:dict(name='Tesla T4',memory_MiB=15360,compute_capability='7.5'))
    assert profile.device_for_plan(plan(),0)['name']=='Tesla T4'
    with pytest.raises(RuntimeError,match='32-GB'): resources.gpu_suitability(0)
    monkeypatch.setattr(resources,'inspect_gpu',lambda i:dict(name='Tesla T4',memory_MiB=15360,compute_capability='7.0'))
    with pytest.raises(RuntimeError): profile.device_for_plan(plan(),0)
