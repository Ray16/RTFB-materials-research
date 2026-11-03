"""Small-GPU admission is confined to a declared fixed-gradient validation pilot."""
from copy import deepcopy
from types import SimpleNamespace

import pytest

from redox.workflows import resources
from redox.workflows.counterion_gradient_pilot import PROFILE,ROLE,METHOD,validate_plan,device_for_plan


def plan():
    return dict(role=ROLE,geometry_resource_profile=PROFILE,fixed_gradient_memory_pilot=True,
        atoms=42,nao=412,charge=1,mult=1,method=deepcopy(METHOD),maximum_host_RSS_GiB=8,
        minimum_free_GiB=32,minimum_GPU_memory_MiB=15360,maximum_peak_GPU_MiB=12288,
        maximum_wall_seconds=3600,maximum_durable_output_bytes=256*2**20,
        Hessian_computed=False,geometry_optimization_performed=False,native_G_qualified=False)


def test_full_molecular_default_still_rejects_a_T4(monkeypatch):
    monkeypatch.setattr(resources.subprocess,'run',lambda *a,**k:SimpleNamespace(stdout='0, UUID, Tesla T4, 15360, 7.5\n'))
    assert resources.inspect_gpu(0)['memory_MiB']==15360
    with pytest.raises(RuntimeError,match='32-GB'):resources.gpu_suitability(0)
    assert device_for_plan(plan(),0)['name']=='Tesla T4'


@pytest.mark.parametrize('field,value',[('role','primary_counterion_thermochemistry'),('atoms',43),('nao',1003),
    ('charge',2),('mult',2),('Hessian_computed',True),('geometry_optimization_performed',True),
    ('native_G_qualified',True),('maximum_peak_GPU_MiB',15360)])
def test_the_pilot_cannot_silently_admit_hessians_optimization_or_other_states(field,value):
    data=plan();data[field]=value
    with pytest.raises(ValueError):validate_plan(data)


def test_the_actual_architecture_must_match_the_prospective_control(monkeypatch):
    monkeypatch.setattr(resources,'inspect_gpu',lambda i:dict(name='Tesla V100',memory_MiB=32768,compute_capability='7.0'))
    with pytest.raises(RuntimeError):device_for_plan(plan(),0)
