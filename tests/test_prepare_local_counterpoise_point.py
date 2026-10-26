import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest
from redox.workflows import resources

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('prepare_local_counterpoise_point',ROOT/'scripts/pipeline/prepare_local_counterpoise_point.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
guard_spec=importlib.util.spec_from_file_location('guarded_counterpoise_point',ROOT/'scripts/pipeline/run_local_counterpoise_point.py')
guard=importlib.util.module_from_spec(guard_spec);guard_spec.loader.exec_module(guard)


def test_incomplete_control_refuses_preparation_without_creating_outputs(tmp_path):
    construction=tmp_path/'construction';construction.mkdir()
    (construction/'manifest.json').write_text(json.dumps({'source_sha256':{},'entries':[]}))
    checked=tmp_path/'construction_review.json'
    checked.write_text(json.dumps({'source_sha256':{},'status':'fixed_parent_local_motif_construction_independently_checked',
        'primary_candidate_admission':False,'models':[{'key':'pose','parent_coordinates_preserved':True}]}))
    unqualified=tmp_path/'control.json';unqualified.write_text(json.dumps({'source_sha256':{},'status':'running'}))
    out=tmp_path/'output'
    with pytest.raises(ValueError,match='checked CPU/GPU'):
        module.prepare(construction,checked,unqualified,'pose',out)
    assert not out.exists()


def test_gpu_resource_gate_rejects_16gb_and_only_delegates_a_free_32gb_device(monkeypatch,tmp_path):
    invoked=[];memory={'value':15360}
    def fake_load(name,path):
        if 'device_probe' in name:return SimpleNamespace(device=lambda gpu:{'memory_MiB':memory['value']})
        return SimpleNamespace(launch=lambda case,gpu,preflight:invoked.append((case,gpu,preflight)) or {'status':'delegated'})
    monkeypatch.setattr(module.api,'load_module',fake_load)
    monkeypatch.setattr(resources,'free_slots',lambda hosts:[])
    with pytest.raises(ValueError,match='32-GB'):module.launch_32gb(tmp_path,0,tmp_path)
    memory['value']=32768
    assert module.launch_32gb(tmp_path,0,tmp_path)['status']=='atomic_placement_unavailable'
    import socket
    monkeypatch.setattr(resources,'free_slots',lambda hosts:[(socket.gethostname().split('.')[0],0)])
    assert module.launch_32gb(tmp_path,0,tmp_path)['status']=='delegated'
    assert len(invoked)==1


def test_xyz_requires_a_finite_complete_geometry(tmp_path):
    p=tmp_path/'a.xyz';p.write_text('2\ncomment\nH 0 0 0\nH 0 0 nan\n')
    with pytest.raises(ValueError,match='Nonfinite'):module.read_xyz(p)
    p.write_text('2\ncomment\nH 0 0 0\n')
    with pytest.raises(ValueError,match='atom count'):module.read_xyz(p)


def test_reserved_worker_rechecks_32gb_and_gas_scope():
    plan={'minimum_GPU_MiB':32768,'primary_candidate_admission':False,'solvent_counterpoise_computed':False}
    with pytest.raises(ValueError):guard.validate_device(plan,{'memory_MiB':15360})
    guard.validate_device(plan,{'memory_MiB':32768})
    with pytest.raises(ValueError):guard.validate_device({**plan,'solvent_counterpoise_computed':True},{'memory_MiB':32768})
