"""A resource probe may change hardware limits, never the selected chemical task."""
import copy
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('crest_point',ROOT/'scripts/validation/qualify_crest_point.py')
MODULE=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(MODULE)


@pytest.fixture
def inputs(tmp_path,monkeypatch):
    root=tmp_path/'project';source=root/'source';inventory=root/'inventory'
    live=source/'live_inputs/c1';live.mkdir(parents=True);(source/'controls').mkdir()
    inventory.mkdir()
    xyz=live/'input.xyz';xyz.write_text('1\nexact source bytes\nH 0 0 0\n')
    manifest=dict(input_xyz=str(xyz),input_sha256={str(xyz):MODULE.api.digest(xyz)})
    (live/'task.json').write_text(json.dumps(manifest))
    task=dict(key='c1',charge=0,mult=1,smiles='source graph',backbone_stereo_tags=['CW','CCW'])
    plan=dict(tasks=[task,dict(task,key='c2')],producer_sha256={},input_sha256=manifest['input_sha256'],
        solvent='acetonitrile',quantum_settings=dict(max_opt_steps=400,geometry_SCF_conv_tol=1e-11),
        thermochemistry=dict(temperature_K=298.15),ordered_DFT_batches=[10,20])
    (source/'plan.json').write_text(json.dumps(plan))
    proof=dict(status='screened_DF_inventory_complete_full_SCF_not_run',source=str(source),
        source_sha256={},task_key='c1',required_host_bytes=121000000000)
    (inventory/'status.json').write_text(json.dumps(proof))
    adapter=root/'adapter.py';adapter.write_text('source\n');checker=root/'checker.py';checker.write_text('source\n')
    monkeypatch.setattr(MODULE,'ROOT',root);monkeypatch.setattr(MODULE,'ADAPTER',adapter);monkeypatch.setattr(MODULE,'CHECKER',checker)
    monkeypatch.setattr(MODULE.api,'load_module',lambda *args:SimpleNamespace(checked_controls=lambda *args:{}))
    return source,inventory,root/'campaign',plan,proof


def test_exact_task_geometry_and_protocol_survive_device_profile(inputs):
    source,inventory,out,old,_=inputs
    MODULE.prepare(source,inventory,out)
    new=json.loads((out/'plan.json').read_text())
    assert new['tasks']==old['tasks'][:1]
    for field in ('solvent','quantum_settings','thermochemistry','ordered_DFT_batches'):
        assert new[field]==old[field]
    assert (out/'live_inputs/c1/input.xyz').read_bytes()==(source/'live_inputs/c1/input.xyz').read_bytes()
    assert new['maximum_RSS_GiB']==80 and new['minimum_free_RAM_GiB']==32
    assert new['minimum_large_GPU_MiB']==32768 and new['maximum_memory_pilot_MiB']==31744
    MODULE.api.verify_sources(new['producer_sha256']);MODULE.api.verify_sources(new['input_sha256'])


@pytest.mark.parametrize('change',[dict(source='other-case'),dict(task_key='missing'),dict(status='inventory_running_no_SCF')])
def test_wrong_or_unfinished_inventory_cannot_prepare(inputs,change):
    source,inventory,out,_,proof=inputs
    (inventory/'status.json').write_text(json.dumps(dict(proof,**change)))
    with pytest.raises(ValueError):MODULE.prepare(source,inventory,out)
    assert not out.exists()


def test_changed_source_geometry_is_rejected_before_preparation(inputs):
    source,inventory,out,_,_=inputs
    (source/'live_inputs/c1/input.xyz').write_text('changed')
    with pytest.raises((ValueError,RuntimeError)):MODULE.prepare(source,inventory,out)
    assert not out.exists()
