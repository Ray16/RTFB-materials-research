import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest


ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


pilot = load('external_uma_pilot_test', ROOT/'scripts/pipeline/coordinate_external_uma_pilot.py')
memory = load('external_uma_memory_test', ROOT/'scripts/fleet/record_uma_gpu_memory.py')


def test_absent_parity_evidence_cannot_launch(tmp_path):
    assert pilot.parity_ready(tmp_path) is False
    (tmp_path/'status.json').write_text(json.dumps({'status':'waiting_for_fixed_device_pair'}))
    assert pilot.parity_ready(tmp_path) is False


def test_rejects_diagnostic_or_mismatched_parity(tmp_path,monkeypatch):
    review = tmp_path/'independent_review.json'
    record = dict(status='independent_fixed_device_pair_passed',all_energy_and_force_checks_pass=True,
                  actual_nitro_and_reference_modes_checked=True,same_device_checked=True,
                  exclusive_residency_checked=True,diagnostic_only=True,source_sha256={})
    review.write_text(json.dumps(record))
    (tmp_path/'status.json').write_text(json.dumps(dict(status='independent_fixed_device_pair_passed',
        review=str(review),review_sha256=pilot.api.digest(review))))
    with pytest.raises(ValueError,match='incomplete'):
        pilot.parity_ready(tmp_path)
    record['diagnostic_only']=False;review.write_text(json.dumps(record))
    with pytest.raises(ValueError,match='changed'):
        pilot.parity_ready(tmp_path)


def test_gpu_memory_parser_counts_only_owned_process_and_rejects_invalid():
    rows = 'GPU-a, 11, 123\nGPU-b, 11, 900\nGPU-a, 12, 700\n'
    assert memory.own_memory(rows,'GPU-a',11) == 123
    assert memory.own_memory(rows,'GPU-b',99) is None
    with pytest.raises(ValueError):
        memory.own_memory('GPU-a, 11, nan','GPU-a',11)


def test_qualified_device_is_source_bound(tmp_path):
    nitro = tmp_path/'nitro'; nitro.mkdir()
    source = nitro/'status.json'
    source.write_text(json.dumps({'mode':'nitro','device_name':'NVIDIA H100 NVL','capability':[9,0]}))
    record = dict(status='independent_fixed_device_pair_passed',all_energy_and_force_checks_pass=True,
                  actual_nitro_and_reference_modes_checked=True,same_device_checked=True,
                  exclusive_residency_checked=True,diagnostic_only=False,
                  source_sha256={str(source):pilot.api.digest(source)})
    review=tmp_path/'independent_review.json';review.write_text(json.dumps(record))
    (tmp_path/'status.json').write_text(json.dumps(dict(status='independent_fixed_device_pair_passed',
        review=str(review),review_sha256=pilot.api.digest(review))))
    assert pilot.qualified_parity_device(tmp_path)=={'name':'NVIDIA H100 NVL','capability':(9,0)}
    source.write_text('{}')
    with pytest.raises(RuntimeError, match='Source or input changed'):
        pilot.qualified_parity_device(tmp_path)


def test_launch_refuses_foreign_preflight_and_partial_output(tmp_path,monkeypatch):
    campaign=tmp_path/'campaign';campaign.mkdir();observer=tmp_path/'observer';observer.mkdir()
    policy=tmp_path/'policy.json';policy.write_text(json.dumps({'nitro_python':'/python'}))
    monkeypatch.setattr(pilot,'POLICY',policy)
    monkeypatch.setattr(pilot,'qualified_parity_device',lambda p:{'name':'device','capability':(9,0)})
    monkeypatch.setattr(pilot.api,'load_module',lambda *a:SimpleNamespace(validate_campaign=lambda c:None,GROUP='g'))
    monkeypatch.setattr(pilot.socket,'gethostname',lambda:'this-host')
    preflight=tmp_path/'preflight.json'
    record={'ok':True,'python':'/python','host':'other-host','GPU_imports_blocked_in_CPU_preflight':True,'GPU_modules_loaded':[]}
    preflight.write_text(json.dumps(record))
    with pytest.raises(ValueError,match='preflight'):
        pilot.launch_local(campaign,observer,tmp_path/'out',1,preflight)
    record['host']='this-host';preflight.write_text(json.dumps(record))
    (campaign/'output/nitro/g/neu').mkdir(parents=True)
    with pytest.raises(ValueError,match='preserved'):
        pilot.launch_local(campaign,observer,tmp_path/'out',1,preflight)
