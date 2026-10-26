import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('extend_manuscript_evidence',ROOT/'manuscript/extend_evidence.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


@pytest.fixture
def snapshot(tmp_path,monkeypatch):
    old=tmp_path/'manuscript/evidence/old';old.mkdir(parents=True)
    figures=tmp_path/'results/figures/manuscript/old';figures.mkdir(parents=True)
    (old/'benchmark.json').write_text('{"value":1}')
    (figures/'benchmark.pdf').write_bytes(b'fixed figure bytes')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    files=[old/'benchmark.json',figures/'benchmark.pdf']
    (old/'manifest.json').write_text(json.dumps({'snapshot_sha256':{str(p.relative_to(tmp_path)):sha(p) for p in files}}))
    pointer=old.parent/'current.json';pointer.write_text(json.dumps({'manifest':'manuscript/evidence/old/manifest.json','figures':'results/figures/manuscript/old'}))
    new=tmp_path/'results/validation/dipea_experimental_anchor_2026-10-26';new.mkdir(parents=True)
    (new/'summary.json').write_text(json.dumps({'source_sha256':{}}));(new/'tests.xml').write_text('<testsuite tests="3"/>')
    script=tmp_path/'scripts/validation/assess_dipea_experimental_anchor.py';script.parent.mkdir(parents=True);script.write_text('# fixture')
    monkeypatch.setattr(module,'ROOT',tmp_path);monkeypatch.setattr(module.subprocess,'check_output',lambda *a,**kw:'test-commit\n')
    return old,pointer,figures


def test_extension_preserves_original_snapshot_and_does_not_switch_current_pointer(snapshot):
    old,pointer,figures=snapshot;before=pointer.read_bytes();result=module.extend('new')
    manifest=json.loads((module.ROOT/result['manifest']).read_text())
    assert manifest['parent_manifest_sha256']==module.sha(old/'manifest.json')
    assert (module.ROOT/'manuscript/evidence/new/benchmark.json').read_bytes()==(old/'benchmark.json').read_bytes()
    assert (module.ROOT/'results/figures/manuscript/new/benchmark.pdf').read_bytes()==(figures/'benchmark.pdf').read_bytes()
    assert pointer.read_bytes()==before
    with pytest.raises(FileExistsError):module.extend('new')


def test_changed_parent_is_rejected_before_any_new_epoch_is_created(snapshot):
    old,pointer,_=snapshot;(old/'benchmark.json').write_text('{"changed":true}')
    with pytest.raises(ValueError,match='Parent snapshot changed'):module.extend('new')
    assert not (module.ROOT/'manuscript/evidence/new').exists()
    assert 'old/manifest.json' in pointer.read_text()
