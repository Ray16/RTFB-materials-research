import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('extend_uv_evidence',ROOT/'manuscript/extend_evidence_with_uv.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


@pytest.fixture
def snapshot(tmp_path,monkeypatch):
    old=tmp_path/'manuscript/evidence/old';old.mkdir(parents=True)
    figures=tmp_path/'results/figures/manuscript/old';figures.mkdir(parents=True)
    (old/'benchmark.json').write_text('{"value":1}');(figures/'benchmark.pdf').write_bytes(b'unchanged figure')
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    (old/'manifest.json').write_text(json.dumps({'snapshot_sha256':{str(p.relative_to(tmp_path)):sha(p) for p in [old/'benchmark.json',figures/'benchmark.pdf']}}))
    pointer=old.parent/'current.json';pointer.write_text(json.dumps({'manifest':'manuscript/evidence/old/manifest.json','figures':'results/figures/manuscript/old'}))
    uv=tmp_path/'results/validation/dipea_uv_anchor_2026-10-26';uv.mkdir(parents=True)
    (uv/'summary.json').write_text(json.dumps({'status':'primary_DIPEA_UV_protonation_constant_source_checked_conditional_comparison','source_sha256':{}}))
    (uv/'independent_review.json').write_text(json.dumps({'status':'independent_UV_source_and_comparison_review_checked','source_sha256':{}}));(uv/'tests.xml').write_text('<testsuite tests="4"/>')
    scripts=tmp_path/'scripts/validation';scripts.mkdir(parents=True)
    for name in ['assess_dipea_uv_anchor.py','review_dipea_uv_anchor.py']:(scripts/name).write_text('# fixture')
    monkeypatch.setattr(module,'ROOT',tmp_path);monkeypatch.setattr(module.subprocess,'check_output',lambda *a,**kw:'test-commit\n')
    return old,pointer,uv


def test_uv_extension_preserves_parent_and_does_not_publish_pointer(snapshot):
    old,pointer,_=snapshot;before=pointer.read_bytes();result=module.extend('new')
    d=json.loads((module.ROOT/result['manifest']).read_text());assert d['parent_manifest_sha256']==module.sha(old/'manifest.json')
    assert (module.ROOT/'manuscript/evidence/new/dipea_uv_anchor.json').exists() and pointer.read_bytes()==before
    with pytest.raises(FileExistsError):module.extend('new')


def test_unreviewed_uv_or_changed_parent_is_rejected_without_publishing(snapshot):
    old,pointer,uv=snapshot;review=uv/'independent_review.json';review.write_text(json.dumps({'status':'running','source_sha256':{}}))
    with pytest.raises(ValueError,match='independent review'):module.extend('new')
    assert not (module.ROOT/'manuscript/evidence/new').exists()
    (old/'benchmark.json').write_text('changed')
    with pytest.raises(ValueError,match='Parent snapshot changed'):module.extend('new')
    assert 'old/manifest.json' in pointer.read_text()
