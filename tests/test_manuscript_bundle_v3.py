import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tarfile
from types import SimpleNamespace

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('bundle_v3_test',ROOT/'scripts/release/prepare_manuscript_bundle_v3.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_current_epoch_must_be_included_and_match_committed_bytes(monkeypatch,tmp_path):
    manifest_path='manuscript/evidence/epoch/manifest.json';item='manuscript/evidence/epoch/value.json'
    contents=b'{"value": 1}\n';digest=hashlib.sha256(contents).hexdigest()
    blobs={'manuscript/evidence/current.json':json.dumps({'manifest':manifest_path}).encode(),
           manifest_path:json.dumps({'snapshot_sha256':{item:digest}}).encode(),item:contents}
    def fake_git(args,**kwargs):return blobs[args[2].split(':',1)[1]]
    monkeypatch.setattr(module.subprocess,'check_output',fake_git)
    closure=SimpleNamespace(included=lambda path,declared:any(path==p or path.startswith(p+'/') for p in declared))
    config={'source_paths':['manuscript/evidence/epoch']}
    assert module.check_epoch(tmp_path,config,'HEAD',closure)==1
    with pytest.raises(ValueError,match='excluded'):
        module.check_epoch(tmp_path,{'source_paths':[]},'HEAD',closure)
    blobs[item]=b'changed'
    with pytest.raises(ValueError,match='differs'):
        module.check_epoch(tmp_path,config,'HEAD',closure)


def test_v3_declares_five_distinct_packages_and_current_tempo_epoch():
    config=json.loads((ROOT/'config/manuscript_reviewer_bundle_v3.json').read_text())
    assert {p['id'] for p in config['packages']}=={'primary_review','first_site','numerical_controls','neighbor_cycle','tempo_addendum'}
    assert 'manuscript/evidence/2026-10-26_tempo_exact_support' in config['source_paths']
    assert config['data_license'] is None and config['author_metadata'] is None


def test_v3_extractor_rejects_escaping_or_symlink_members(tmp_path):
    path=ROOT/'scripts/release/validate_manuscript_bundle_v3.py'
    s=importlib.util.spec_from_file_location('bundle_v3_validator_test',path)
    validator=importlib.util.module_from_spec(s);s.loader.exec_module(validator)
    for index,kind in enumerate(('escape','symlink')):
        data=io.BytesIO()
        with tarfile.open(fileobj=data,mode='w') as tar:
            member=tarfile.TarInfo('../outside' if kind=='escape' else 'root/link')
            if kind=='symlink':member.type=tarfile.SYMTYPE;member.linkname='../outside';tar.addfile(member)
            else:member.size=1;tar.addfile(member,io.BytesIO(b'x'))
        data.seek(0);out=tmp_path/str(index);out.mkdir()
        with tarfile.open(fileobj=data,mode='r:') as tar,pytest.raises(ValueError):
            validator.safe_extract(tar,out,expected_root='root')
    assert not (tmp_path/'outside').exists()
