import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('historical_snapshot_review',ROOT/'scripts/release/reconcile_historical_snapshot.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_verified_archive_reconstructs_changed_publication_input(tmp_path,monkeypatch):
    monkeypatch.setattr(module,'ROOT',tmp_path)
    original=tmp_path/'active.txt';original.write_text('new')
    archived=tmp_path/'old.txt';archived.write_text('original')
    h=hashlib.sha256(archived.read_bytes()).hexdigest()
    record=tmp_path/'record.json';record.write_text(json.dumps({'source_sha256':{'active.txt':h}}))
    mapping=tmp_path/'mapping.json';mapping.write_text(json.dumps({'active.txt':'old.txt'}))
    out=tmp_path/'review.json';r=module.reconcile([record],mapping,out)
    assert r['checked_inputs']==1 and r['records'][0]['retained']=='old.txt'
    with pytest.raises(FileExistsError):module.reconcile([record],mapping,out)


def test_changed_archive_does_not_reconstitute_an_old_source(tmp_path,monkeypatch):
    monkeypatch.setattr(module,'ROOT',tmp_path)
    archived=tmp_path/'old.txt';archived.write_text('changed')
    record=tmp_path/'record.json';record.write_text(json.dumps({'source_sha256':{'active.txt':hashlib.sha256(b'original').hexdigest()}}))
    mapping=tmp_path/'mapping.json';mapping.write_text(json.dumps({'active.txt':'old.txt'}));out=tmp_path/'review.json'
    with pytest.raises(ValueError,match='No exact historical source'):module.reconcile([record],mapping,out)
    assert not out.exists()
