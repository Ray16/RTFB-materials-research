import importlib.util
import json
from pathlib import Path

import pytest


@pytest.fixture
def stream_module(monkeypatch):
    folder=Path(__file__).resolve().parents[1]/'scripts/release'
    monkeypatch.syspath_prepend(str(folder))
    spec=importlib.util.spec_from_file_location('archive_stream',folder/'stream_source_archive.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def inputs(tmp_path,module):
    raw=tmp_path/'raw.dat';raw.write_bytes(b'original scientific bytes\n')
    inv=tmp_path/'inventory.json'
    inv.write_text(json.dumps(dict(status='declared_dependency_inventory_verified',unresolved=[],
        files=[dict(path=str(raw),sha256=module.sha256_file(raw),bytes=raw.stat().st_size,category='research_evidence_or_source_snapshot')],
        explicit_archive_relocations=[],traversal_scope='test')))
    return raw,inv,tmp_path/'package.tar.zst',tmp_path/'metadata'


def test_archive_roundtrip_without_original_data(tmp_path,stream_module):
    raw,inv,archive,meta=inputs(tmp_path,stream_module)
    result=stream_module.create(inv,archive,meta)
    raw.rename(tmp_path/'source_unavailable.dat')
    check=stream_module.verify_archive(archive,meta/'manifest.json')
    assert check['objects']==1 and check['members_verified']==6
    assert not result['uncompressed_copy_created']
    assert not check['numerical_replay_performed']
    with pytest.raises(FileExistsError):stream_module.create(inv,archive,meta)


def test_changed_source_cannot_enter_archive(tmp_path,stream_module):
    raw,inv,archive,meta=inputs(tmp_path,stream_module)
    raw.write_bytes(b'X'*raw.stat().st_size)
    with pytest.raises(ValueError,match='Source bytes changed'):
        stream_module.create(inv,archive,meta)


def test_changed_manifest_rejected_at_archive_verification(tmp_path,stream_module):
    raw,inv,archive,meta=inputs(tmp_path,stream_module)
    stream_module.create(inv,archive,meta)
    path=meta/'manifest.json';data=json.loads(path.read_text());data['inventory_scope']='altered'
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError,match='Archive member bytes differ'):
        stream_module.verify_archive(archive,path)


def test_existing_metadata_directory_preserves_unrelated_records(tmp_path,stream_module):
    raw,inv,archive,meta=inputs(tmp_path,stream_module)
    meta.mkdir();record=meta/'earlier_review.json';record.write_text('{"preserved":true}')
    stream_module.create(inv,archive,meta,existing_metadata=True)
    assert record.read_text()=='{"preserved":true}'
    manifest=json.loads((meta/'manifest.json').read_text())
    assert 'earlier_review.json' not in manifest['support_files_sha256']
    assert stream_module.verify_archive(archive,meta/'manifest.json')['members_verified']==6


def test_existing_metadata_conflict_refuses_before_archive_write(tmp_path,stream_module):
    raw,inv,archive,meta=inputs(tmp_path,stream_module)
    meta.mkdir();(meta/'README.md').write_text('earlier evidence')
    with pytest.raises(FileExistsError,match='conflicting'):
        stream_module.create(inv,archive,meta,existing_metadata=True)
    assert not archive.exists() and (meta/'README.md').read_text()=='earlier evidence'
