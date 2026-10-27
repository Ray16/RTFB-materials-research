import hashlib
import importlib.util
import json
from pathlib import Path

import pytest


script = Path(__file__).resolve().parents[1]/'scripts/release/local_motif_raw_archive.py'
spec = importlib.util.spec_from_file_location('local_motif_raw_archive', script)
archive = importlib.util.module_from_spec(spec)
spec.loader.exec_module(archive)


def package(tmp_path, relative, content=b'state'):
    project = tmp_path/'project'
    path = project/relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    manifest = {'original_root':'/original/project', 'files_sha256':{relative:hashlib.sha256(content).hexdigest()}}
    (tmp_path/'manifest.json').write_text(json.dumps(manifest))
    return path


def test_modified_packaged_source_is_rejected(tmp_path):
    path = package(tmp_path, 'records/energy.json')
    path.write_bytes(b'changed')
    with pytest.raises(ValueError, match='hash differs'):
        archive.replay(tmp_path)


def test_external_symlink_is_rejected(tmp_path):
    path = package(tmp_path, 'records/energy.json')
    outside = tmp_path/'outside.json'
    outside.write_bytes(b'state')
    path.unlink()
    path.symlink_to(outside)
    with pytest.raises(ValueError, match='escapes'):
        archive.replay(tmp_path)


def test_parent_traversal_is_rejected(tmp_path):
    (tmp_path/'project').mkdir()
    (tmp_path/'manifest.json').write_text(json.dumps({'original_root':'/original/project',
                                                     'files_sha256':{'../escape':'0'*64}}))
    with pytest.raises(ValueError, match='Unpackaged'):
        archive.replay(tmp_path)
