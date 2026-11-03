"""Portable raw-object verification rejects corruption and invalid mappings."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / 'scripts/release/verify_source_objects.py'
spec = importlib.util.spec_from_file_location('source_object_verifier', SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


@pytest.fixture
def package(tmp_path):
    payload = b'preserved raw evidence\n'
    digest = hashlib.sha256(payload).hexdigest()
    (tmp_path / 'objects').mkdir()
    (tmp_path / 'objects' / digest).write_bytes(payload)
    manifest = dict(objects={digest: len(payload)},
                    files=[dict(path='/unavailable/old/source', sha256=digest)],
                    relocations=[dict(original='/old/source', sha256=digest)])
    (tmp_path / 'manifest.json').write_text(json.dumps(manifest))
    return tmp_path, manifest, digest


def test_relocated_package_needs_no_original_paths(package):
    folder, _, _ = package
    result = module.verify(folder)
    assert result['verified_objects'] == 1
    assert not result['quantum_validation_repeated']
    assert not result['external_dependencies_verified']


def test_changed_raw_bytes_rejected(package):
    folder, _, digest = package
    (folder / 'objects' / digest).write_bytes(b'changed raw evidence!\n')
    with pytest.raises(ValueError, match='Object differs'):
        module.verify(folder)


def test_unbacked_mapping_rejected(package):
    folder, manifest, _ = package
    manifest['relocations'][0]['sha256'] = '0' * 64
    (folder / 'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='Missing relocated'):
        module.verify(folder)


def test_path_traversal_identity_rejected(package):
    folder, manifest, _ = package
    manifest['objects'] = {'../outside': 1}
    (folder / 'manifest.json').write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match='Invalid object identity'):
        module.verify(folder)


def test_symlink_substitution_rejected(package):
    folder, _, digest = package
    obj = folder / 'objects' / digest
    target = folder / 'external'; obj.rename(target); obj.symlink_to(target)
    with pytest.raises(ValueError, match='Object differs'):
        module.verify(folder)


def test_staging_keeps_external_relocations_separate(tmp_path):
    stage_spec = importlib.util.spec_from_file_location(
        'source_object_stager', SCRIPT.with_name('stage_declared_sources.py'))
    stager = importlib.util.module_from_spec(stage_spec)
    stage_spec.loader.exec_module(stager)
    raw = tmp_path / 'raw'; raw.write_bytes(b'raw')
    external = tmp_path / 'runtime'; external.write_bytes(b'external')
    raw_hash = module.sha256(raw); ext_hash = module.sha256(external)
    inv = tmp_path / 'inventory.json'
    inv.write_text(json.dumps(dict(files=[
        dict(path=str(raw), sha256=raw_hash, bytes=3, category='versioned_project_artifact'),
        dict(path=str(external), sha256=ext_hash, bytes=8, category='runtime_dependency')],
        explicit_archive_relocations=[dict(original='/old/runtime', path=str(external), sha256=ext_hash)],
        traversal_scope='test declared mappings')))
    content = tmp_path / 'content.json'
    content.write_text(json.dumps(dict(inventory_sha256=module.sha256(inv), resolved=[])))
    metadata = tmp_path / 'metadata.json'
    metadata.write_text(json.dumps(dict(previous_resolution_sha256=module.sha256(content),
                                        unresolved=[], resolved=[])))
    output = tmp_path / 'package'
    stager.stage(inv, content, metadata, output)
    manifest = json.loads((output / 'manifest.json').read_text())
    assert manifest['relocations'] == []
    assert len(manifest['external_relocations']) == 1
    assert not (output / 'objects' / ext_hash).exists()
    assert module.verify(output)['verified_objects'] == 1
    with pytest.raises(FileExistsError):
        stager.stage(inv, content, metadata, output)
    (output / 'manifest.json').unlink()
    (output / 'verify.py').unlink()
    (output / 'README.md').unlink()
    (output / 'objects' / raw_hash).write_bytes(b'bad')
    with pytest.raises(ValueError, match='Copied evidence differs'):
        stager.stage(inv, content, metadata, output, resume=True)
