"""A real retained state must survive relocation and reject altered raw inputs."""
import importlib.util
import io
import json
from pathlib import Path
import tarfile

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
definition = importlib.util.spec_from_file_location('portable_native_cache', ROOT/'scripts/release/native_cache.py')
cache = importlib.util.module_from_spec(definition)
definition.loader.exec_module(cache)


@pytest.fixture
def package(tmp_path):
    root = tmp_path/'moved'
    root.mkdir()
    with tarfile.open(ROOT/'tests/data/native_cache_one_state.tar.gz', 'r:gz') as archive:
        for item in archive:
            assert item.isfile()
            target = cache.local(root, item.name)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(archive.extractfile(item).read())
    return root


def save(root, manifest):
    (root/'manifest.json').write_text(json.dumps(manifest))


def replace_object(root, manifest, row, role, data):
    """Deliberately rehash an altered dataset to test the numerical gates."""
    digest = __import__('hashlib').sha256(data).hexdigest()
    relative = 'objects/'+digest
    (root/relative).write_bytes(data)
    manifest['files_sha256'][relative] = digest
    row['files'][role] = relative
    if role == 'review':
        row['original_review_sha256'] = digest
    else:
        original = row['original_selected_sources'][role]
        original['sha256'] = digest
        report = json.loads((root/row['files']['review']).read_text())
        report['source_sha256'][original['path']] = digest
        replace_object(root, manifest, row, 'review', json.dumps(report).encode())


def test_real_raw_state_reconstructs_without_original_sources(package):
    result = cache.replay(package)
    assert result['states'] == 1
    row = result['records'][0]
    assert row['key'] == 'all_same__occ_none__conf001'
    assert row['frequency_error_cm'] < 1e-5
    assert max(row['nuclear_errors_eV'].values()) < 1e-8
    assert not result['full_historical_source_closure_revalidated']
    assert not result['live_execution_revalidated'] and not result['quantum_calculations_run']


def test_changed_raw_bytes_refused(package):
    manifest = json.loads((package/'manifest.json').read_text())
    (package/manifest['states'][0]['files']['gradient']).write_bytes(b'changed')
    with pytest.raises(ValueError, match='artifact differs'):
        cache.replay(package)


@pytest.mark.parametrize('mutation', ['hessian', 'gradient', 'nan_mode', 'nan_G', 'mass', 'method', 'parent_task'])
def test_rehashed_but_physically_inconsistent_inputs_refused(package, mutation):
    manifest = json.loads((package/'manifest.json').read_text())
    row = manifest['states'][0]
    if mutation in ('hessian', 'gradient'):
        array = np.load(package/row['files'][mutation], allow_pickle=False)
        array *= 2
        buffer = io.BytesIO(); np.save(buffer, array)
        replace_object(package, manifest, row, mutation, buffer.getvalue())
    elif mutation == 'parent_task':
        parent = json.loads((package/row['files']['parent']).read_text())
        parent['task']['key'] = 'another exact state'
        replace_object(package, manifest, row, 'parent', json.dumps(parent).encode())
    else:
        native = json.loads((package/row['files']['partition']).read_text())
        if mutation == 'nan_mode':native['nuclear']['frequencies_cm'][0] = float('nan')
        if mutation == 'nan_G':native['G_by_cutoff_eV']['50'] = float('nan')
        if mutation == 'mass':native['nuclear']['atomic_masses_amu'][0] = -1
        if mutation == 'method':native['experimental_method_id'] = 'different solvent model'
        replace_object(package, manifest, row, 'partition', json.dumps(native).encode())
    save(package, manifest)
    with pytest.raises(ValueError):
        cache.replay(package)


def test_missing_displaced_gradient_cannot_be_hidden_by_manifest(package):
    manifest = json.loads((package/'manifest.json').read_text())
    del manifest['states'][0]['files']['evaluation_0001']
    save(package, manifest)
    with pytest.raises((KeyError, ValueError)):
        cache.replay(package)


def test_duplicate_exact_state_cannot_inflate_coverage(package):
    manifest = json.loads((package/'manifest.json').read_text())
    manifest['states'].append(manifest['states'][0])
    save(package, manifest)
    with pytest.raises(ValueError, match='identity'):
        cache.replay(package)


@pytest.mark.parametrize('path', ['/outside', '../outside', 'nested/../../outside'])
def test_artifacts_cannot_open_original_or_parent_paths(package, path):
    with pytest.raises(ValueError):
        cache.local(package, path)
