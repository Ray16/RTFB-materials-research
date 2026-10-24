import importlib.util
import json
import hashlib
from pathlib import Path

import pytest


spec = importlib.util.spec_from_file_location('restore_environment',
    Path(__file__).resolve().parents[1] / 'scripts/release/restore_quantum_environment.py')
restore = importlib.util.module_from_spec(spec)
spec.loader.exec_module(restore)


def test_local_project_and_cuda_index_are_separate():
    other, torch = restore.split_lock('redox==0.1.0\ntorch==2.8.0+cu128\nnumpy==2.2.6\n')
    assert other == ['numpy==2.2.6']
    assert torch == 'torch==2.8.0+cu128'


@pytest.mark.parametrize('extra', ['numpy>=2', 'numpy @ https://example.org/a.whl',
    'numpy==2.2.6; python_version>="3"', 'my_pkg==1\nmy-pkg==2'])
def test_unpinned_or_ambiguous_records_are_refused(extra):
    with pytest.raises(ValueError):
        restore.split_lock('redox==0.1.0\ntorch==2.8.0+cu128\n' + extra)


@pytest.mark.parametrize('text', ['torch==2.8.0+cu128', 'redox==0.1.0\ntorch==2.8.0',
    'redox==9\ntorch==2.8.0+cu128'])
def test_missing_local_source_or_changed_cuda_build_is_refused(text):
    with pytest.raises(ValueError):
        restore.split_lock(text)


@pytest.mark.parametrize('status', ['third_party_install',
    'fresh_environment_installation_failed_preserved', 'not_started'])
def test_unfinished_or_dependency_failed_environment_cannot_enter_quantum_validation(tmp_path, status):
    path = Path(__file__).resolve().parents[1] / 'scripts/validation/qualify_restored_quantum.py'
    spec = importlib.util.spec_from_file_location('qualify_restored', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    (tmp_path / 'status.json').write_text(json.dumps(dict(status=status)))
    with pytest.raises(ValueError, match='finish successfully'):
        module.environment(tmp_path)


def closure_fixture(tmp_path, relative='evidence/control.json'):
    path = Path(__file__).resolve().parents[1] / 'scripts/release/certificate_closure.py'
    spec = importlib.util.spec_from_file_location('certificate_closure', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    source, checkout = tmp_path / 'source', tmp_path / 'checkout'
    (source / 'evidence').mkdir(parents=True)
    (source / 'evidence/control.json').write_text('{"verified":true}\n')
    expected = hashlib.sha256((source / 'evidence/control.json').read_bytes()).hexdigest()
    (checkout / 'config').mkdir(parents=True)
    (checkout / 'config/projected_thermal_validation.json').write_text(json.dumps(
        dict(evidence_sha256={relative: expected})))
    (checkout / 'config/native_thermal_validation.json').write_text(json.dumps(dict(controls=[])))
    return module, source, checkout, tmp_path / 'closure.json'


def test_certificate_backing_records_are_included_and_rechecked(tmp_path):
    module, source, checkout, manifest = closure_fixture(tmp_path)
    proof = module.closure(source, checkout, manifest)
    assert (checkout / 'evidence/control.json').read_bytes() == (source / 'evidence/control.json').read_bytes()
    assert proof['scientific_certificate_changed'] is False
    assert module.closure(source, checkout, manifest) == proof


def test_existing_changed_certificate_dependency_is_never_overwritten(tmp_path):
    module, source, checkout, manifest = closure_fixture(tmp_path)
    (checkout / 'evidence').mkdir()
    target = checkout / 'evidence/control.json'
    target.write_text('changed')
    with pytest.raises(ValueError, match='qualified hash'):
        module.closure(source, checkout, manifest)
    assert target.read_text() == 'changed'
    assert not manifest.exists()


@pytest.mark.parametrize('relative', ['../escaped.json', '/tmp/escaped.json'])
def test_certificate_dependency_cannot_escape_snapshot(tmp_path, relative):
    module, source, checkout, manifest = closure_fixture(tmp_path, relative)
    with pytest.raises(ValueError, match='escapes'):
        module.closure(source, checkout, manifest)
