import importlib.util
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_error_capture_preserves_original_controller_template():
    adapter = load('capacity_v2_test', ROOT / 'scripts/pipeline/coordinate_with_capacity_preflight_v2.py')
    upstream = load('electronic_controller_test', ROOT / 'scripts/pipeline/coordinate_smooth_electronic_batch.py')
    patched = adapter.coordinate_with_durable_errors(upstream)
    assert callable(patched)
    assert 'placement_error.json' in patched.__code__.co_consts
    assert 'placement_response_unverified_preserved' in patched.__code__.co_consts
