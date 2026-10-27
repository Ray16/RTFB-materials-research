import importlib.util
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def load(name,path):
    s=importlib.util.spec_from_file_location(name,ROOT/path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def test_native_adapter_compiles_against_current_queue_without_gpu_imports():
    module=load('bounded_native','scripts/pipeline/coordinate_native_bounded.py')
    queue=module.api.load_module('native_queue_actual_template',module.ORIGINAL)
    original=queue.cpu_preflight
    module.adapt(queue)
    assert queue.SCRIPT==Path(module.__file__).resolve() and queue.cpu_preflight is original


def test_native_reconciliation_preserves_attempts_and_raw_data(tmp_path):
    m=load('native_intents','scripts/pipeline/reconcile_native_intents.py')
    case=tmp_path/'case';case.mkdir();raw=tmp_path/'raw'
    assert m.no_attempt(case,raw)
    for name in m.FORBIDDEN:
        p=case/name;p.mkdir();assert not m.no_attempt(case,raw);p.rmdir()
    raw.mkdir();assert not m.no_attempt(case,raw)
