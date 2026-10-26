import importlib.util
from pathlib import Path

ROOT=Path(__file__).parents[1]


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def test_resource_adapter_preserves_the_scientific_program_identity():
    adapter=load('local_resource_adapter_test','scripts/pipeline/dispatch_local_motif_energy.py')
    worker=load('original_local_scientific_worker_test','scripts/pipeline/local_motif_energy.py')
    modified=adapter.adapted_launch(worker)
    assert modified.__globals__['__file__']==worker.__file__
    assert modified.__globals__['PROFILE']==worker.PROFILE
    assert modified.__globals__['SETTINGS']==worker.SETTINGS


def test_only_confirmed_pre_spawn_cpu_capacity_errors_can_be_retried(tmp_path):
    controller=load('local_CPU_capacity_classification_test','scripts/pipeline/coordinate_local_motif_energy.py')
    message='RuntimeError: Requested CPU allocation is unavailable or overlaps an existing worker'
    assert controller.pre_spawn_cpu_refusal(tmp_path,message)
    assert not controller.pre_spawn_cpu_refusal(tmp_path,'unattributed RPC failure')
    (tmp_path/'residency_guard').mkdir()
    assert not controller.pre_spawn_cpu_refusal(tmp_path,message)
