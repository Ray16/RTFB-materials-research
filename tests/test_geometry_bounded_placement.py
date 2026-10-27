import fcntl
import importlib.util
import json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('bounded_geometry',ROOT/'scripts/pipeline/coordinate_smooth_geometry_bounded.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)


def test_busy_allocation_lock_is_immediate_resource_wait(tmp_path):
    lock_path=tmp_path/'host.lock';out=tmp_path/'phase.json';record={}
    with lock_path.open('a+') as other:
        fcntl.flock(other,fcntl.LOCK_EX|fcntl.LOCK_NB)
        with pytest.raises(mod.HostAllocationBusy):
            with mod.bounded_allocation(lock_path,record,out):
                pytest.fail('Busy allocation guard entered')
    assert json.loads(out.read_text())['phase']=='host_allocation_guard_busy'
    with mod.bounded_allocation(lock_path,record,out):pass
    assert json.loads(out.read_text())['phase']=='host_allocation_guard_acquired'


def test_allocation_lock_released_after_failure(tmp_path):
    path=tmp_path/'lock';out=tmp_path/'phase.json'
    with pytest.raises(ValueError):
        with mod.bounded_allocation(path,{},out):raise ValueError('stop')
    with mod.bounded_allocation(path,{},out):pass


def test_existing_queue_adapts_with_unchanged_scientific_functions():
    base=mod.api.load_module('test_bounded_inventory_base',mod.BASE)
    module=mod.api.load_module('test_bounded_queue',mod.ROOT/'scripts/pipeline/coordinate_geometry_batch.py')
    worker_check=module.owner_live
    mod.adapt(module,{'python':'/python','library_path':'/library'},base.adapt)
    assert module.SCRIPT==Path(mod.__file__).resolve()
    assert module.owner_live is worker_check
    assert callable(module.launch_local) and callable(module.coordinate)
