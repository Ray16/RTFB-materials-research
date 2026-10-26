import importlib.util
import os
from pathlib import Path

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('portable_thermal_entrypoint',ROOT/'scripts/pipeline/run_portable_thermal.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_thread_caps_follow_mode_before_numerical_imports(monkeypatch):
    keys=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
    for key in keys:monkeypatch.setenv(key,'before')
    assert module.configure_thread_pools('projected',['--threads','3'])==3
    assert all(os.environ[k]=='3' for k in keys)
    assert module.configure_thread_pools('probe',[])==1
    assert all(os.environ[k]=='1' for k in keys)
    assert module.configure_thread_pools('native',['--threads=2'])==2
    assert all(os.environ[k]=='2' for k in keys)


def test_nonpositive_thermal_worker_threads_are_rejected():
    with pytest.raises(ValueError,match='Positive'):module.configure_thread_pools('native',['--threads','0'])
