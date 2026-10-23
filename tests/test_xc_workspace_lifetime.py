"""The retained final AO batch is released after rho is copied, not before."""
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import weakref

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("XC_lifetime", ROOT / "scripts/pipeline/xc_workspace_execution.py")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def test_final_density_batch_can_be_released_without_changing_rho():
    body = '''def get_exc_full_response():
    rho = np.empty(3)
    g1 = 0
    ngrids = 3
    for ao, idx, weight, _ in block_loop():
        dms_masked = np.ones(1)
        rho[:] = ao * 2
        g1 += len(weight)
    assert g1 == ngrids

    exc, vxc = evaluate(rho)
    return exc, vxc
'''
    alive = []
    seen = []

    def block_loop():
        ao = np.array([1., 2., 3.])
        alive.append(weakref.ref(ao))
        yield ao, np.arange(3), np.ones(3), None

    def evaluate(rho):
        seen.append(alive[-1]() is not None)
        return rho.sum(), rho.copy()

    ns = dict(np=np, block_loop=block_loop, evaluate=evaluate,
        cupy=SimpleNamespace(cuda=SimpleNamespace(get_current_stream=lambda:SimpleNamespace(synchronize=lambda:None))),
        _workspace_trace=lambda stage:None)
    exec(body, ns)
    before = ns["get_exc_full_response"]()
    exec(MODULE.release_density_batch(body), ns)
    after = ns["get_exc_full_response"]()
    assert seen == [True, False]
    assert before[0] == after[0] and np.array_equal(before[1], after[1])


def test_changed_density_loop_template_is_rejected():
    with pytest.raises(ValueError, match="template changed"):
        MODULE.release_density_batch("def get_exc_full_response(): return None")
