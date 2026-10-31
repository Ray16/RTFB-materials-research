"""Saved solutions must satisfy the freshly evaluated original response equation."""
import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from redox.qm.response_resume import checked_blocks


@pytest.fixture
def checkpoint(tmp_path):
    path=Path(__file__).parent/'fixtures/response_solver_reference.py'
    spec=importlib.util.spec_from_file_location('frozen_resumable_solver',path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    matrix=np.array([[.2,.1],[.1,.3]]);rhs=np.array([[1.,2.],[2.,-1.]])
    solution=module.residual_krylov(lambda x:x@matrix,rhs,tol=1e-10)
    source=tmp_path/'source';source.mkdir();saved=source/'block_000.npy';np.save(saved,solution)
    metadata=dict(rhs_sha256=hashlib.sha256(rhs.tobytes()).hexdigest(),original_solver_source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        checkpoint=dict(path=str(saved),file_sha256=hashlib.sha256(saved.read_bytes()).hexdigest(),array_sha256=hashlib.sha256(solution.tobytes()).hexdigest(),shape=list(solution.shape)))
    (source/'block_000.json').write_text(json.dumps(metadata))
    return module,matrix,rhs,solution,source,tmp_path


def test_matching_candidate_is_reused_only_after_a_fresh_operator_call(checkpoint):
    module,matrix,rhs,solution,source,root=checkpoint;rows=[];calls=[];target=SimpleNamespace(residual_krylov=module.residual_krylov)
    def operator(x):calls.append(x.copy());return x@matrix
    with checked_blocks(target,source,root/'replay') as results:
        actual=target.residual_krylov(operator,rhs,tol=1e-10,record=rows.append)
    np.testing.assert_array_equal(actual,solution)
    assert len(calls)==1 and results[0]['reused'] and rows[0]['true_residual_norm']<1e-12


def test_changed_operator_forces_a_real_solve_instead_of_using_stale_solution(checkpoint):
    module,matrix,rhs,_,source,root=checkpoint;rows=[];target=SimpleNamespace(residual_krylov=module.residual_krylov)
    operator=lambda x:x@(matrix*2)
    with checked_blocks(target,source,root/'replay') as results:
        actual=target.residual_krylov(operator,rhs,tol=1e-10,record=rows.append)
    expected=module.residual_krylov(operator,rhs,tol=1e-10)
    np.testing.assert_allclose(actual,expected,rtol=0,atol=1e-13)
    assert not results[0]['reused'] and rows[-1]['true_residual_norm']<1e-12


def test_changed_payload_is_rejected_without_borrowing_a_solution(checkpoint):
    module,matrix,rhs,solution,source,root=checkpoint
    np.save(source/'block_000.npy',solution*2);target=SimpleNamespace(residual_krylov=module.residual_krylov)
    with pytest.raises(ValueError):
        with checked_blocks(target,source,root/'replay'):
            target.residual_krylov(lambda x:x@matrix,rhs,tol=1e-10,record=lambda row:None)
