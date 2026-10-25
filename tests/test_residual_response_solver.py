import importlib.util
from pathlib import Path

import numpy as np
import pytest

PATH = Path(__file__).parents[1] / 'scripts/pipeline/residual_response_solver.py'
spec = importlib.util.spec_from_file_location('response_solver_test', PATH)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


@pytest.mark.parametrize('multiple', [False, True])
def test_large_operator_eigenvalues_solve_original_equations(multiple):
    # A perfectly solvable response can have ||A|| > 1. The norm of an
    # unnormalized Krylov basis is not the original-equation residual.
    operator = np.array([[30., 5., 0.], [0., 10., 2.], [1., 0., -2.]])
    rhs = np.array([[1., 2., 3.], [-2., 1., 4.]]) if multiple else np.array([1., 2., 3.])
    records = []
    result = module.residual_krylov(lambda x: x @ operator.T, rhs,
        tol=1e-9, max_cycle=30, record=records.append)
    expected = np.linalg.solve(np.eye(3) + operator, rhs.T).T
    assert np.allclose(result, expected, rtol=0, atol=1e-11)
    assert records[0]['true_residual_norm'] <= 1e-11
    assert records[0]['level_shift_added'] is False


def test_nonconverged_response_cannot_be_accepted():
    rng = np.random.default_rng(103)
    operator = rng.normal(size=(20, 20))
    with pytest.raises(RuntimeError, match='residual check'):
        module.residual_krylov(lambda x: x @ operator, np.ones(20), max_cycle=1)


@pytest.mark.parametrize('rhs', [np.ones(3, dtype=np.float32), np.array([np.nan]), np.ones((2, 2, 2))])
def test_invalid_response_inputs_fail_closed(rhs):
    with pytest.raises(ValueError, match='FP64'):
        module.residual_krylov(lambda x: x, rhs)


def test_operator_cannot_return_wrong_dimensions():
    with pytest.raises(ValueError, match='dimensions'):
        module.residual_krylov(lambda x: np.zeros(1), np.ones(3))
