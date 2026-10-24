import importlib.util
from pathlib import Path

import numpy as np
import pytest


source = Path(__file__).resolve().parents[1] / 'scripts/pipeline/smd_weighted_response.py'
spec = importlib.util.spec_from_file_location('weighted_response_test', source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


@pytest.mark.parametrize('transpose', [False, True])
@pytest.mark.parametrize('side', ['none', 'left', 'right', 'both'])
def test_signed_nonsymmetric_weighted_contractions(transpose, side):
    rng = np.random.default_rng(24)
    matrix = rng.normal(size=(9, 7))
    left = rng.normal(size=9) if side in ('left', 'both') else None
    right = rng.normal(size=7) if side in ('right', 'both') else None
    stats = dict(matrix_tiles=0, maximum_matrix_tile_bytes=0, bounded_response_contracts=0)
    operator = module.ScaledMatrix(matrix, np, left, right, tile=2, metrics=stats)
    dense = matrix.copy()
    if left is not None:
        dense *= left[:, None]
    if right is not None:
        dense *= right
    if transpose:
        operator, dense = operator.T, dense.T
    vector = rng.normal(size=dense.shape[1])
    batch = rng.normal(size=(dense.shape[1], 5))
    np.testing.assert_allclose(operator @ vector, dense @ vector, atol=2e-14)
    np.testing.assert_allclose(operator @ batch, dense @ batch, atol=2e-14)
    response = rng.normal(size=(4, 3, dense.shape[1]))
    fallback = lambda *_: pytest.fail('The response path must use bounded matrix products')
    got = module.response_contract('ij,Adj->Adi', operator, response,
        xp=np, fallback=fallback, metrics=stats)
    np.testing.assert_allclose(got, np.einsum('ij,Adj->Adi', dense, response), atol=2e-14)
    response = rng.normal(size=(4, 3, dense.shape[0]))
    got = module.response_contract('Adi,ij->Adj', response, operator,
        xp=np, fallback=fallback, metrics=stats)
    np.testing.assert_allclose(got, np.einsum('Adi,ij->Adj', response, dense), atol=2e-14)
    assert stats['maximum_matrix_tile_bytes'] <= 2 * max(matrix.shape) * 8


def test_shape_precision_and_nonmatrix_contract_guards():
    with pytest.raises(ValueError, match='FP64'):
        module.ScaledMatrix(np.ones((3, 3), dtype=np.float32), np)
    with pytest.raises(ValueError, match='weights'):
        module.ScaledMatrix(np.ones((3, 3)), np, left=np.ones(2))
    operator = module.ScaledMatrix(np.eye(3), np)
    with pytest.raises(ValueError, match='matching surface'):
        operator @ np.ones(4)
    answer = module.response_contract('i,i->', np.ones(3), np.ones(3),
        xp=np, fallback=np.einsum, metrics={})
    assert answer == 3


def test_source_transform_is_scoped_and_fail_closed():
    body = ('        DA = D*A\n' + "        AS = (A * S.T).T # It's just diag(A) @ S\n") * 2
    changed = module.transformed_source(body)
    assert changed.count('DA = _scaled_matrix(D, right=A)') == 1
    assert changed.count('AS = _scaled_matrix(S, left=A)') == 1
    assert changed.count('DA = D*A') == 1
    with pytest.raises(ValueError, match='template'):
        module.transformed_source(body + '        DA = D*A')
