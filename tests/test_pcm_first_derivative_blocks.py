"""Bounded first-derivative rectangles preserve dense installed PCM values."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("PCM_first_derivative_blocks",
    ROOT / "scripts/pipeline/pcm_first_derivative_blocks.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


@pytest.fixture
def tensors():
    from pyscf.solvent.grad import pcm
    rng = np.random.default_rng(20261023)
    surface = dict(charge_exp=rng.uniform(.5, 2., 19), grid_coords=rng.normal(size=(19, 3)),
        norm_vec=rng.normal(size=(19, 3)))
    def primitive(local, **kwargs):
        n = len(local["charge_exp"])
        D, S, _ = pcm.get_dD_dS(dict(local, switch_fun=np.ones(n)),
            np.zeros((n, 1, 3)), **kwargs)
        return (D.transpose(2,0,1) if D is not None else None, S.transpose(2,0,1))
    expected = primitive(surface, with_D=True, with_S=True)
    metrics = dict(primitive_calls=0, maximum_output_bytes=0, maximum_primitive_points=0,
        maximum_primitive_tensor_bytes=0, rectangles=0)
    blocks = m.FirstDerivativeBlocks(surface, primitive, np, metrics, tile_points=5)
    return expected, blocks, metrics


@pytest.mark.parametrize("rows,columns", [(slice(0,4),slice(10,19)),
    (slice(3,12),slice(4,16)),(slice(3,3),slice(1,7))])
@pytest.mark.parametrize("kind", [0,1])
def test_rectangles_preserve_overlaps_diagonals_transpose_and_sign(tensors, rows, columns, kind):
    expected, blocks, metrics = tensors
    lazy = m.FirstDerivativeTensor(blocks, kind)
    np.testing.assert_allclose(lazy[:,rows,columns], expected[kind][:,rows,columns], atol=1e-12, rtol=1e-12)
    np.testing.assert_allclose((-lazy.transpose(0,2,1))[:,rows,columns],
        -expected[kind].transpose(0,2,1)[:,rows,columns], atol=1e-12, rtol=1e-12)
    assert metrics["maximum_primitive_points"] <= 10


@pytest.mark.parametrize("kind", [0,1])
def test_atom_responses_and_translation_invariance_match_dense_primitive(tensors, kind):
    expected, blocks, _ = tensors
    dense, lazy = expected[kind], m.FirstDerivativeTensor(blocks, kind)
    q = np.linspace(-.3, .8, 19)
    actual, reference = [], []
    for lower, upper in ((0, 4),(4, 13),(13, 19)):
        response = -(lazy[:,:,lower:upper] @ q[lower:upper])
        response[:,lower:upper] += lazy[:,lower:upper,:] @ q
        reference_response = -(dense[:,:,lower:upper] @ q[lower:upper])
        reference_response[:,lower:upper] += dense[:,lower:upper,:] @ q
        actual.append(response)
        reference.append(reference_response)
    np.testing.assert_allclose(actual, reference, atol=1e-12, rtol=1e-12)
    np.testing.assert_allclose(np.sum(actual,axis=0),0.,atol=1e-12)


def test_unbounded_or_unsupported_access_fails_closed(tensors):
    _, blocks, _ = tensors
    blocks.maximum_output_bytes = 100
    lazy = m.FirstDerivativeTensor(blocks, 0)
    with pytest.raises(MemoryError):
        lazy[:,:,:]
    with pytest.raises(ValueError):
        lazy.transpose(1,0,2)
    with pytest.raises(ValueError):
        lazy[:,0,:]
