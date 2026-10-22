"""Bounded contractions reproduce independent dense installed PCM primitives."""

import math
import numpy as np
import pytest

from redox.qm.pcm_bounded import (
    diagonal_contraction,
    CurvatureBlocks,
    OffDiagonalTensor,
)


@pytest.mark.parametrize("atoms,points", [(1, 1), (3, 37), (7, 129)])
def test_diagonal_fusion_matches_dense_formula_for_strided_tensors(atoms, points):
    rng = np.random.default_rng(123)
    first = rng.normal(size=(3, points, atoms))
    second = rng.normal(size=(points, atoms, atoms, 3, 3)).transpose(1, 2, 3, 4, 0)
    switch = rng.uniform(0.2, 1.0, size=points)
    charge = rng.uniform(0.2, 2.0, size=points)
    vector = rng.normal(size=points)
    derivative = first.transpose(2, 0, 1)
    dense = (
        (
            2
            * derivative[:, None, :, None, :]
            * derivative[None, :, None, :, :]
            / switch**3
            - second / switch**2
        )
        * math.sqrt(2 / math.pi)
        * charge
    )
    actual = diagonal_contraction(
        dict(switch_fun=switch, charge_exp=charge), first, second, vector
    )
    np.testing.assert_allclose(actual, dense @ vector, atol=1e-11, rtol=1e-12)


@pytest.fixture
def example():
    from pyscf.solvent.hessian import pcm

    rng = np.random.default_rng(1234)
    surface = dict(
        charge_exp=rng.uniform(0.5, 2.0, 19),
        grid_coords=rng.normal(size=(19, 3)),
        norm_vec=rng.normal(size=(19, 3)),
    )
    metrics = dict(
        uploaded_blocks=0,
        maximum_uploaded_block_bytes=0,
        primitive_calls=0,
        peak_cached_host_bytes=0,
    )
    return pcm, surface, metrics


@pytest.mark.parametrize(
    "rows,columns",
    [
        (slice(0, 4), slice(10, 14)),
        (slice(1, 6), slice(3, 8)),
        (slice(5, 5), slice(2, 4)),
    ],
)
def test_CPU_primitive_blocks_keep_diagonal_and_transpose_semantics(
    example, rows, columns
):
    pcm, surface, metrics = example
    expected = pcm.get_d2D_d2S(surface, with_S=True, with_D=True)
    cache = CurvatureBlocks(surface, pcm.get_d2D_d2S, np.asarray, metrics)
    for kind in (0, 1):
        tensor = OffDiagonalTensor(cache, kind)
        np.testing.assert_allclose(
            tensor[:, :, rows, columns],
            expected[kind][:, :, rows, columns],
            atol=1e-12,
            rtol=1e-12,
        )
        np.testing.assert_allclose(
            tensor.transpose(0, 1, 3, 2)[:, :, rows, columns],
            expected[kind].transpose(0, 1, 3, 2)[:, :, rows, columns],
            atol=1e-12,
            rtol=1e-12,
        )
    assert metrics["uploaded_blocks"] == 4


def test_cache_limit_does_not_change_scientific_block_values(example):
    pcm, surface, metrics = example
    expected = pcm.get_d2D_d2S(surface, with_S=True, with_D=True)[0]
    cache = CurvatureBlocks(
        surface, pcm.get_d2D_d2S, np.asarray, metrics, cache_bytes=4000
    )
    tensor = OffDiagonalTensor(cache, 0)
    for i in (0, 5, 10, 0):
        block = tensor[:, :, i : i + 5, 2:7]
        np.testing.assert_allclose(
            block, expected[:, :, i : i + 5, 2:7], atol=1e-12, rtol=1e-12
        )
        assert cache.bytes <= 4000


def test_invalid_layouts_and_unbounded_slices_fail_closed(example):
    pcm, surface, metrics = example
    cache = CurvatureBlocks(
        surface, pcm.get_d2D_d2S, np.asarray, metrics, cache_bytes=1
    )
    with pytest.raises(MemoryError):
        OffDiagonalTensor(cache, 0)[:, :, :5, :5]
    with pytest.raises(ValueError):
        OffDiagonalTensor(cache, 0).transpose(1, 0, 2, 3)
