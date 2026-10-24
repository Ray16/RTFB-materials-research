"""Streaming must preserve tensor axes, weighting and transposed cross blocks."""
from pathlib import Path
import importlib.util

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('streamed_pcm_math',ROOT/'scripts/pipeline/streamed_pcm_execution.py')
MODULE=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


def metrics():
    return dict(diagonal_contractions=0,offdiagonal_primitive_calls=0,
        maximum_offdiagonal_primitive_bytes=0,peak_GPU_block_cache_bytes=0,
        offdiagonal_cache_hits=0,offdiagonal_block_requests=0)


def test_diagonal_contraction_matches_dense_reference_with_signed_weights():
    from redox.qm.pcm_bounded import diagonal_contraction
    rng=np.random.default_rng(61);n,points=3,11
    first=rng.normal(size=(3,points,n));second=rng.normal(size=(n,n,3,3,points))
    surface=dict(switch_fun=rng.uniform(.2,1,points),charge_exp=rng.uniform(.3,2,points))
    class Source:
        atoms=n
        def __init__(self):self.points=points
        def contract(self,kind,vector):return second@vector
    operator=MODULE.StreamedDiagonal(surface,first,MODULE.WeightedSwitch(Source(),'F'),np,metrics())
    vector=rng.normal(size=points)
    np.testing.assert_allclose(operator@vector,diagonal_contraction(surface,first,second,vector),atol=1e-12,rtol=1e-12)


def primitive(surface,**kwargs):
    x=surface['grid_coords'][:,0]
    base=(100*x[:,None]+x[None,:])[None,None,:,:]
    # An intentionally nonsymmetric tensor exposes incorrect D transposes.
    D=np.broadcast_to(base,(3,3,len(x),len(x))).copy()
    return D,2*D


def test_cross_block_indexing_and_surface_transpose_preserve_nonsymmetry():
    x=np.arange(6.)
    surface=dict(grid_coords=np.column_stack((x,x,x)),charge_exp=x+1,norm_vec=np.ones((6,3)))
    source=MODULE.GPUCurvatureBlocks(surface,primitive,np,metrics(),cache_bytes=4096)
    D,S=primitive(surface)
    operator=MODULE.GPUOffDiagonal(source,0)
    np.testing.assert_array_equal(operator[:,:,1:4,3:6],D[:,:,1:4,3:6])
    np.testing.assert_array_equal(operator.transpose(0,1,3,2)[:,:,1:4,3:6],D.transpose(0,1,3,2)[:,:,1:4,3:6])
    np.testing.assert_array_equal(MODULE.GPUOffDiagonal(source,1)[:,:,2:5,2:5],S[:,:,2:5,2:5])
    assert operator[:,:,0:0,2:4].shape==(3,3,0,2)


def test_unsupported_dense_union_is_refused_before_allocation():
    n=2000;surface=dict(charge_exp=np.ones(n),grid_coords=np.zeros((n,3)),norm_vec=np.zeros((n,3)))
    def forbidden(*args,**kwargs):raise AssertionError('A dense primitive must not run')
    source=MODULE.GPUCurvatureBlocks(surface,forbidden,np,metrics())
    with pytest.raises(MemoryError):source.block(np.arange(n),np.arange(n))


def test_switch_transformation_guards_the_installed_template():
    with pytest.raises(ValueError):MODULE.vectorized_switch_source('def other(): pass')
