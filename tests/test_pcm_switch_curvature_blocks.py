"""Grid tiling retains full-atom switch curvature and bounded contractions."""
import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('host_switch_curvature',ROOT/'scripts/pipeline/pcm_switch_curvature_blocks.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)


@pytest.fixture
def surface():
    from pyscf import gto
    from pyscf.solvent import pcm
    mol=gto.M(atom='O 0 0 0;H 0 0 1;H 0 1 0',basis='sto-3g',verbose=0)
    return pcm.gen_surface(mol,ng=38)


def test_all_atom_curvature_and_area_contractions_match_dense_CPU_primitive(surface):
    from pyscf.solvent.hessian import pcm
    expected=pcm.get_d2F_d2A(surface)
    metrics=dict(host_tensor_bytes=0,primitive_calls=0,maximum_tile_tensor_bytes=0,
        maximum_uploaded_bytes=0,uploads=0)
    actual=m.assemble(surface,pcm.get_d2F_d2A,np,metrics,tile_points=7)
    vector=np.linspace(-.3,.8,len(surface['grid_coords']))
    for tensor,dense in zip(actual,expected):
        np.testing.assert_allclose(tensor.array,dense,atol=1e-12,rtol=1e-12)
        np.testing.assert_allclose(tensor@vector,dense@vector,atol=1e-12,rtol=1e-12)
        np.testing.assert_allclose(tensor.flatten(),dense.flatten(),atol=1e-12,rtol=1e-12)
        np.testing.assert_allclose(np.sum(tensor.array,axis=0),0.,atol=1e-12)
    assert metrics['primitive_calls']>1
    assert metrics['maximum_tile_tensor_bytes']<=2*7*3*3*9*8


def test_surface_tiles_keep_all_atoms_and_original_ownership(surface):
    local=m.local_surface(surface,5,20)
    assert local['atom_coords'] is surface['atom_coords']
    assert sum(b-a for a,b in local['gslice_by_atom'])==15
    assert len(local['gslice_by_atom'])==len(surface['atom_coords'])
    np.testing.assert_array_equal(local['grid_coords'],surface['grid_coords'][5:20])


def test_oversized_host_storage_and_dense_uploads_are_refused(surface):
    from pyscf.solvent.hessian import pcm
    metrics=dict(host_tensor_bytes=0,primitive_calls=0,maximum_tile_tensor_bytes=0,
        maximum_uploaded_bytes=0,uploads=0)
    with pytest.raises(MemoryError):
        m.assemble(surface,pcm.get_d2F_d2A,np,metrics,maximum_host_bytes=1)
    data=np.zeros((10,3,3,3,3))
    tensor=m.HostSwitchTensor(data,np,metrics,maximum_upload_bytes=1)
    with pytest.raises(MemoryError):tensor.flatten()
