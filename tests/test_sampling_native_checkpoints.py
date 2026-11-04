"""Actual CPU DFT/Hessian equality and exact center/CPHF-grid retention."""
import json
from types import SimpleNamespace
import numpy as np
import pytest
from pyscf import dft,gto

from redox.qm import scf_checkpoints
from redox.qm.sampling_native_checkpoints import recording_native


def make_mf():
    mol=gto.M(atom='H 0 0 0; H 0 0 .75',basis='sto-3g',verbose=0)
    mf=dft.RKS(mol);mf.xc='lda,vwn';mf.grids.level=0;mf.conv_tol=1e-11
    return mf


def test_recording_preserves_actual_DFT_hessian_and_center_grids(tmp_path):
    reference=make_mf();energy=reference.kernel();expected=reference.Hessian().kernel()
    api=SimpleNamespace(_build_smd_mf=make_mf);before=api._build_smd_mf
    with recording_native(api,scf_checkpoints,tmp_path,{'pid':1}):
        target=api._build_smd_mf();assert abs(target.kernel()-energy)<1e-12
        target.cphf_grids=SimpleNamespace(coords=target.grids.coords.copy(),weights=target.grids.weights.copy())
        actual=target.Hessian().kernel()
    assert api._build_smd_mf is before
    np.testing.assert_allclose(actual,expected,atol=1e-11,rtol=0)
    frame=json.loads((tmp_path/'restart_frame.json').read_text())
    assert frame['CPHF_grids_matched'] and not frame['native_G_qualified']
    with np.load(tmp_path/'grids.npz') as saved:
        np.testing.assert_array_equal(saved['coords'],target.grids.coords)
        np.testing.assert_array_equal(saved['weights'],target.grids.weights)
    assert (tmp_path/'SCF_frames/0000/0000/SCF.h5').exists()


def test_mismatched_response_grids_are_rejected_before_hessian(tmp_path):
    api=SimpleNamespace(_build_smd_mf=make_mf)
    with recording_native(api,scf_checkpoints,tmp_path,{'pid':1}):
        target=api._build_smd_mf();target.kernel()
        target.cphf_grids=SimpleNamespace(coords=target.grids.coords.copy(),weights=target.grids.weights.copy())
        target.cphf_grids.weights[0]+=.1
        with pytest.raises(ValueError,match='matched finite SCF/CPHF grids'):target.Hessian()
    assert not (tmp_path/'restart_frame.json').exists()
