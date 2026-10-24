"""SCF logging preserves solver settings, callbacks and restoration on failure."""
from pathlib import Path
from types import SimpleNamespace
import importlib.util
import json

import numpy as np
import pytest

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('checkpoint_scf_worker',ROOT/'scripts/fleet/checkpoint_scf_worker.py')
MODULE=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(MODULE)


def test_checkpoint_callback_keeps_solver_and_existing_callback(tmp_path,monkeypatch):
    from redox.qm import dft
    called=[]
    molecule=SimpleNamespace(charge=0,spin=0,basis='def2-svp',atom_coords=lambda:np.zeros((1,3)))
    mf=SimpleNamespace(mol=molecule,xc='r2scan',disp='d4',with_solvent=SimpleNamespace(solvent='acetonitrile'),
        conv_tol=1e-12,callback=lambda env:called.append(env['e_tot']))
    def factory(*args,**kwargs):return mf
    monkeypatch.setattr(dft,'_build_smd_mf',factory)
    with MODULE.checkpoints(tmp_path/'checkpoints'):
        assert dft._build_smd_mf(None,'r2scan') is mf
        mf.callback(dict(cycle=2,e_tot=-10.))
        assert mf.conv_tol==1e-12 and called==[-10.]
        assert Path(mf.chkfile).parent==tmp_path/'checkpoints/SCF_0000'
    assert dft._build_smd_mf is factory
    record=json.loads((tmp_path/'checkpoints/SCF_0000/progress.json').read_text())
    assert record['E_Ha']==-10 and not record['qualified_result']


def test_factory_is_restored_after_an_interrupted_worker(tmp_path,monkeypatch):
    from redox.qm import dft
    original=dft._build_smd_mf
    with pytest.raises(RuntimeError):
        with MODULE.checkpoints(tmp_path/'checkpoints'):
            raise RuntimeError('interrupted before calculation')
    assert dft._build_smd_mf is original
