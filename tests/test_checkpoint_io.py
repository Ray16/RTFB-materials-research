import errno
import json
from pathlib import Path
import pytest
from redox.workflows.checkpoint_io import retry_checkpoint


def test_only_transient_locks_are_retried_with_same_arguments(tmp_path):
    calls=[];sleeps=[]
    def operation(*args,**kwargs):
        calls.append((args,kwargs))
        if len(calls)<3:raise BlockingIOError(errno.EAGAIN,'unable to lock checkpoint')
        return 'complete'
    result=retry_checkpoint(operation,('molecule','exact.h5'),{'state':'unchanged'},tmp_path/'retry.jsonl',sleep=sleeps.append)
    assert result=='complete' and len(calls)==3 and calls[0]==calls[1]==calls[2]
    assert sleeps==[.1,.2]


def test_other_checkpoint_errors_are_preserved_immediately(tmp_path):
    def invalid(*args,**kwargs):raise ValueError('bad checkpoint')
    with pytest.raises(ValueError):retry_checkpoint(invalid,(),{},tmp_path/'retry.jsonl',sleep=lambda s:pytest.fail('Unexpected retry'))
    assert not (tmp_path/'retry.jsonl').exists()


def test_persistent_lock_has_bounded_attempts(tmp_path):
    def locked(*args,**kwargs):raise BlockingIOError(errno.EAGAIN,'checkpoint lock busy')
    with pytest.raises(BlockingIOError):retry_checkpoint(locked,('m','x.h5'),{},tmp_path/'retry.jsonl',sleep=lambda s:None,delays=(0,0))
    rows=[json.loads(s) for s in (tmp_path/'retry.jsonl').read_text().splitlines()]
    assert len(rows)==3 and rows[-1]['action']=='exhausted_preserved'


def test_checkpoint_context_preserves_cpu_hdf5_contents_and_restores_writer(tmp_path):
    import numpy as np
    from pyscf import gto
    from pyscf.scf import chkfile
    from redox.workflows.checkpoint_io import checkpoint_lock_retries
    mol=gto.M(atom='He 0 0 0',basis='sto-3g',verbose=0)
    path=tmp_path/'SCF.h5';original=chkfile.dump_scf
    with checkpoint_lock_retries(tmp_path/'retry.jsonl'):
        chkfile.dump_scf(mol,str(path),-2.,np.array([-.5]),np.array([[1.]]),np.array([2.]))
    assert chkfile.dump_scf is original
    saved,state=chkfile.load_scf(str(path))
    assert saved.nelectron==2 and float(state['e_tot'])==-2.
    np.testing.assert_array_equal(state['mo_coeff'],[[1.]])
