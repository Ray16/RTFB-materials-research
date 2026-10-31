"""Checkpoint-frame rejection and value-preserving diagnostic instrumentation."""
import hashlib
import importlib.util
import json
import time
from pathlib import Path
from types import SimpleNamespace

import h5py
import numpy as np
import pytest

from redox.qm.response_diagnostic import read_checkpoint,trace_solver,response_memory_cap,save_array


@pytest.fixture
def frame(tmp_path):
    raw=tmp_path/'raw';execution=raw/'execution';scf=raw/'SCF_evaluations/0000';execution.mkdir(parents=True);scf.mkdir(parents=True)
    arrays=dict(mo_coeff=np.eye(2),mo_energy=np.array([-1.,1.]),mo_occ=np.array([2.,0.]))
    with h5py.File(scf/'SCF.h5','w') as f:
        f['scf/e_tot']=-1.
        for k,v in arrays.items():f['scf/'+k]=v
    matrix=execution/'rks_make_h1_0.npy';np.save(matrix,np.ones((1,3,2,1)))
    identity=dict(stage='rks_make_h1',basis="'def2-svp'",charge=0,mult=1,xc='r2scan',SCF_converged=True,SCF_E_Ha=-1.,coordinates_Bohr=[[0.,0.,0.]],arrays={str(matrix):hashlib.sha256(matrix.read_bytes()).hexdigest()})
    for k,v in arrays.items():identity[k]=dict(shape=list(v.shape),dtype='float64',sha256=hashlib.sha256(v.tobytes()).hexdigest())
    (execution/'rks_make_h1_state_identity.json').write_text(json.dumps(identity))
    (scf/'evaluation.json').write_text(json.dumps(dict(SCF_converged=True,E_Ha=-1.,coordinates_Bohr=[[0.,0.,0.]])))
    (scf/'settings.json').write_text(json.dumps(dict(charge=0,mult=1,basis='def2-svp',xc='r2scan',disp='d4',solvent='acetonitrile',grid_level=3)))
    (raw/'geometry.xyz').write_text('1\nNumeric checkpoint fixture\nHe 0 0 0\n');return raw


def test_checkpoint_retains_exact_orbital_frame_and_readonly_first_derivative(frame):
    result,arrays,h1=read_checkpoint(frame)
    assert result['energy_Ha']==-1. and np.array_equal(arrays['mo_coeff'],np.eye(2))
    assert h1.shape==(1,3,2,1) and not h1.flags.writeable


@pytest.mark.parametrize('change',['orbital','energy','matrix','geometry','solvent','spin','compressed'])
def test_changed_orbital_frame_or_numerical_model_is_rejected(frame,change):
    if change in ('orbital','energy'):
        with h5py.File(frame/'SCF_evaluations/0000/SCF.h5','r+') as f:
            key='scf/mo_coeff' if change=='orbital' else 'scf/e_tot';f[key][...]=0.
    if change=='matrix':np.save(frame/'execution/rks_make_h1_0.npy',np.zeros((1,3,2,1)))
    if change=='compressed':
        with h5py.File(frame/'SCF_evaluations/0000/SCF.h5','r+') as f:
            del f['scf/mo_coeff'];f.create_dataset('scf/mo_coeff',data=np.eye(2),compression='gzip')
    if change=='geometry':
        p=frame/'SCF_evaluations/0000/evaluation.json';r=json.loads(p.read_text());r['coordinates_Bohr'][0][0]=.1;p.write_text(json.dumps(r))
    if change=='solvent':
        p=frame/'SCF_evaluations/0000/settings.json';r=json.loads(p.read_text());r['solvent']='water';p.write_text(json.dumps(r))
    if change=='spin':
        p=frame/'execution/rks_make_h1_state_identity.json';r=json.loads(p.read_text());r['mult']=3;p.write_text(json.dumps(r))
    with pytest.raises(ValueError):read_checkpoint(frame)


def test_trace_retains_frozen_solver_outputs_and_operator_calls(tmp_path):
    path=Path(__file__).parent/'fixtures/response_solver_reference.py'
    reference=json.loads(path.with_suffix('.json').read_text())
    assert hashlib.sha256(path.read_bytes()).hexdigest()==reference['sha256']
    spec=importlib.util.spec_from_file_location('frozen_diagnostic_residual',path);frozen=importlib.util.module_from_spec(spec);spec.loader.exec_module(frozen)
    matrix=np.array([[.2,.1],[.1,.3]]);rhs=np.array([[1.,2.],[2.,-1.]])
    seen=[]
    def operator(x):seen.append(x.copy());return x@matrix
    expected=frozen.residual_krylov(operator,rhs);original_calls=[x.copy() for x in seen];seen.clear()
    module=SimpleNamespace(residual_krylov=frozen.residual_krylov)
    with trace_solver(module,tmp_path/'trace',trace_seconds=30):actual=module.residual_krylov(operator,rhs)
    assert np.array_equal(actual,expected) and len(seen)==len(original_calls)
    assert all(np.array_equal(a,b) for a,b in zip(seen,original_calls))
    assert module.residual_krylov is frozen.residual_krylov
    assert json.loads((tmp_path/'trace/trace.json').read_text())['calls'][0]['status']=='returned'


def test_failed_operator_restores_solver_and_keeps_last_event(tmp_path):
    def solver(aop,b):return aop(b)
    def operator(x):raise RuntimeError('intentional diagnostic failure')
    module=SimpleNamespace(residual_krylov=solver)
    with pytest.raises(RuntimeError):
        with trace_solver(module,tmp_path/'trace'):module.residual_krylov(operator,np.ones((2,2)))
    assert module.residual_krylov is solver
    assert json.loads((tmp_path/'trace/latest.json').read_text())['kind']=='operator_entered'


def test_prospective_batch_cap_never_increases_available_memory_and_restores_on_error():
    module=SimpleNamespace(get_avail_mem=lambda:8*2**30);original=module.get_avail_mem
    with pytest.raises(RuntimeError):
        with response_memory_cap(module,3*2**30):
            assert module.get_avail_mem()==3*2**30
            raise RuntimeError('test cleanup')
    assert module.get_avail_mem is original
    module.get_avail_mem=lambda:2**30
    with response_memory_cap(module,3*2**30):assert module.get_avail_mem()==2**30


def test_periodic_sampler_never_uses_the_GIL_free_C_watchdog(tmp_path,monkeypatch):
    import redox.qm.response_diagnostic as module
    monkeypatch.setattr(module.faulthandler,'dump_traceback_later',lambda *a,**k:pytest.fail('Unsafe C watchdog must not be started'))
    solver=SimpleNamespace(residual_krylov=lambda aop,b:aop(b))
    with module.trace_solver(solver,tmp_path/'trace',trace_seconds=.001):time.sleep(.015)
    assert (tmp_path/'trace/stacks.log').stat().st_size>0
    assert not json.loads((tmp_path/'trace/trace.json').read_text())['asynchronous_C_watchdog_used']


def test_solved_response_block_is_checkpointed_without_changing_its_value(tmp_path):
    def solver(aop,b):return aop(b)
    module=SimpleNamespace(residual_krylov=solver);rhs=np.array([[1.,2.]])
    with trace_solver(module,tmp_path/'trace',checkpoint_root=tmp_path/'blocks',checkpoint_identity={'frame':'test'}):
        value=module.residual_krylov(lambda x:x*2,rhs)
    assert np.array_equal(value,rhs*2) and np.array_equal(np.load(tmp_path/'blocks/block_000.npy'),value)
    row=json.loads((tmp_path/'blocks/block_000.json').read_text());assert row['identity']=={'frame':'test'}
    assert row['rhs_sha256']==hashlib.sha256(rhs.tobytes()).hexdigest() and not row['native_G_qualified']
    assert not list((tmp_path/'blocks').glob('*.partial'))
    with pytest.raises(FileExistsError):save_array(tmp_path/'blocks/block_000.npy',value)
    with pytest.raises(ValueError):save_array(tmp_path/'invalid.npy',np.array([np.nan]))
