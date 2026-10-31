"""Byte-exact orbital/derivative binding and immutable component reuse."""
import hashlib
import json
import sys
from types import SimpleNamespace

import h5py
import numpy as np
import pytest

from redox.qm.hessian_restart import read_frame,reuse_components,array_digest,response_ceiling_gib,capture_solution_h1,read_continuation_frame


@pytest.fixture
def frame_case(tmp_path):
    root=tmp_path/'case';execution=root/'execution';point=root/'points/0000';scf=root/'SCF_evaluations/0000/0000'
    for p in (execution,point,scf):p.mkdir(parents=True)
    arrays=dict(mo_coeff=np.eye(2),mo_occ=np.array([2.,0.]),mo_energy=np.array([-1.,1.]))
    with h5py.File(scf/'SCF.h5','w') as f:
        f['scf/e_tot']=-1.
        for k,v in arrays.items():f['scf/'+k]=v
    identity=dict(coordinates_Bohr=[[0.,0.,0.]],basis="'sto-3g'",charge=0,mult=1,xc='r2scan',SCF_E_Ha=-1.,SCF_converged=True)
    identity.update({k:dict(shape=list(v.shape),dtype='float64',sha256=array_digest(v)) for k,v in arrays.items()})
    for name,value in [('partial_hess_elec',np.eye(3).reshape(1,1,3,3)),('make_h1',np.ones((1,3,2,1)))]:
        path=execution/f'rks_{name}_0.npy';np.save(path,value)
        row=dict(identity,stage='rks_'+name,arrays={str(path):hashlib.sha256(path.read_bytes()).hexdigest()})
        (execution/f'rks_{name}_state_identity.json').write_text(json.dumps(row))
    (point/'evaluation.json').write_text(json.dumps(dict(SCF_converged=True,coordinates_Bohr=identity['coordinates_Bohr'],E_Ha=-1.)))
    (scf/'settings.json').write_text(json.dumps(dict(coordinates_Bohr=identity['coordinates_Bohr'],charge=0,mult=1,xc='r2scan')))
    np.savez(root/'grids.npz',coords=np.zeros((2,3)),weights=np.ones(2))
    return root


def test_exact_orbitals_grid_and_completed_components_share_one_frame(frame_case):
    frame,arrays=read_frame(frame_case)
    assert frame['natom']==1 and frame['nao']==2 and frame['nocc']==1
    assert arrays['h1'].shape==(1,3,2,1) and not arrays['h1'].flags.writeable


def test_response_cap_accounts_for_the_installed_minimum_aligned_block():
    assert response_ceiling_gib(1037,206)[0]==3
    cap,minimum=response_ceiling_gib(1193,238)
    assert cap==4 and 3*2**30<minimum<=4*2**30
    with pytest.raises(ValueError):response_ceiling_gib(10000,5000)


@pytest.mark.parametrize('change',['orbitals','frame','H1','grids'])
def test_changed_partial_orbital_or_grid_data_cannot_resume(frame_case,change):
    if change=='orbitals':
        with h5py.File(frame_case/'SCF_evaluations/0000/0000/SCF.h5','r+') as f:f['scf/mo_coeff'][...]=0.
    if change=='frame':
        p=frame_case/'execution/rks_make_h1_state_identity.json';r=json.loads(p.read_text());r['coordinates_Bohr'][0][0]=.1;p.write_text(json.dumps(r))
    if change=='H1':np.save(frame_case/'execution/rks_make_h1_0.npy',np.zeros((1,3,2,1)))
    if change=='grids':np.savez(frame_case/'grids.npz',coords=np.full((2,3),np.nan),weights=np.ones(2))
    with pytest.raises(ValueError):read_frame(frame_case)


def test_reuse_returns_independent_copies_and_restores_only_the_target_object(tmp_path,monkeypatch):
    monkeypatch.setitem(sys.modules,'cupy',SimpleNamespace(asarray=np.asarray))
    class Hessian:
        def partial_hess_elec(self):return 'original partial'
        def make_h1(self):return np.zeros((2,2))
    class WithSolvent:
        def undo_solvent(self):pass
        def make_h1(self):return self.__class__.__bases__[1].make_h1(self)+2
    class SolutionHessian(WithSolvent,Hessian):pass
    target=SolutionHessian();other=SolutionHessian();a=np.ones((2,2));a.flags.writeable=False
    with reuse_components(target,dict(partial=a,h1=a),tmp_path/'reuse'):
        partial=target.partial_hess_elec();partial[0,0]=42
        assert target.make_h1()[0,0]==3 and a[0,0]==1
        assert other.partial_hess_elec()=='original partial'
        assert other.make_h1()[0,0]==2
    assert target.partial_hess_elec()=='original partial' and target.make_h1()[0,0]==2


def test_captured_solution_checkpoint_keeps_the_solvent_addition(tmp_path,monkeypatch):
    monkeypatch.setitem(sys.modules,'cupy',SimpleNamespace(asarray=np.asarray))
    class Base:
        def __init__(self):self.base=SimpleNamespace(with_solvent=SimpleNamespace(equilibrium_solvation=True))
        def partial_hess_elec(self):return np.zeros((2,2))
        def make_h1(self):return np.zeros((2,2))
    class Solvent:
        def undo_solvent(self):pass
        def make_h1(self):return self.__class__.__bases__[1].make_h1(self)+2.
    class Solution(Solvent,Base):pass
    target=Solution();arrays=dict(partial=np.ones((2,2)),h1=np.ones((2,2)))
    with reuse_components(target,arrays,tmp_path/'reuse'),capture_solution_h1(target,tmp_path/'full',dict(frame='same')):
        target.partial_hess_elec();value=target.make_h1()
    np.testing.assert_array_equal(value,np.full((2,2),3.))
    np.testing.assert_array_equal(np.load(tmp_path/'full/full_H1.npy'),value)
    assert json.loads((tmp_path/'full/identity.json').read_text())['scope']=='full SMD H1 including PCM contribution'
    target.base.with_solvent.equilibrium_solvation=False
    with pytest.raises(ValueError,match='equilibrium'):
        with capture_solution_h1(target,tmp_path/'frozen',dict(frame='same')):target.make_h1()


def test_continuation_resolves_original_SCF_but_checks_its_own_center_and_grids(frame_case,tmp_path):
    import shutil
    frame,_=read_frame(frame_case);warm=tmp_path/'warm';(warm/'points/0000').mkdir(parents=True)
    (warm/'plan.json').write_text(json.dumps(dict(resume_components=dict(source_case=str(frame_case),frame=frame))))
    shutil.copyfile(frame_case/'points/0000/evaluation.json',warm/'points/0000/evaluation.json')
    shutil.copyfile(frame_case/'grids.npz',warm/'grids.npz')
    restored,_,source=read_continuation_frame(warm)
    assert restored==frame and source==frame_case
    np.savez(warm/'grids.npz',coords=np.zeros((2,3)),weights=np.full(2,2.))
    with pytest.raises(ValueError,match='grids'):read_continuation_frame(warm)


def response_source(frame_case):
    from redox.qm.response_diagnostic import save_array
    from redox.qm.response_resume import digest
    frame,_=read_frame(frame_case);center=frame_case/'points/0000/evaluation.json'
    (frame_case/'plan.json').write_text(json.dumps(dict(role='screening_polymer_Hessian_control',screening_variant='corrected',full_polymer_Hessian_control=True,original_namespace='/original.json')))
    identity=dict(center_evaluation=str(center),center_evaluation_sha256=digest(center),grids_sha256=digest(frame_case/'grids.npz'),source_namespace='/original.json',screening_variant='corrected')
    full=frame_case/'solution_H1';full.mkdir();payload=save_array(full/'full_H1.npy',np.ones((1,3,2,1)))
    (full/'identity.json').write_text(json.dumps(dict(payload,scope='full SMD H1 including PCM contribution',identity=identity)))
    blocks=frame_case/'response_blocks';blocks.mkdir();payload=save_array(blocks/'block_000.npy',np.ones((3,2)))
    (blocks/'block_000.json').write_text(json.dumps(dict(checkpoint=payload,identity=identity,rhs_sha256='test',original_solver_source_sha256='solver',kind='raw_Krylov_solution_before_occupied_block_reconstruction')))
    return frame


def test_only_explicit_full_solvent_H1_can_supply_raw_blocks(frame_case):
    from redox.qm.response_resume import inspect_source
    frame=response_source(frame_case);bound=inspect_source(frame_case,frame)
    assert bound['completed_blocks']==1 and bound['completed_RHS']==3
    assert bound['full_H1']['scope']=='full SMD H1 including PCM contribution'


@pytest.mark.parametrize('change',['H1_scope','H1_payload','block_identity','block_payload','center','gapped_prefix','wrong_role'])
def test_mismatched_full_response_backing_is_not_admitted(frame_case,change):
    from redox.qm.response_resume import inspect_source
    frame=response_source(frame_case)
    if change=='H1_scope':
        p=frame_case/'solution_H1/identity.json';r=json.loads(p.read_text());r['scope']='solute contribution';p.write_text(json.dumps(r))
    if change=='H1_payload':np.save(frame_case/'solution_H1/full_H1.npy',np.zeros((1,3,2,1)))
    if change=='block_identity':
        p=frame_case/'response_blocks/block_000.json';r=json.loads(p.read_text());r['identity']['grids_sha256']='changed';p.write_text(json.dumps(r))
    if change=='block_payload':np.save(frame_case/'response_blocks/block_000.npy',np.zeros((3,2)))
    if change=='center':
        p=frame_case/'points/0000/evaluation.json';r=json.loads(p.read_text());r['coordinates_Bohr'][0][0]=.1;p.write_text(json.dumps(r))
    if change=='gapped_prefix':(frame_case/'response_blocks/block_000.json').rename(frame_case/'response_blocks/block_001.json')
    if change=='wrong_role':
        p=frame_case/'plan.json';r=json.loads(p.read_text());r['role']='frozen_solvent_response_probe';p.write_text(json.dumps(r))
    with pytest.raises(ValueError):inspect_source(frame_case,frame)
