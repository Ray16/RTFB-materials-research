import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
from ase import Atoms
from ase.io import write

ROOT=Path(__file__).resolve().parents[1]


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,ROOT/path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


reviewer=load('scripts/validation/review_molecular_native.py','tested_saved_molecular_derivatives')
observer=load('scripts/pipeline/observe_molecular_native.py','tested_molecular_observer')


@pytest.fixture
def derivatives(tmp_path):
    xyz=tmp_path/'geometry.xyz';atoms=Atoms('H2',positions=[[0,0,0],[0,0,.8]]);write(xyz,atoms)
    coordinates=atoms.positions/.52917721067
    matrix=np.arange(36,dtype=float).reshape(6,6)/100;matrix=(matrix+matrix.T)/2+np.eye(6)
    base=(matrix@coordinates.reshape(-1)).reshape(2,3)
    checks=[];points=[coordinates]
    for atom,axis in [(0,1),(1,2)]:
        for step in [.0005,.001]:
            for sign in [-1,1]:
                x=coordinates.copy();x[atom,axis]+=sign*step;points.append(x)
            checks.append(dict(atom_index=atom,axis=axis,step_Bohr=step,relative_error=0.))
    for i,x in enumerate(points):
        folder=tmp_path/'SCF_evaluations'/f'{i:04d}'/'0000';folder.mkdir(parents=True)
        (folder/'settings.json').write_text(json.dumps(dict(coordinates_Bohr=x.tolist(),charge=0,mult=1)))
        (folder/'SCF_result.json').write_text(json.dumps(dict(converged=True)))
        np.save(folder/'gradient.npy',(matrix@x.reshape(-1)).reshape(2,3))
    record=dict(identity=dict(charge=0,mult=1),gradient_max=float(np.abs(base).max()),derivative_checks=checks)
    return tmp_path,xyz,record,matrix


def test_sampled_columns_reconstruct_from_actual_saved_displacements(derivatives):
    checks,paths=reviewer.finite_difference_review(*derivatives)
    assert len(checks)==4 and max(c['relative_error'] for c in checks)<1e-11
    assert len(paths)==27


@pytest.mark.parametrize('change',['gradient','missing','charge','coverage'])
def test_unreconstructable_hessian_evidence_cannot_be_admitted(derivatives,change):
    root,xyz,record,matrix=derivatives
    folder=root/'SCF_evaluations/0001/0000'
    if change=='gradient':np.save(folder/'gradient.npy',np.ones((2,3)))
    if change=='missing':(folder/'gradient.npy').unlink()
    if change=='charge':
        p=folder/'settings.json';r=json.loads(p.read_text());r['charge']=1;p.write_text(json.dumps(r))
    if change=='coverage':record['derivative_checks']=record['derivative_checks'][:-1]
    with pytest.raises(ValueError):reviewer.finite_difference_review(root,xyz,record,matrix)


def test_guard_and_runtime_closure_are_required_before_collecting(tmp_path):
    assert not observer.execution_closed(tmp_path)
    for name,record in [('residency_guard/owner.json',dict(guard_stopped=True)),('host_guard/status.json',dict(finished=1)),('tensor_runtime/status.json',dict(status='running'))]:
        p=tmp_path/name;p.parent.mkdir();p.write_text(json.dumps(record))
    assert not observer.execution_closed(tmp_path)
    (tmp_path/'tensor_runtime/status.json').write_text(json.dumps(dict(finished=1)))
    assert observer.execution_closed(tmp_path)


def test_another_namespace_cannot_reuse_the_same_state_review(tmp_path):
    namespace=tmp_path/'namespace.json';namespace.write_text('{}')
    p=tmp_path/'review.json';p.write_text(json.dumps(dict(namespace=str(namespace),namespace_sha256='wrong')))
    with pytest.raises(ValueError):observer.checked_cached(p,namespace)


def test_delayed_shared_review_is_pending_and_does_not_start_another_producer(tmp_path):
    from redox.validation.molecular_native import completed_review_record,ReviewPending
    path=tmp_path/'delayed.json'
    with pytest.raises(ReviewPending):completed_review_record(path,timeout_seconds=0)
    path.write_text('{')
    with pytest.raises(ReviewPending):completed_review_record(path,timeout_seconds=0)
    path.write_text('{"status":"finished"}')
    assert completed_review_record(path,timeout_seconds=0)=={'status':'finished'}
