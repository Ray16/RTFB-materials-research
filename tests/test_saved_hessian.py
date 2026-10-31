"""A complete checkpoint is usable only at its fresh matched Cartesian state."""
import hashlib
import json
from types import SimpleNamespace

import numpy as np
import pytest

from redox.qm.saved_hessian import builder


@pytest.fixture
def source(tmp_path):
    coordinates=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
    row=dict(SCF_converged=True,E_Ha=-10.,coordinates_Bohr=coordinates.tolist(),gradient_Ha_Bohr=np.zeros((3,3)).tolist())
    case=tmp_path/'case';(case/'points/0000').mkdir(parents=True)
    (case/'points/0000/evaluation.json').write_text(json.dumps(row))
    old=tmp_path/'source.json';old.write_text(json.dumps(row));matrix=tmp_path/'Hessian.npy'
    dense=np.arange(81,dtype=np.float64).reshape(9,9)/100;np.save(matrix,dense.T@dense+np.eye(9))
    plan=dict(symbols=['O','H','H'],saved_full_Hessian=dict(center=str(old),matrix=str(matrix),
        source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [old,matrix]}))
    mol=SimpleNamespace(charge=0,spin=0,basis='def2-svp',atom_coords=lambda:coordinates.copy())
    mf=SimpleNamespace(mol=mol,xc='r2scan',Hessian=lambda:SimpleNamespace(kernel=lambda:'original'))
    original=lambda *a,**k:mf;module=SimpleNamespace(_build_smd_mf=original)
    return module,original,mf,plan,case,matrix


def test_exact_source_matrix_reaches_unchanged_worker_tensor_order(source):
    module,original,mf,plan,case,matrix=source
    with builder(module,plan,case) as calls:
        actual=module._build_smd_mf().Hessian().kernel()
    np.testing.assert_array_equal(actual.transpose(0,2,1,3).reshape(9,9),np.load(matrix))
    assert len(calls)==1 and not calls[0]['full_response_recomputed']
    assert module._build_smd_mf is original


@pytest.mark.parametrize('change',['geometry','energy','force','basis','source_bytes'])
def test_wrong_state_or_checkpoint_never_borrows_the_matrix(source,change):
    module,original,mf,plan,case,matrix=source
    p=case/'points/0000/evaluation.json';row=json.loads(p.read_text())
    if change=='geometry':row['coordinates_Bohr'][0][0]=.01
    if change=='energy':row['E_Ha']+=.01
    if change=='force':row['gradient_Ha_Bohr'][0][0]=.001
    if change=='basis':mf.mol.basis='def2-svpd'
    if change=='source_bytes':np.save(matrix,np.eye(9)*2)
    p.write_text(json.dumps(row))
    with pytest.raises(ValueError):
        with builder(module,plan,case):module._build_smd_mf().Hessian().kernel()


def test_solvent_factory_shallow_copy_does_not_inherit_instance_override(source):
    module,original,mf,plan,case,matrix=source
    class Vacuum:
        def Hessian(self):return SimpleNamespace(kernel=lambda:'vacuum')
    class Solvent(Vacuum):
        def undo_solvent(self):
            copy=Vacuum();copy.__dict__.update(self.__dict__);return copy
        def Hessian(self):return self.undo_solvent().Hessian()
    target=Solvent();target.__dict__.update({k:v for k,v in mf.__dict__.items() if k!='Hessian'})
    module._build_smd_mf=lambda *a,**k:target
    with builder(module,plan,case):actual=module._build_smd_mf().Hessian().kernel()
    np.testing.assert_array_equal(actual.transpose(0,2,1,3).reshape(9,9),np.load(matrix))


def test_saved_validation_uses_bounded_operators_without_claiming_a_new_response_solve():
    import importlib.util
    from pathlib import Path
    root=Path(__file__).resolve().parents[1];p=root/'scripts/validation/run_screening_polymer_hessian.py'
    spec=importlib.util.spec_from_file_location('saved_phase_operator_scope',p);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    assert module.workspace_file({},root).name=='residual_response_solver.py'
    assert module.workspace_file(dict(Hessian_validation_continuation=True),root).name=='df_occupied_workspace.py'
