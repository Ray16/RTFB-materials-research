"""Actual SCF/gradient/Hessian equality and failed-optimizer return retention."""
import gc
import json
import weakref

import numpy as np
from pyscf import gto,scf

from redox.qm.scf_checkpoints import attach,recording_optimizer


def molecule():
    return gto.M(atom='H 0 0 0; H 0 0 .75',basis='sto-3g',verbose=0)


def test_durable_checkpoint_wrappers_preserve_actual_electronic_derivatives(tmp_path):
    mol=molecule();reference=scf.RHF(mol);reference.conv_tol=1e-12
    energy=reference.kernel();gradient=reference.nuc_grad_method().kernel();hessian=reference.Hessian().kernel()
    target=scf.RHF(mol);target.conv_tol=1e-12;before=target.conv_tol
    assert attach(target,tmp_path/'solve',dict(pid=123)) is target
    assert target.conv_tol==before
    assert abs(target.kernel()-energy)<1e-13
    np.testing.assert_allclose(target.nuc_grad_method().kernel(),gradient,atol=1e-13,rtol=0)
    np.testing.assert_allclose(target.Hessian().kernel(),hessian,atol=1e-13,rtol=0)
    assert (tmp_path/'solve/0000/SCF.h5').is_file()
    np.testing.assert_allclose(np.load(tmp_path/'solve/0000/density.npy'),target.make_rdm1(),atol=0,rtol=0)
    assert json.loads((tmp_path/'solve/0000/SCF_result.json').read_text())['converged']


def test_checkpoint_wrappers_do_not_retain_a_completed_mean_field(tmp_path):
    target=attach(scf.RHF(molecule()),tmp_path/'solve',{})
    reference=weakref.ref(target);del target;gc.collect()
    assert reference() is None


def test_orbitals_only_checkpoint_reconstructs_the_same_density(tmp_path):
    from pyscf import lib
    target=attach(scf.RHF(molecule()),tmp_path/'solve',{},save_density=False)
    energy=target.kernel()
    assert not (tmp_path/'solve/0000/density.npy').exists()
    path=str(tmp_path/'solve/0000/SCF.h5')
    orbitals=lib.chkfile.load(path,'scf/mo_coeff');occupations=lib.chkfile.load(path,'scf/mo_occ')
    rebuilt=(orbitals*occupations)@orbitals.T
    np.testing.assert_allclose(rebuilt,target.make_rdm1(),atol=1e-13,rtol=0)
    assert abs(lib.chkfile.load(path,'scf/e_tot')-energy)<1e-13


def test_nonconverged_optimizer_coordinates_survive_before_caller_rejection(tmp_path):
    mol=molecule();result=(False,mol)
    run=recording_optimizer(lambda:result,tmp_path,{})
    assert run() is result
    record=json.loads((tmp_path/'0000.json').read_text())
    assert record['converged'] is False and record['qualified'] is False
    np.testing.assert_allclose(record['coordinates_A'],mol.atom_coords(unit='Angstrom'),atol=0,rtol=0)


def test_checkpointed_geometry_scanner_evaluates_the_updated_geometry(tmp_path):
    mol=molecule();mf=attach(scf.RHF(mol),tmp_path/'solve',{})
    scanner=mf.nuc_grad_method().as_scanner()
    moved=mol.copy().set_geom_('H 0 0 0; H 0 0 1.05')
    energy,gradient=scanner(moved)
    reference=scf.RHF(moved);expected=reference.kernel()
    assert abs(energy-expected)<1e-12
    np.testing.assert_allclose(gradient,reference.nuc_grad_method().kernel(),atol=1e-12,rtol=0)
    settings=json.loads((tmp_path/'solve/0000/settings.json').read_text())
    np.testing.assert_allclose(settings['coordinates_Bohr'],moved.atom_coords(),atol=0,rtol=0)


def test_unrestricted_checkpoint_wrappers_preserve_energy_and_derivatives(tmp_path):
    mol=gto.M(atom='Li 0 0 0; H 0 0 1.6',charge=1,spin=1,basis='sto-3g',verbose=0)
    reference=scf.UHF(mol);reference.conv_tol=1e-12
    energy=reference.kernel();gradient=reference.nuc_grad_method().kernel();hessian=reference.Hessian().kernel()
    target=scf.UHF(mol);target.conv_tol=1e-12;attach(target,tmp_path/'solve',{})
    assert abs(target.kernel()-energy)<1e-13
    np.testing.assert_allclose(target.nuc_grad_method().kernel(),gradient,atol=1e-12,rtol=0)
    np.testing.assert_allclose(target.Hessian().kernel(),hessian,atol=1e-12,rtol=0)
    assert np.load(tmp_path/'solve/0000/density.npy').shape[0]==2


def test_unrestricted_geometry_scanner_preserves_the_new_geometry(tmp_path):
    mol=gto.M(atom='Li 0 0 0; H 0 0 1.6',charge=1,spin=1,basis='sto-3g',verbose=0)
    target=attach(scf.UHF(mol),tmp_path/'solve',{}).nuc_grad_method().as_scanner()
    moved=mol.copy().set_geom_('Li 0 0 0; H 0 0 1.8')
    energy,gradient=target(moved)
    reference=scf.UHF(moved);expected=reference.kernel()
    assert abs(energy-expected)<1e-12
    np.testing.assert_allclose(gradient,reference.nuc_grad_method().kernel(),atol=1e-12,rtol=0)
    settings=json.loads((tmp_path/'solve/0000/settings.json').read_text())
    np.testing.assert_allclose(settings['coordinates_Bohr'],moved.atom_coords(),atol=0,rtol=0)
