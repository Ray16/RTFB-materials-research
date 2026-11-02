"""Actual determinant, SCF restart and callback retention under compact storage."""
import gc
import weakref

import h5py
import numpy as np
import pytest
from pyscf import gto, scf

from redox.qm.occupied_checkpoints import attach, density, read, review, save


def molecule():
    return gto.M(atom='O 0 0 0; H 0 0 .97; H .94 0 -.24', basis='6-31g', verbose=0)


def test_actual_solver_determinant_density_and_restart_are_unchanged(tmp_path):
    mol = molecule(); reference = scf.RHF(mol); reference.conv_tol = 1e-11
    energy = reference.kernel()
    actual = attach(scf.RHF(mol), tmp_path / 'solves'); actual.conv_tol = reference.conv_tol
    assert abs(actual.kernel() - energy) < 1e-12
    path = tmp_path / 'solves/0000/occupied.h5'
    np.testing.assert_allclose(density(path), actual.make_rdm1(), atol=1e-14, rtol=0)
    saved_mol, saved = read(path)
    np.testing.assert_array_equal(saved['mo_coeff_occupied'], actual.mo_coeff[:, actual.mo_occ > 0])
    np.testing.assert_array_equal(saved['mo_energy'], actual.mo_energy)
    assert saved_mol.dumps() == mol.dumps()
    assert review(path)['converged']
    assert abs(scf.RHF(saved_mol).kernel(dm0=density(path)) - energy) < 1e-11
    with h5py.File(path) as f:
        assert 'scf/mo_coeff' not in f and f['scf/mo_coeff_occupied'].shape[1] < mol.nao_nr()
    actual.kernel(dm0=density(path))
    assert (tmp_path / 'solves/0001/occupied.h5').exists()


def test_interruption_keeps_the_last_actual_iteration_density(tmp_path):
    target = scf.RHF(molecule())
    def stop(env):
        if env['cycle'] == 2:
            raise KeyboardInterrupt('deliberate interruption')
    target.callback = stop; attach(target, tmp_path / 'solves')
    with pytest.raises(KeyboardInterrupt):
        target.kernel()
    path = tmp_path / 'solves/0000/occupied.h5'
    _, saved = read(path)
    assert saved['cycle'] == 1 and not saved['converged']
    assert np.isfinite(density(path)).all()
    assert review(path)['density_restart_supported']


def test_no_completed_mean_field_is_kept_alive(tmp_path):
    target = attach(scf.RHF(molecule()), tmp_path / 'solves')
    pointer = weakref.ref(target); del target; gc.collect()
    assert pointer() is None


def test_wrong_electron_count_and_corrupted_occupied_space_are_rejected(tmp_path):
    mf = scf.RHF(molecule()); mf.kernel(); path = tmp_path / 'state.h5'
    bad = mf.mo_occ.copy(); bad[0] = 0
    with pytest.raises(ValueError, match='occupations'):
        save(path, mf.mol, mf.mo_coeff, bad, mf.mo_energy, mf.e_tot, converged=True, cycle=-1)
    save(path, mf.mol, mf.mo_coeff, mf.mo_occ, mf.mo_energy, mf.e_tot, converged=True, cycle=-1)
    with h5py.File(path, 'r+') as f:
        f['scf/mo_coeff_occupied'][0, 0] += .01
    with pytest.raises(ValueError, match='AO-metric'):
        review(path)


def test_multiple_guesses_keep_each_final_determinant_and_restore_solver(tmp_path):
    from redox.qm import dft
    from redox.qm.occupied_checkpoints import scf_checkpoints
    before = dft._lowest_scf, dft._kernel_robust
    with scf_checkpoints(tmp_path / 'SCF'):
        energy, converged, _, guesses = dft._lowest_scf(lambda: scf.RHF(molecule()), ('minao', 'atom'))
    assert (dft._lowest_scf, dft._kernel_robust) == before
    paths = sorted(tmp_path.glob('SCF/*/final/occupied.h5'))
    assert len(paths) == len(guesses) == 2 and converged
    assert abs(min(review(p)['energy_Ha'] for p in paths) - energy) < 1e-12
