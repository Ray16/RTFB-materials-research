"""Physical identities for the diagnostic frozen-slow energy functional."""
import numpy as np
import pytest

from redox.qm.nonequilibrium_pcm import (
    attach_frozen_slow, frozen_slow_terms, independent_hcore_reference,
    pcm_object, surface_potential,
)


@pytest.fixture
def molecule():
    from pyscf import gto
    return gto.M(atom="O 0 0 0; H 0 0 .96; H .92 0 -.24",
                 basis="sto-3g", verbose=0, max_memory=1000)


def test_equilibrium_recovery_and_variational_derivative(molecule):
    from pyscf import scf
    mf = scf.RHF(molecule)
    dm = mf.get_init_guess()
    static = pcm_object(molecule, 37.5, order=11)
    eq_e, eq_v = static._get_vind(dm)
    noneq = attach_frozen_slow(scf.RHF(molecule), dm, 37.5, 1.806, order=11)
    ne_e, ne_v = noneq.with_solvent._get_vind(dm)
    assert abs(ne_e - eq_e) < 1e-11
    np.testing.assert_allclose(ne_v, eq_v, atol=1e-11)
    direction = np.diag(np.linspace(-.2, .2, dm.shape[0]))
    step = 1e-5
    plus = noneq.with_solvent._get_vind(dm + step * direction)[0]
    minus = noneq.with_solvent._get_vind(dm - step * direction)[0]
    expected = np.einsum("ij,ji", ne_v, direction)
    assert abs((plus - minus) / (2 * step) - expected) < 1e-9


def test_frozen_density_penalty_equals_pekar_contraction(molecule):
    from pyscf import scf
    dm0 = scf.RHF(molecule).get_init_guess()
    dm1 = dm0 * .91
    static = pcm_object(molecule, 37.5, order=11)
    optical = pcm_object(molecule, 1.806, order=11)
    ne = attach_frozen_slow(scf.RHF(molecule), dm0, 37.5, 1.806, order=11)
    penalty = ne.with_solvent._get_vind(dm1)[0] - static._get_vind(dm1)[0]
    delta = surface_potential(static, dm1) - surface_potential(static, dm0)
    # Independently contract the explicitly assembled symmetrized operators.
    def q(cm, v):
        K, R = (cm._intermediates[k] for k in ("K", "R"))
        return .5 * (np.linalg.solve(K, R @ v)
                     + R.T @ np.linalg.solve(K.T, v))
    reference = .5 * delta @ (q(optical, delta) - q(static, delta))
    assert penalty > 0
    assert abs(penalty - reference) < 1e-11


def test_no_slow_response_limit(molecule):
    from pyscf import scf
    dm = scf.RHF(molecule).get_init_guess()
    optical, q, v, self_work = frozen_slow_terms(molecule, dm, 2., 2., order=11)
    np.testing.assert_array_equal(q, 0.)
    np.testing.assert_array_equal(v, 0.)
    assert self_work == 0.


def test_independent_hcore_energy_assembly(molecule):
    from pyscf import scf
    dm0 = scf.RHF(molecule).get_init_guess()
    dm1 = dm0 * .9
    ne = attach_frozen_slow(scf.RHF(molecule), dm0, 37.5, 1.806, order=11)
    independent, constant = independent_hcore_reference(
        scf.RHF(molecule), dm0, 37.5, 1.806, order=11)
    assert abs(ne.energy_tot(dm1) - independent.energy_tot(dm1) - constant) < 1e-10


def test_invalid_dielectrics_and_basis_fail_closed(molecule):
    from pyscf import scf
    dm = scf.RHF(molecule).get_init_guess()
    with pytest.raises(ValueError):
        frozen_slow_terms(molecule, dm, 1.8, 37.5, order=11)
    with pytest.raises(ValueError):
        frozen_slow_terms(molecule, np.eye(2), 37.5, 1.8, order=11)


def test_dense_reference_preserves_actual_discretized_grid(molecule):
    import importlib.util
    from redox.core.common import ROOT
    spec=importlib.util.spec_from_file_location('dense_reference_test',ROOT/'scripts/validation/physical_reorg.py')
    worker=importlib.util.module_from_spec(spec);spec.loader.exec_module(worker)
    cm=pcm_object(molecule,37.5,order=11)
    K,R=worker.dense_pcm_operators_on_actual_surface(cm)
    assert K.shape[0]==len(cm.surface['grid_coords'])
    np.testing.assert_allclose(K,cm._intermediates['K'],atol=1e-12)
    np.testing.assert_allclose(R,cm._intermediates['R'],atol=1e-12)
