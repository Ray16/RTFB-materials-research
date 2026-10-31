"""Native-SMD identity checks; no inferred experimental validation."""
import numpy as np
import pytest

from redox.qm.nonequilibrium_smd import (
    attach_frozen_slow_smd,smd_object,smd_dielectrics,cds_relaxation_difference,
)


@pytest.fixture
def molecule():
    from pyscf import gto
    return gto.M(atom='O 0 0 0; H 0 0 .96; H .92 0 -.24',basis='sto-3g',verbose=0)


def test_smd_reference_recovers_native_energy_including_cds(molecule):
    from pyscf import scf
    dm=scf.RHF(molecule).get_init_guess()
    native=scf.RHF(molecule).SMD(smd_object(molecule,order=11))
    neq=attach_frozen_slow_smd(scf.RHF(molecule),dm,order=11)
    assert abs(native.energy_tot(dm)-neq.energy_tot(dm))<1e-11
    assert abs(native.with_solvent.get_cds()-neq.with_solvent.get_cds())<1e-14
    np.testing.assert_array_equal(native.with_solvent.surface['grid_coords'],
                                  neq.with_solvent.surface['grid_coords'])


def test_cds_is_charge_and_density_independent_at_same_coordinates(molecule):
    from pyscf import scf,gto
    other=gto.M(atom=molecule.atom,basis='sto-3g',charge=1,spin=1,verbose=0)
    a=smd_object(molecule,order=11);b=smd_object(other,order=11)
    assert abs(a.get_cds()-b.get_cds())<1e-14
    assert abs(cds_relaxation_difference(a,b))<1e-14
    dm=scf.RHF(molecule).get_init_guess()
    neq=attach_frozen_slow_smd(scf.RHF(molecule),dm,order=11)
    neq.energy_tot(dm*.9)
    assert abs(neq.with_solvent.get_cds()-a.get_cds())<1e-14


def test_cds_cancels_in_symmetric_directional_mean(molecule):
    from pyscf import gto
    displaced=gto.M(atom='O 0 0 0; H 0 0 1.03; H .98 0 -.24',basis='sto-3g',verbose=0)
    a=smd_object(molecule,order=11);b=smd_object(displaced,order=11)
    forward=cds_relaxation_difference(a,b)
    reverse=cds_relaxation_difference(b,a)
    assert abs(forward)>1e-8
    assert forward+reverse==0.


def test_native_mecn_descriptor_is_used_and_invalid_overrides_rejected(molecule):
    from pyscf import scf
    constants=smd_dielectrics()
    assert constants['eps_s']==35.688
    assert constants['eps_op']==pytest.approx(1.3442**2)
    dm=scf.RHF(molecule).get_init_guess()
    with pytest.raises(ValueError):
        attach_frozen_slow_smd(scf.RHF(molecule),dm,eps_s=1.,eps_op=2.)


def test_open_shell_reference_uses_the_same_total_electrostatic_density(molecule):
    from pyscf import scf,gto
    ion=gto.M(atom=molecule.atom,basis='sto-3g',charge=1,spin=1,verbose=0)
    spin_density=scf.UHF(ion).get_init_guess()
    target=scf.RHF(molecule).get_init_guess()
    separate=attach_frozen_slow_smd(scf.RHF(molecule),spin_density,order=11)
    total=attach_frozen_slow_smd(scf.RHF(molecule),spin_density.sum(axis=0),order=11)
    assert abs(separate.energy_tot(target)-total.energy_tot(target))<1e-11
    np.testing.assert_allclose(separate.with_solvent._get_vind(target)[1],
                               total.with_solvent._get_vind(target)[1],atol=1e-11,rtol=0)
