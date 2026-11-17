"""Partition-function boundary tests, including a real internal mode below 1 cm^-1."""
import numpy as np
import pytest

from redox.qm import qrrho as Q


def test_installed_kernel_matches_upstream_thermo_unit_test():
    moments = np.array([15649.5560990126, 484991.786527778, 485024.382510520])
    frequencies = [117.29518854, 161.91929949, 164.05816261, 218.07506352,
                   402.96868229, 558.86387616, 1523.62725135, 1561.31279119,
                   3460.43559680, 3633.65239963, 3636.81476002, 3665.34271369]
    rotational = .5 * Q.AUTORCM / moments[::-1]
    average = moments.mean() / ((1 / 9.10938356e-31) * (1 / .52917726)**2 * 1e20)
    values = Q.kernel_thermal(frequencies, rotational, average,
                              2 * 15.99940492 + 4 * 1.00794075, 1, cutoff_cm=100)
    assert values['G_Ha'] == pytest.approx(.017417501635698, abs=1e-12)
    assert values['H_Ha'] == pytest.approx(.050275771916811, abs=1e-12)
    assert values['TS_Ha'] == pytest.approx(.032858270281112, abs=1e-12)


def test_soft_stretch_survives_geometric_projection_and_contributes_to_partition():
    mass = Q.atomic_masses([1, 1])
    stretch = np.array([-1., 0., 0., 1., 0., 0.]) / np.sqrt(2)
    hessian = (.391 / Q.AUTORCM)**2 * mass[0] * Q.AMUTOAU * np.outer(stretch, stretch)
    coordinates = np.array([[-.5, 0., 0.], [.5, 0., 0.]])
    result = Q.nuclear_thermals(hessian, coordinates, [1, 1], 2)
    assert result['rigid_rank'] == 5 and result['internal_mode_count'] == 1
    assert result['frequencies_cm'] == pytest.approx([.391], abs=1e-10)
    assert result['retained_soft_mode_indices'] == [0]
    rot, average, _ = Q.inertia_properties(coordinates, mass)
    empty = Q.kernel_thermal([], rot, average, sum(mass), 2, linear=True)['G_Ha'] * Q.HARTREE_EV
    assert abs(result['nuclear_by_cutoff_eV']['50'] - empty) > .001


def test_projection_does_not_hide_an_internal_instability():
    stretch = np.array([-1., 0., 0., 1., 0., 0.]) / np.sqrt(2)
    mass = Q.atomic_masses([1, 1])
    hessian = -(20 / Q.AUTORCM)**2 * mass[0] * Q.AMUTOAU * np.outer(stretch, stretch)
    coordinates = np.array([[-.5, 0., 0.], [.5, 0., 0.]])
    frequency, rank = Q.internal_vibrations(hessian, coordinates, mass)
    assert rank == 5 and frequency[0] == pytest.approx(-20)
    with pytest.raises(ValueError, match='unstable'):
        Q.nuclear_thermals(hessian, coordinates, [1, 1], 2)


@pytest.mark.parametrize('frequency', [[0.], [-.1], [np.nan]])
def test_kernel_never_clamps_or_discards_bad_modes(frequency):
    with pytest.raises(ValueError):
        Q.kernel_thermal(frequency, [1., 2., 3.], 1e-45, 28., 1)


def test_exact_zero_internal_mode_is_not_given_an_artificial_frequency():
    with pytest.raises(ValueError, match='zero internal'):
        Q.nuclear_thermals(np.zeros((6, 6)), [[-.5, 0., 0.], [.5, 0., 0.]], [1, 1], 2)


def test_different_binary_is_refused_before_loading(monkeypatch, tmp_path):
    path = tmp_path/'lib.so';path.write_bytes(b'unqualified replacement')
    Q._checked_library.cache_clear()
    monkeypatch.setattr(Q, 'LIBRARY', path)
    try:
        with pytest.raises(RuntimeError, match='qualified build'):
            Q.atomic_masses([1])
    finally:
        Q._checked_library.cache_clear()


def test_kernel_process_is_separate_and_uses_capped_native_pools():
    import os
    Q.atomic_masses([1])
    worker=Q.kernel_execution()
    assert worker['pid']!=os.getpid()
    assert set(worker['thread_limits'].values())=={'1'}


def test_kernel_still_works_with_pyscf_D4_loaded_in_the_parent():
    from pyscf.dispersion import dftd4
    assert dftd4 is not None
    Q._close_worker()
    assert Q.atomic_masses([8])==pytest.approx([15.99940492])
    result=Q.kernel_thermal([100.],[1.,2.,3.],1e-45,28.,1)
    assert np.isfinite(result['G_Ha'])
