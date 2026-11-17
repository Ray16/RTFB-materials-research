"""Synthetic (no-QM) unit tests for the solution-phase inner-sphere lambda_i machinery:
the 4-point Nelsen arithmetic, matched-conformer pairing + Boltzmann weighting, and the
structural-nonnegativity contract. All inputs are hand-made energies; no DFT is run."""
import math

import pytest

from redox.properties import reorg


# --------------------------------------------------------------- 4-point Nelsen arithmetic
def test_four_point_lambda_is_sum_of_two_distortions():
    # Own minima at -10.0 (O) and -12.0 (R); each distorts UPWARD onto the partner geometry.
    E_O_at_O, E_O_at_R = -10.0, -9.7      # O costs +0.3 to sit on R's geometry
    E_R_at_R, E_R_at_O = -12.0, -11.5     # R costs +0.5 to sit on O's geometry
    lam, relax_ox, relax_red = reorg._four_point_lambda(E_O_at_O, E_O_at_R, E_R_at_R, E_R_at_O)
    assert relax_ox == pytest.approx(0.3)
    assert relax_red == pytest.approx(0.5)
    assert lam == pytest.approx(0.8)      # lambda_i = sum, NOT mean


def test_four_point_lambda_nonnegative_when_own_points_are_minima():
    # On a self-consistent surface the own point is the minimum => both halves >= 0 => lam >= 0.
    for d_ox, d_red in [(0.0, 0.0), (0.05, 0.2), (0.4, 0.4), (1.0, 0.01)]:
        lam, ro, rr = reorg._four_point_lambda(-5.0, -5.0 + d_ox, -6.0, -6.0 + d_red)
        assert ro >= 0 and rr >= 0 and lam >= 0


def test_four_point_lambda_negative_half_signals_non_stationary_point():
    # If a "cross" energy is BELOW the claimed own minimum, a half goes negative -> SCF/geometry
    # fault (exactly the artifact the all-SMD surface removes). The arithmetic exposes it.
    lam, ro, rr = reorg._four_point_lambda(-5.0, -5.3, -6.0, -6.0)
    assert ro < 0 and lam < 0


# --------------------------------------------------------------- conformer matching
def test_match_conformers_pairs_shared_families():
    O = [dict(key="conf_00"), dict(key="conf_01")]
    R = [dict(key="conf_00"), dict(key="conf_01")]
    pairs, repaired = reorg.match_conformers(O, R)
    assert not repaired
    assert sorted((a, b) for a, b, _ in pairs) == [(0, 0), (1, 1)]


def test_match_conformers_repairs_disjoint_families():
    # O and R landed in DIFFERENT families (no shared key): re-pair lowest-to-lowest instead
    # of discarding the couple (the old conformer_jump behaviour).
    O = [dict(key="A0"), dict(key="A1")]
    R = [dict(key="B0"), dict(key="B1")]
    pairs, repaired = reorg.match_conformers(O, R)
    assert repaired
    assert [(a, b) for a, b, _ in pairs] == [(0, 0), (1, 1)]


def test_match_conformers_single_minimum_fallback():
    pairs, repaired = reorg.match_conformers([dict(key="min")], [dict(key="min")])
    assert pairs == [(0, 0, "min")] and not repaired


# --------------------------------------------------------------- Boltzmann weighting
def test_boltzmann_weights_single_family_is_unity():
    assert reorg._boltzmann_weights([-123.4]) == [1.0]


def test_boltzmann_weights_favor_lower_energy_and_sum_to_one():
    w = reorg._boltzmann_weights([0.0, reorg.KT_EV * math.log(2)])  # second is 2x less likely
    assert w[0] == pytest.approx(2 * w[1])
    assert sum(w) == pytest.approx(1.0)


def test_boltzmann_average_of_lambda_lies_between_family_values():
    fams = [(0.0, 0.5), (0.10, 0.9)]   # (relative energy eV, lambda eV)
    w = reorg._boltzmann_weights([f[0] for f in fams])
    lam = sum(wi * f[1] for wi, f in zip(w, fams))
    assert 0.5 <= lam <= 0.9
    assert lam == pytest.approx(0.5 * w[0] + 0.9 * w[1])


# --------------------------------------------------------------- cache contract (v4 = SMD)
def test_cache_schema_is_smd_and_requires_smd_energy():
    assert reorg.CACHE_SCHEMA == 4
    from redox.core.protocol import ACTIVE_SP
    base = dict(cache_schema=4, sp_xc=ACTIVE_SP["xc"], sp_basis=ACTIVE_SP["basis"],
                sp_nlc=ACTIVE_SP["nlc"], sp_disp=ACTIVE_SP["disp"], charge=0, mult=1,
                geom_sha1="g", protocol=ACTIVE_SP, converged_smd=True,
                **{k: 0.0 for k in reorg.REQUIRED_DIAG})
    assert reorg._cache_ok(dict(base, e_smd_eV=-1.0), 0, 1, "g")
    assert not reorg._cache_ok(dict(base, e_smd_eV=None, e_gas_eV=-1.0), 0, 1, "g")  # gas-only
    assert not reorg._cache_ok(dict(base, e_smd_eV=float("nan")), 0, 1, "g")
