"""Synthetic (no-QM) unit tests for the solution-phase inner-sphere lambda_i machinery:
the 4-point Nelsen arithmetic, single-conformer scope, and the
structural-nonnegativity contract. All inputs are hand-made energies; no DFT is run."""
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


# --------------------------------------------------------------- unsupported ensembles
@pytest.mark.parametrize("counts", [(2, 2), (1, 2), (2, 1)])
@pytest.mark.parametrize("operation", [reorg.lambda_for_couple, reorg._couple_cross_jobs])
def test_ensemble_refused_even_if_only_one_family_matches(monkeypatch, counts, operation):
    def conformers(gid, state):
        n = counts[0 if state == "O" else 1]
        return [dict(key=f"conf_{i}") for i in range(n)], True

    monkeypatch.setattr(reorg, "_state_conformers", conformers)
    with pytest.raises(ValueError, match="state-population convention.*free energies"):
        operation("test", ("O", 0, 1), ("R", -1, 2))


def test_group_preflight_refuses_later_ensemble_before_any_calculation(monkeypatch):
    O, R, R2 = ("O", 0, 1), ("R", -1, 2), ("R2", -2, 1)
    monkeypatch.setattr(reorg, "_couples", lambda gid: [(O, R), (R, R2)])
    monkeypatch.setattr(reorg, "_state_conformers", lambda gid, state: (
        [dict(key=f"conf_{i}", geom=reorg.DFT / "unused")
         for i in range(2 if state == "R2" else 1)], True))
    monkeypatch.setattr(reorg, "_cross_energy", lambda *a, **kw: pytest.fail("started QM"))
    with pytest.raises(ValueError, match="conformer averaging is unsupported"):
        reorg.compute_group("test")


def test_single_conformer_keeps_both_distortions_and_gas_diagnostic(monkeypatch, tmp_path):
    monkeypatch.setattr(reorg, "DFT", tmp_path)
    monkeypatch.setattr(reorg, "_state_conformers", lambda gid, state: ([dict(
        key="min", geom=tmp_path / state,
        e_smd_eV=-10.0 if state == "O" else -12.0,
        e_gas_eV=-9.0 if state == "O" else -11.0)], False))
    monkeypatch.setattr(reorg, "_read_cross", lambda gid, state, *a: dict(
        e_smd_eV=-9.7 if state == "O" else -11.5,
        e_gas_eV=-8.9 if state == "O" else -10.8,
        smd_s_squared=0.0 if state == "O" else 0.75, anion_unbound=False))
    monkeypatch.setattr(reorg, "_res", lambda *a: {})
    monkeypatch.setattr(reorg, "_heavy_rmsd", lambda *a: None)
    monkeypatch.setattr(reorg, "_sym_heavy_rmsd", lambda *a, **kw: None)
    result = reorg.lambda_for_couple("test", ("O", 0, 1), ("R", -1, 2))
    assert result["lambda_i_eV"] == pytest.approx(0.8)
    assert result["relax_ox_meV"] == pytest.approx(300)
    assert result["relax_red_meV"] == pytest.approx(500)
    assert result["lambda_i_gas_eV"] == pytest.approx(0.3)
    assert result["n_conf_matched"] == 1
    assert result["conformer_averaged"] is False
    assert result["flag"] == ""


@pytest.mark.parametrize("energy", [None, float("nan"), float("inf")])
def test_invalid_own_energy_is_not_reported_as_a_result(monkeypatch, energy):
    monkeypatch.setattr(reorg, "_state_conformers", lambda *a: (
        [dict(key="min", e_smd_eV=energy)], True))
    assert reorg.lambda_for_couple("test", ("O", 0, 1), ("R", -1, 2)) is None


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
