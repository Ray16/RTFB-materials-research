"""Data-contract tests: thermodynamic identities, state ordering, cache/protocol invalidation,
contiguous capacity, Pareto completeness/sigma, outer-sphere limits.

  PYTHONPATH=src python -m pytest -q tests
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import pytest

import redox.core.common as common
from redox.core import protocol
from redox.screening import scorecard
from redox.properties import stability
from redox.properties.lambda_outer import born_lambda_outer, sphere_sanity_check, _pekar_lambda_ha
from redox.screening.pareto import _dominates, _pareto_front


# ----------------------------------------------------------------------- helpers / fixtures
def _write_state(root: Path, gid, st, q, m, e_smd=None, g_th=1.0, xyz="H 0 0 0\nH 0 0 0.74",
                 n_imag=0):
    sd = root / gid / st
    sd.mkdir(parents=True, exist_ok=True)
    lines = xyz.strip().splitlines()
    (sd / "opt.xyz").write_text(f"{len(lines)}\n\n" + "\n".join(lines) + "\n")
    res = dict(charge=q, mult=m, optimized=True, g_thermal_eV=g_th, n_imag=n_imag,
               e_smd_eV=e_smd, e_gas_eV=e_smd, sp_xc="wb97m-v", sp_basis="def2-tzvp")
    (sd / "result.json").write_text(json.dumps(res))
    return sd


def _write_record(sd: Path, q, m, e_smd, e_gas=None, converged=True):
    rp = protocol.record_path(sd, q, m)
    rp.parent.mkdir(parents=True, exist_ok=True)
    rp.write_text(json.dumps(dict(e_smd_eV=e_smd, e_gas_eV=e_gas if e_gas is not None else e_smd,
                                  converged_smd=converged, converged_gas=converged,
                                  protocol_hash=rp.stem)))
    return rp


@pytest.fixture
def tree(tmp_path, monkeypatch):
    monkeypatch.setattr(common, "DFT", tmp_path)
    return tmp_path


# ------------------------------------------------------------------ protocol / cache contract
def test_protocol_hash_changes_with_every_determinant():
    base = protocol.protocol_hash("g1", 0, 1)
    assert base == protocol.protocol_hash("g1", 0, 1)
    assert base != protocol.protocol_hash("g2", 0, 1)              # geometry
    assert base != protocol.protocol_hash("g1", -1, 2)             # charge / spin
    assert base != protocol.protocol_hash("g1", 0, 1, dict(protocol.ACTIVE_SP, basis="def2-tzvp"))
    assert base != protocol.protocol_hash("g1", 0, 1, dict(protocol.ACTIVE_SP, conv_tol=1e-6))


def test_active_protocol_is_uniform_diffuse():
    from redox.qm import dft
    assert dft.SP_BASIS == dft.SP_BASIS_ANION == protocol.ACTIVE_SP["basis"] == "def2-tzvpd"


def test_read_result_never_mixes_levels(tree):
    sd = _write_state(tree, "m", "neu", 0, 1, e_smd=-100.0)       # legacy energy in result.json
    r = common.read_result("m", "neu", root=tree)
    assert r["sp_status"] == "missing" and r["e_smd_eV"] is None   # legacy value NOT used
    assert r["raw_e_smd_eV"] == -100.0
    _write_record(sd, 0, 1, -101.0)
    r = common.read_result("m", "neu", root=tree)
    assert r["sp_status"] == "ok" and r["e_smd_eV"] == -101.0


def test_record_invalidated_by_geometry_change_and_unconverged(tree):
    sd = _write_state(tree, "m", "neu", 0, 1)
    _write_record(sd, 0, 1, -101.0)
    (sd / "opt.xyz").write_text("2\n\nH 0 0 0\nH 0 0 0.80\n")         # geometry moved
    assert common.read_result("m", "neu", root=tree)["e_smd_eV"] is None
    _write_record(sd, 0, 1, -102.0, converged=False)
    assert common.read_result("m", "neu", root=tree)["e_smd_eV"] is None


def test_state_names_skip_helper_dirs(tree):
    for st in ("neu", "red1"):
        _write_state(tree, "m", st, 0, 1)
    for helper in ("_spincheck_red1_m3", "reorg", "sp", "neu.bak_pretscan"):
        _write_state(tree, "m", helper, -1, 3)
    assert common.state_names("m", root=tree) == ["neu", "red1"]


def test_reorg_cache_rejects_other_level():
    from redox.properties.reorg import _cache_ok, CACHE_SCHEMA
    d = dict(e_gas_eV=-1.0, cache_schema=CACHE_SCHEMA, gas_homo_eV=-5, anion_unbound=False,
             charge=0, mult=1, geom_sha1="abc", sp_xc="wb97m-v", sp_basis="def2-tzvpd",
             sp_nlc="vv10", sp_disp=None, converged_gas=True)
    assert _cache_ok(d, 0, 1, "abc")
    assert not _cache_ok(dict(d, sp_basis="def2-tzvp"), 0, 1, "abc")      # basis changed
    assert not _cache_ok(dict(d, sp_xc="b3lyp"), 0, 1, "abc")             # functional changed
    assert not _cache_ok(d, 0, 1, "xyz")                                   # geometry changed
    assert not _cache_ok(dict(d, converged_gas=False), 0, 1, "abc")       # unconverged


# ------------------------------------------------------------- thermal gate / identities
def test_free_energy_never_substitutes_zero_thermal():
    assert common.free_energy(dict(e_smd_eV=-1.0, g_thermal_eV=None)) == (None, "missing")
    assert common.free_energy(dict(e_smd_eV=-1.0, g_thermal_eV=0.5, n_imag=0)) == (-0.5, "ok")
    assert common.free_energy(dict(e_smd_eV=-1.0, g_thermal_eV=0.5, n_imag=2))[1] == "imag_modes"


def test_disproportionation_equals_wave_spacing(monkeypatch):
    # G ladder: neu (0) -> red1 (-1) -> red2 (-2); E_k = -(G_R - G_O)
    G = {"neu": -10.0, "red1": -13.2, "red2": -15.9}
    monkeypatch.setattr(stability, "G", lambda gid, s: G[s])
    monkeypatch.setattr(stability, "_states_by_charge",
                        lambda gid: [("neu", 0), ("red1", -1), ("red2", -2)])
    (row,) = stability.disproportionation("x")
    E1 = -(G["red1"] - G["neu"]); E2 = -(G["red2"] - G["red1"])
    assert row["intermediate"] == "red1"
    assert row["dG_disp_eV"] == pytest.approx(E1 - E2, abs=1e-4)          # dG_disp = F(E1-E2)


def test_disproportionation_incomplete_when_thermal_missing(monkeypatch):
    monkeypatch.setattr(stability, "G", lambda gid, s: None if s == "red2" else -1.0)
    monkeypatch.setattr(stability, "_states_by_charge",
                        lambda gid: [("neu", 0), ("red1", -1), ("red2", -2)])
    assert stability.disproportionation("x") == []


def test_integrity_couples_ordered_by_charge(tree, monkeypatch):
    from redox.properties import integrity
    monkeypatch.setattr(integrity, "state_names", lambda gid: ["red2", "neu", "red1"])
    q = {"neu": 0, "red1": -1, "red2": -2}
    monkeypatch.setattr(integrity, "_res", lambda gid, s, raw=False: {"charge": q[s]})
    assert integrity._couples("x") == [(("neu", 0), ("red1", -1)), (("red1", -1), ("red2", -2))]


# ------------------------------------------------------------------ contiguous capacity
def _c(couple, E, integ="intact_bound", status="ok", lo=-3.0, hi=2.0):
    return dict(couple=couple, E=E, status=status, integrity=integ, in_window=(E is not None and lo <= E <= hi))


def test_second_electron_needs_first():
    couples = {0: _c("neu->red1", -3.5), -1: _c("red1->red2", -1.0)}     # 1st outside window
    paths, stop = scorecard.contiguous_paths(couples, 0, 0.0)
    assert paths["anolyte"] == [] and "outside window" in stop["anolyte"]


def test_ambipolar_split_not_double_counted():
    couples = {1: _c("ox->rad", 0.3), 0: _c("rad->red", -1.9)}           # TEMPO-like, rest q=0
    paths, _ = scorecard.contiguous_paths(couples, 0, 0.0)
    assert [c["couple"] for c in paths["anolyte"]] == ["rad->red"]
    assert [c["couple"] for c in paths["catholyte"]] == ["ox->rad"]


def test_missing_data_is_incomplete_not_rejected():
    couples = {0: _c("neu->red1", None, status="INCOMPLETE:neu(no_energy)")}
    paths, stop = scorecard.contiguous_paths(couples, 0, 0.0)
    assert paths["anolyte"] == [] and stop["anolyte"].startswith("INCOMPLETE")
    couples = {0: _c("neu->red1", -1.0, integ=None)}
    _, stop = scorecard.contiguous_paths(couples, 0, 0.0)
    assert stop["anolyte"].startswith("INCOMPLETE")


def test_bond_change_stops_path():
    couples = {0: _c("neu->red1", -1.0), -1: _c("red1->red2", -1.5, integ="bond_change")}
    paths, stop = scorecard.contiguous_paths(couples, 0, 0.0)
    assert len(paths["anolyte"]) == 1 and "bond_change" in stop["anolyte"]


# ------------------------------------------------------------------------------ Pareto
class _SC:
    SIGMA_CAPACITY = 0.0; SIGMA_DISP_EV = 0.17; SIGMA_LAMBDA_EV = 0.10


def _row(i, E, cap, disp, lam, sl=0.2, n=2):
    return dict(id=i, E_V=E, sigma_E_V=0.1, specific_capacity_mAh_g=cap, dG_disp_kJmol=disp,
                sigma_disp_eV=0.17, disp_applicable=str(n >= 2), lambda_i_eV=lam,
                sigma_lambda_eV=(sl if lam is not None else None))


def test_incomplete_candidate_cannot_dominate():
    good = _row("complete", -1.0, 100, 50, 0.40)
    better_but_no_lambda = _row("nolam", -1.6, 160, 80, None)
    front, _, dom, missing = _pareto_front([good, better_but_no_lambda], "anolyte", _SC)
    assert front == ["complete"]
    assert dom["complete"] == []                           # NOT dominated by the incomplete row
    assert missing["nolam"] == ["kinetics"]


def test_pareto_uses_row_sigma_not_config():
    a = _row("a", -1.0, 100, 50, 0.20, sl=0.21)
    b = _row("b", -1.0, 100, 50, 0.50, sl=0.21)
    front, objs, dom, _ = _pareto_front([a, b], "anolyte", _SC)
    # dominance is decided with the ROW sigma (0.21), not the 0.10 config fallback
    assert _dominates(objs["a"], objs["b"]) == (0.30 > math.hypot(0.21, 0.21))
    assert objs["a"]["kinetics"][1] == 0.21                # the row's own sigma


def test_single_electron_stability_not_applicable():
    r = _row("one_e", -1.0, 100, None, 0.4, n=1)
    _, _, _, missing = _pareto_front([r], "anolyte", _SC)
    assert missing["one_e"] == []


# --------------------------------------------------------------------- outer sphere
def test_pcm_point_charge_sphere_equals_born():
    for el in ("He", "Ne"):
        d = sphere_sanity_check(element=el)
        assert abs(d["rel_err"]) < 1e-4
        assert d["pcm_eV"] > 0


def test_pekar_lambda_sign_error_is_raised():
    import numpy as np
    K = np.eye(3); R = -0.5 * np.eye(3); dV = np.ones(3)
    # physical order: the static (eps_s) response is the stronger one -> lambda > 0
    assert _pekar_lambda_ha(K, 2 * R, K, R, dV) > 0
    # swapped dielectrics (op <-> s) flip the sign: must be REJECTED, not abs()'d
    with pytest.raises(ValueError):
        _pekar_lambda_ha(K, R, K, 2 * R, dV)


def test_self_exchange_limits():
    lo1, a = 0.8, 4.0
    assert scorecard.lambda_o_self_exchange(lo1, a, 2 * a) == pytest.approx(lo1)
    assert scorecard.lambda_o_self_exchange(lo1, a, 1e9) == pytest.approx(2 * lo1)
    with pytest.raises(ValueError):
        scorecard.lambda_o_self_exchange(lo1, a, 1.5 * a)


def test_born_scaling():
    assert born_lambda_outer(2.0, 37.5, 1.806) == pytest.approx(2 * born_lambda_outer(4.0, 37.5, 1.806))


def test_sigma_uses_only_grounded_points_and_never_the_reference():
    """sigma_E must come from grounded experimental points only: the Fc reference (0 = 0 by
    construction) and anchors with grounded=False must not enter the RMSE."""
    rows = [dict(id="ferrocene", event="ox->neu", E_vs_Fc_V="0.0"),
            dict(id="tempo_parent", event="ox->rad", E_vs_Fc_V="0.449"),       # grounded, +0.200
            dict(id="phenothiazine_parent", event="ox->neu", E_vs_Fc_V="5.0")]  # ungrounded
    pts = scorecard.experimental_e_points(rows)
    assert [(p["id"], p["event"]) for p in pts] == [("tempo_parent", "ox->rad")]   # once
    assert pts[0]["res"] == pytest.approx(0.200)   # 0.449 - 0.249 (Gerken & Stahl 2015)
    _, pooled = scorecard.derive_family_sigma(rows)
    assert pooled == pytest.approx(0.200)


def test_reorg_cache_requires_explicit_convergence():
    """A cross point without converged_gas=True is unverified and must not be used."""
    from redox.properties import reorg
    from redox.core.protocol import ACTIVE_SP
    d = dict(e_gas_eV=-1.0, sp_xc=ACTIVE_SP["xc"], sp_basis=ACTIVE_SP["basis"],
             sp_nlc=ACTIVE_SP["nlc"], sp_disp=ACTIVE_SP["disp"], cache_schema=reorg.CACHE_SCHEMA,
             charge=0, mult=1, geom_sha1="g", **{k: None for k in reorg.REQUIRED_DIAG})
    assert not reorg._cache_ok(dict(d), 0, 1, "g")            # flag missing
    assert reorg._cache_ok(dict(d, converged_gas=True), 0, 1, "g")


def test_adopted_record_needs_density_fit_provenance():
    """An energy adopted from result.json must prove density fitting + SCF tolerance; anion
    states optimized before RI-J became the default were adopted without it (7-12 meV off)."""
    P = protocol.ACTIVE_SP
    assert protocol.record_provenance_ok(dict(source="computed:gpu"))
    assert not protocol.record_provenance_ok(dict(source="adopted:result.json"))
    assert protocol.record_provenance_ok(dict(source="adopted:result.json", density_fit=True,
                                              conv_tol=P["conv_tol"]))
    from redox.qm.sp import _adoptable
    res = dict(optimized=True, sp_xc=P["xc"], sp_basis=P["basis"], sp_nlc=P["nlc"],
               sp_disp=P["disp"], solvent=P["solvent"], converged_smd=True, converged_gas=True,
               e_smd_eV=-1.0, e_gas_eV=-1.0, mult=2)
    assert not _adoptable(res)                                         # no provenance
    assert not _adoptable(dict(res, density_fit=True, conv_tol=P["conv_tol"],
                               scf_guesses=["minao"]))                  # open shell, no sweep
    assert _adoptable(dict(res, density_fit=True, conv_tol=P["conv_tol"],
                           scf_guesses=["minao", "atom", "huckel"]))


def test_lowest_scf_keeps_lowest_converged_guess():
    from redox.qm import dft

    class FakeMF:
        E = {"minao": -1.0, "atom": -1.5, "huckel": -9.0}
        CONV = {"minao": True, "atom": True, "huckel": False}

        def __init__(self):
            self.init_guess = None

    calls = []

    def fake_kernel(mf):
        calls.append(mf.init_guess)
        return FakeMF.E[mf.init_guess], FakeMF.CONV[mf.init_guess]

    orig = dft._kernel_robust
    dft._kernel_robust = fake_kernel
    try:
        e, conv, _, per = dft._lowest_scf(FakeMF, ("minao", "atom", "huckel"))
    finally:
        dft._kernel_robust = orig
    assert calls == ["minao", "atom", "huckel"]
    assert (e, conv) == (-1.5, True)          # lowest CONVERGED (huckel did not converge)
    assert len(per) == 3
    dft._kernel_robust = fake_kernel
    try:
        _, _, post, _ = dft._lowest_scf(FakeMF, ("minao", "atom", "huckel"),
                                        post=lambda m: m.init_guess)
    finally:
        dft._kernel_robust = orig
    assert post == "atom"                     # post() result belongs to the kept guess


def test_discovery_rules_rediscover_every_current_candidate():
    """The systematic identification rules (redox.screening.discovery), applied to the
    precursor of each grafted candidate in the configs, must regenerate exactly that
    candidate and classify its electron count like the config does."""
    from redox.screening.discovery import rediscover
    groups = [g for cfg in ("starting_candidates", "merrifield_multielectron", "redox_groups")
              for g in common.load_config(cfg).GROUPS if "frag" in g]
    rows = rediscover(groups, faraday=96485.33212)
    assert rows and all(r["rediscovered"] for r in rows), [r["id"] for r in rows if not r["rediscovered"]]
    multi = {r["id"]: r["multi_e"] for r in rows}
    for gid in ("viologen", "ethylviologen", "bisviologen", "pmdi", "ndi_ammonium",
                "mophquinone", "dmophquinone", "aq_benzyloxy", "nq_benzyloxy", "aq_benzylamino"):
        assert multi[gid], gid
    for gid in ("pyridinium", "phenothiazine", "tempo", "dtbc_phenol"):
        assert not multi[gid], gid
