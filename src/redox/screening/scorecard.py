"""Unified scorecard — the input contract for the selection engine.

Joins the per-axis tables (redox potentials, state integrity, reorganization energy,
disproportionation, capacity/proxies, outer-sphere lambda) and emits ONE ROW PER (molecule,
electrode pool), with a sigma per axis for sigma-aware domination.

STATUS: candidate | REJECTED (a real negative: no accessible electron) | INCOMPLETE (missing
data on the path — never scored as if the data were absent or zero).

ACCESSIBLE ELECTRONS = a CONTIGUOUS path from the declared resting state (manifest n_e=0):
  anolyte   : successive reductions  rest -> rest-1e -> ...  while each couple is
              intact_bound (redox.properties.integrity), inside WINDOW_V_VS_FC and below the divider;
  catholyte : successive oxidations  rest -> rest+1e -> ...  above the divider.
A second electron counts only if the first is accessible. An ambipolar molecule gets two pool
rows, each with ONLY its own side's electrons (no double counting across electrodes).
`intact_bound` is a necessary-condition screen, not proof of electrochemical reversibility.

CAPACITY = n F / mass on that path. The counter-ion-inclusive variants are SCENARIO
conventions (which ions count in the material mass depends on the device and state of charge),
not exact physics.

REORGANIZATION ENERGY — explicit conventions (Nelsen 4-point, gas, QC-clean couples only):
  lambda_i_ox_eV   = lambda_O : E_O(geom_R) - E_O(geom_O)   (distortion on the O surface)
  lambda_i_red_eV  = lambda_R : E_R(geom_O) - E_R(geom_R)
  lambda_i_eV      = lambda_O + lambda_R  (= inner-sphere lambda of a SELF-EXCHANGE pair)
  lambda_het_eV    = heterogeneous (electrode) ET: (lambda_O + lambda_R)/2  +  lambda_o,1-body
                     (one molecule reorganizes; outer sphere = 1-body molecular-cavity PCM,
                     image term neglected)
  lambda_se_contact_eV = self-exchange / polymer hopping at contact: lambda_O + lambda_R +
                     lambda_o,SE(d = 2a); lambda_o_self_exchange() gives any separation d.
  The Pareto kinetics objective is lambda_i_eV (monotone with lambda_i/2); lambda_het /
  lambda_se are reported, not ranked, because lambda_o is a continuum proxy (sigma ~0.3 eV).

  PYTHONPATH=src python -m redox.screening.scorecard
"""
from __future__ import annotations
import csv
import json
import statistics

from redox.core.common import RESULTS, UMA
from redox.core.common import load_config as _cfg
from redox.core.common import to_float as _f

# not real candidates: the Fc/Fc+ internal reference
REFERENCE_IDS = {"ferrocene"}


# Screening pools, in the order they were added. Each entry is a config module whose GROUPS
# are grafted screening candidates; `batch` is carried onto every scorecard row so downstream
# plots/tables can distinguish the cohorts.
CANDIDATE_CONFIGS = [
    ("starting", "starting_candidates"),
    ("merrifield_multi", "merrifield_multielectron"),
]


def _candidate_batches():
    """{id: batch} over every grafted screening candidate.

    `viologen` (the sheet's methylviologen) is deduped into config/redox_groups.py rather than
    config/starting_candidates.py, so it is added to the 'starting' batch explicitly.

    The scorecard is scoped to these ids — the old exploratory redox_groups, the validation
    cores, and the standalone `*_sa` reference forms all have data in the per-axis tables but
    are NOT screening candidates and are excluded here."""
    out = {}
    for batch, mod in CANDIDATE_CONFIGS:
        for g in _cfg(mod).GROUPS:
            out[g["id"]] = batch
    out["viologen"] = "starting"
    return out


def _candidate_chem_flags():
    """{id: note} for candidates whose modelled chemistry departs from the source sheet
    (see the `flag` key in config/merrifield_multielectron.py). Empty string when clean."""
    out = {}
    for _, mod in CANDIDATE_CONFIGS:
        for g in _cfg(mod).GROUPS:
            if g.get("flag"):
                out[g["id"]] = g["flag"]
    return out


def _unrankable_ids():
    """Candidates that are computed and REPORTED but must not enter the Pareto ranking,
    because their electrochemistry is not a reversible outer-sphere couple (config key
    `rankable=False`). Ranking an EC / bond-making mechanism against clean outer-sphere
    couples would not be like-for-like — the same reason the metal-oxo rows were excluded."""
    out = set()
    for _, mod in CANDIDATE_CONFIGS:
        for g in _cfg(mod).GROUPS:
            if g.get("rankable") is False:
                out.add(g["id"])
    return out


def _candidate_ids():
    return set(_candidate_batches())


SPINCHECK = UMA.parent / "spincheck"     # calcs/spincheck/<id>/<state>_m<mult>/result.json


def _min_spin_gap_eV(gid):
    """(smallest |spin-state gap| over the molecule's states, source).

    Per state the gap comes from a DFT check when one exists — calcs/spincheck/<id>/
    <state>_m<k>/result.json, the alternative multiplicity optimized at the production DFT
    level; gap = G(alt) - G(chosen) — otherwise from the UMA multiplicity scan. UMA is only
    the pre-optimizer and can badly underestimate gaps (aq_benzyloxy dianion S-T: UMA 0.19 eV,
    DFT 0.81 eV). None when no state had an alternative multiplicity. A small value means
    the spin ground state is near-degenerate -> lower-confidence potentials."""
    from redox.core.common import free_energy, read_result
    gaps, srcs = {}, {}
    gdir = UMA / gid
    if gdir.is_dir():
        for rj in gdir.glob("*/result.json"):
            try:
                g = json.loads(rj.read_text()).get("spin_gap_eV")
            except (OSError, ValueError):
                continue
            if g is not None:
                gaps[rj.parent.name] = abs(float(g)); srcs[rj.parent.name] = "uma"
    for rj in (SPINCHECK / gid).glob("*_m*/result.json"):
        st = rj.parent.name.rsplit("_m", 1)[0]
        alt = json.loads(rj.read_text())
        Galt = (alt.get("e_smd_eV") + alt["g_thermal_eV"]
                if alt.get("e_smd_eV") is not None and alt.get("g_thermal_eV") is not None else None)
        chosen = read_result(gid, st, raw=True)
        Gch = (chosen["e_smd_eV"] + chosen["g_thermal_eV"]
               if chosen and chosen.get("g_thermal_eV") is not None else None)
        if Galt is not None and Gch is not None:
            gaps[st] = Galt - Gch            # signed: < 0 would mean the WRONG ground state
            srcs[st] = "dft"
    if not gaps:
        return None, ""
    st = min(gaps, key=lambda k: gaps[k])
    return round(gaps[st], 4), f"{srcs[st]}:{st}"


def _resting_states():
    """{id: (state, charge)} — the as-prepared resting state (manifest n_e == 0) from which
    the contiguous redox paths are traversed."""
    from redox.core.common import read_manifest
    return {r["id"]: (r["state"], int(r["charge"])) for r in read_manifest()
            if str(r.get("n_e")) == "0"}


def contiguous_paths(couples_by_qox, q0, divider):
    """Accessible electrons per electrode pool by traversing the redox graph from the resting
    charge q0. `couples_by_qox` maps q_ox -> couple dict (keys: couple, E, status, integrity,
    in_window). Anolyte walks reductions q0 -> q0-1 -> ...; catholyte walks oxidations
    q0 -> q0+1 -> ... . A step is taken only if the couple is complete, intact_bound, inside
    the window and on the pool's side of the divider; the walk STOPS at the first step that
    is not, so an electron after an inaccessible one never counts. Returns
    (paths{pool: [couples]}, stop{pool: reason}); an INCOMPLETE reason means missing data."""
    paths = {"anolyte": [], "catholyte": []}
    stop = {}
    for pool, step in (("anolyte", -1), ("catholyte", +1)):
        q = q0
        while True:
            c = couples_by_qox.get(q if step < 0 else q + 1)
            if c is None:
                stop[pool] = "end of computed ladder"; break
            if c["status"] != "ok" or c["E"] is None:
                stop[pool] = f"INCOMPLETE {c['couple']}: {c['status']}"; break
            if c["integrity"] is None or c["integrity"] == "incomplete":
                stop[pool] = f"INCOMPLETE {c['couple']}: integrity not assessed"; break
            side_ok = (c["E"] < divider) if pool == "anolyte" else (c["E"] >= divider)
            if not (c["integrity"] == "intact_bound" and c["in_window"] and side_ok):
                why = ("integrity=" + c["integrity"] if c["integrity"] != "intact_bound"
                       else ("outside window" if not c["in_window"] else "other side of divider"))
                stop[pool] = f"{c['couple']} inaccessible ({why})"; break
            paths[pool].append(c)
            q += step
    return paths, stop


def lambda_o_self_exchange(lam_o1_eV, a_A, d_A):
    """Marcus two-sphere outer-sphere lambda for self-exchange between two equal sites of
    Born-equivalent radius a at centre separation d, from the 1-body value lambda_o1:
        lambda_o,SE(d) = e^2/(4 pi eps0) Pekar (1/a - 1/d) = 2 lambda_o1 (1 - a/d)
    (-> lambda_o1 at contact d = 2a; -> 2 lambda_o1 as d -> inf). a_eff is taken from the
    molecular-cavity lambda_o1 itself (a = e^2 Pekar / (8 pi eps0 lambda_o1))."""
    if d_A < 2.0 * a_A - 1e-9:
        raise ValueError("d < 2a: cavities overlap")
    return 2.0 * lam_o1_eV * (1.0 - a_A / d_A)


def _states_of(couple):
    """State-name pair of a couple, order-free: 'ox2->ox1' and 'ox1/ox2' -> {'ox1','ox2'}."""
    for sep in ("->", "/"):
        if sep in couple:
            return frozenset(s.strip() for s in couple.split(sep))
    return frozenset([couple])


def _ion_mass(smiles):
    if not smiles:
        return None
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    return Descriptors.MolWt(Chem.MolFromSmiles(smiles))


def _max_counterion_mass(acc, m_anion, m_cation):
    """Heaviest counter-ion load (g/mol) over all states of the accessible couples:
    q>0 needs q anions, q<0 needs |q| cations. None if a needed ion mass is unknown."""
    loads = []
    for q in {q for c in acc for q in (c["q_ox"], c["q_red"])}:
        if q > 0:
            if m_anion is None:
                return None
            loads.append(q * m_anion)
        elif q < 0:
            if m_cation is None:
                return None
            loads.append(-q * m_cation)
        else:
            loads.append(0.0)
    return max(loads) if loads else 0.0


def _load(name):
    p = RESULTS / f"{name}.csv"
    if not p.exists():
        return []
    with p.open() as f:
        return list(csv.DictReader(f))


# Benchmark family label (data/raw/validation/two_wave_mecn) -> candidate family label.
_BENCH_FAMILY = (("quinone", "quinone (n-type)"), ("aromatic imide", "imide (n-type)"),
                 ("viologen", "pyridine-multi-e"), ("phenothiazine", "amine (p-type)"),
                 ("nitroxide", "nitroxide"))


def _bench_family(label: str) -> str | None:
    for prefix, fam in _BENCH_FAMILY:
        if label.startswith(prefix):
            return fam
    return None


def experimental_e_points(redox_rows):
    """Every GROUNDED experimental E deg we can score, ONE point per (molecule, couple):
      * config/benchmark_mecn.py events of kind "E" in tier A (E deg' / E1/2 vs Fc, direct or
        authors' own calibration); several tier-A sources for one couple are averaged;
      * config/validation.py events with grounded=True for couples the benchmark does not
        cover, EXCEPT the ferrocene reference (E = 0 by construction on both sides — a
        residual of 0 that would flatter the RMSE).
    Returns [dict(id, event, family, exp, calc, res, n_sources, source)]; couples without a
    computed E_vs_Fc (INCOMPLETE) are skipped."""
    comp = {(r["id"], r["event"]): _f(r.get("E_vs_Fc_V")) for r in redox_rows}
    exp = {}
    for b in getattr(_cfg("benchmark_mecn"), "BENCHMARK", []):
        if b.get("excluded"):
            continue
        for ev in b["events"]:
            if ev["kind"] == "E" and ev["tier"] == "A":
                k = (b["id"], ev["event"])
                d = exp.setdefault(k, dict(family=_bench_family(b["bench_family"]), vals=[],
                                           src=[]))
                d["vals"].append(ev["exp_V_vs_Fc"]); d["src"].append(ev["doi"])
    for e in getattr(_cfg("validation"), "VALIDATION", []):
        if e["id"] == "ferrocene":
            continue
        for ev in e.get("events", []):
            k = (e["id"], ev["event"])
            if ev.get("grounded") and k not in exp:
                exp[k] = dict(family=e.get("family"), vals=[ev["exp_V_vs_Fc"]],
                              src=["config/validation.py"])
    pts = []
    for (gid, event), d in exp.items():
        c = comp.get((gid, event))
        if c is None:
            continue
        x = sum(d["vals"]) / len(d["vals"])
        pts.append(dict(id=gid, event=event, family=d["family"], exp=x, calc=c, res=c - x,
                        n_sources=len(d["vals"]), source="; ".join(d["src"])))
    return pts


def derive_family_sigma(redox_rows):
    """Per-family redox sigma (V) + pooled sigma = RMSE of computed - experimental over the
    GROUNDED experimental points (experimental_e_points). RMSE keeps any systematic bias in the
    sigma. A family with <2 points is omitted (callers fall back to the pooled RMSE).
    Returns (family_sigma, pooled_sigma)."""
    import math
    pts = experimental_e_points(redox_rows)
    by_fam = {}
    for p in pts:
        if p["family"]:
            by_fam.setdefault(p["family"], []).append(p["res"])

    def _rmse(xs):
        return math.sqrt(sum(x * x for x in xs) / len(xs))

    pooled = round(_rmse([p["res"] for p in pts]), 3) if pts else None
    family_sigma = {f: round(_rmse(xs), 3) for f, xs in by_fam.items() if len(xs) >= 2}
    return family_sigma, pooled


def _rmse(xs):
    import math
    xs = [x for x in xs if x is not None]
    return math.sqrt(sum(x * x for x in xs) / len(xs)) if xs else None


def _mad_sigma(xs):
    """Robust (outlier-resistant) scale: 1.4826 x median-absolute-deviation, which equals the
    std for Gaussian data but ignores a heavy tail. Used where a known outlier population
    (e.g. D3TaLES unbound anions) would inflate a plain RMSE."""
    xs = sorted(x for x in xs if x is not None)
    if not xs:
        return None
    def _median(v):
        n = len(v)
        return v[n // 2] if n % 2 else 0.5 * (v[n // 2 - 1] + v[n // 2])
    med = _median(xs)
    mad = _median(sorted(abs(x - med) for x in xs))
    return 1.4826 * mad


def derive_uncertainties(redox_rows):
    """ALL model uncertainties, COMPUTED from benchmark result files at build time — nothing
    hard-coded, each tightens as its benchmark grows. Each falls back to the documented prior
    in config/scorecard_config only if its benchmark file is missing. Returns a dict:
      family_redox : per-family redox sigma (V) from validation residuals (>=2 anchors)
      pooled_redox : pooled validation RMSE (V) — used for a not-yet-validated family
      charge_redox : by-charge-class RMSE (V) from the OROP benchmark (deep fallback)
      lambda       : reorg sigma (eV) = spread of our production lambda vs D3TaLES (0-2 eV subset)
      disp         : disproportionation sigma (eV) from stability_validation residuals
    """
    sc = _cfg("scorecard_config")
    family_redox, pooled_redox = derive_family_sigma(redox_rows)

    # by-charge-class prior, derived from the OROP benchmark rather than hard-coded
    cls = {"cation": [], "anion": [], "multi": []}
    for r in _load("orop_benchmark"):
        e, q = _f(r.get("err")), _f(r.get("charge_ox"))
        if e is None or q is None:
            continue
        k = "multi" if (q >= 2 or q < 0) else ("cation" if q >= 1 else "anion")
        cls[k].append(e)
    charge_redox = {k: round(_rmse(v), 3) for k, v in cls.items() if v} or dict(sc.SIGMA_REDOX_V)

    # lambda: spread of our production lambda vs D3TaLES on the reliable 0-2 eV subset. Use a
    # ROBUST scale (1.4826 x MAD) not RMSE: D3TaLES's diffuse-free def2-SVP leaves a heavy
    # unbound-anion outlier tail (see reorg QC) that would otherwise inflate the sigma ~5x.
    comp = RESULTS / "d3tales_reorg_validation_prod" / "comparison.csv"
    lam_res = []
    if comp.exists():
        with comp.open() as f:
            for r in csv.DictReader(f):
                a, b = _f(r.get("our_smd_eV")), _f(r.get("d3tales_eV"))
                if a is not None and b is not None and 0 <= a <= 2 and 0 <= b <= 2:
                    lam_res.append(a - b)
    sigma_lambda = round(_mad_sigma(lam_res), 3) if lam_res else sc.SIGMA_LAMBDA_EV

    # disproportionation: stability_validation residuals (kJ/mol -> eV)
    KJ_EV = 1.0 / 96.485
    disp_res = [_f(r.get("err_kJmol")) for r in _load("stability_validation")
                if (r.get("tier") or "A") == "A"]          # grounded tier-A spacings only
    disp_res = [x * KJ_EV for x in disp_res if x is not None]
    sigma_disp = round(_rmse(disp_res), 3) if disp_res else sc.SIGMA_DISP_EV

    return dict(family_redox=family_redox, pooled_redox=pooled_redox, charge_redox=charge_redox,
                sigma_lambda=sigma_lambda, sigma_disp=sigma_disp)


def build():
    ele = _cfg("electrolyte")
    sc = _cfg("scorecard_config")
    wlo, whi = ele.WINDOW_V_VS_FC
    divider = ele.ANOLYTE_CATHOLYTE_DIVIDER_V

    redox = _load("redox_potentials")
    unc = derive_uncertainties(redox)   # all sigmas computed from benchmarks (not hard-coded)
    # lambda comes with a QC verdict from redox.properties.reorg (negative_lambda / negative_half /
    # conformer_jump / anion_unbound / large_half>1eV). A flagged value is NOT a usable
    # inner-sphere lambda, so it must not silently enter the mean or the Pareto ranking.
    reorg_rows = {(r["id"], r["couple"]): r for r in _load("reorganization")}
    integ = {(r["id"], r["couple"]): r["verdict"] for r in _load("state_integrity")}
    disp = {(r["id"], r["intermediate"]): _f(r["dG_disp_kJmol"])
            for r in _load("stability_disproportionation")}
    cap = {r["id"]: r for r in _load("capacity_and_proxies")}
    # outer-sphere lambda_o (molecular-cavity PCM + Born, 1-body), keyed by state pair
    lam_o = {(r["id"], _states_of(r["couple"])): (_f(r["lambda_o_pcm_eV"]),
                                                  _f(r["lambda_o_born_eV"]),
                                                  r.get("method", ""))
             for r in _load("lambda_outer")}
    pekar = 1.0 / ele.SOLVENT["eps_optical"] - 1.0 / ele.SOLVENT["eps_r"]
    m_anion = _ion_mass(getattr(sc, "COUNTERION_ANION_SMILES", None))
    m_cations = {k: _ion_mass(s) for k, s in getattr(sc, "SUPPORTING_CATIONS", {}).items()}
    candidate_ids = _candidate_ids()
    batches, chem_flags = _candidate_batches(), _candidate_chem_flags()
    unrankable = _unrankable_ids()
    resting = _resting_states()

    # couples by molecule (scoped to the current screening candidates)
    by_mol = {}
    for r in redox:
        if r["id"] not in candidate_ids or r["id"] in REFERENCE_IDS:
            continue
        e = _f(r.get("E_vs_Fc_V")); qo = int(r["q_ox"]); qr = int(r["q_red"])
        couple = r["event"]
        sO, sR = couple.split("->")
        rr = reorg_rows.get((r["id"], couple), {})
        by_mol.setdefault(r["id"], {"name": r["name"], "family": r["family"], "couples": {}})
        by_mol[r["id"]]["couples"][qo] = dict(
            couple=couple, sO=sO, sR=sR, E=e, q_ox=qo, q_red=qr,
            status=(r.get("status") or "ok"), thermal_qc=r.get("thermal_qc", ""),
            integrity=integ.get((r["id"], couple)),          # None = never assessed
            in_window=(e is not None and wlo <= e <= whi),
            lam=_f(rr.get("lambda_i_eV")), lam_flag=(rr.get("flag") or "").strip(),
            lam_ox=_f(rr.get("relax_ox_meV")), lam_red=_f(rr.get("relax_red_meV")),
            lam_o=lam_o.get((r["id"], frozenset((sO, sR)))),
            sigma_E=sc.sigma_redox(qo, qr, family=r["family"],
                                   family_sigma=unc["family_redox"],
                                   pooled_sigma=unc["pooled_redox"]),
        )

    rows = []
    for gid, m in by_mol.items():
        base = dict(id=gid, name=m["name"], family=m["family"], batch=batches.get(gid, ""),
                    chem_flag=chem_flags.get(gid, ""), rankable=(gid not in unrankable))
        rest = resting.get(gid)
        if rest is None:
            rows.append(dict(base, status="INCOMPLETE", reason="no declared resting state (n_e=0)",
                             n_accessible=0)); continue
        # CONTIGUOUS redox paths from the resting state (charge q0): an electron counts only
        # if every earlier step on the same side is accessible. Anolyte = successive
        # reductions below the divider; catholyte = successive oxidations above it. Each
        # electrode pool gets ONLY its own side's electrons (an ambipolar molecule is two
        # separate pool rows, never double-counted).
        paths, stop = contiguous_paths(m["couples"], rest[1], divider)
        incomplete = [v for v in stop.values() if v.startswith("INCOMPLETE")]
        if not paths["anolyte"] and not paths["catholyte"]:
            rows.append(dict(base, status=("INCOMPLETE" if incomplete else "REJECTED"),
                             reason="; ".join(stop.values()), n_accessible=0))
            continue
        c_ = cap.get(gid, {})
        mw = _f(c_.get("MW")); mw_rep = _f(c_.get("MW_repeat_unit"))
        spin_gap, spin_src = _min_spin_gap_eV(gid)
        for pool, acc in paths.items():
            if not acc:
                continue
            n_acc = len(acc)
            Es = [c["E"] for c in acc]
            # --- capacity (bookkeeping on the contiguous path; counter-ion variants are
            # scenario conventions, not exact physics: whether supporting ions count in a
            # material-level gravimetric capacity depends on the device mass convention) ---
            def _cap(mass):
                return round(n_acc * ele.FARADAY / (mass * 3.6), 1) if mass else None
            cap_maxload = {}
            for salt, m_cat in m_cations.items():
                m_ion = _max_counterion_mass(acc, m_anion, m_cat)
                cap_maxload[salt] = _cap(mw_rep + m_ion) if (mw_rep and m_ion is not None) else None
            # --- inner-sphere lambda, explicit conventions (QC-clean couples only) ---
            clean = [c for c in acc if c["lam"] is not None and not c["lam_flag"]]
            flags_seen = sorted({c["lam_flag"] for c in acc if c["lam_flag"]})
            n_lam_flagged = sum(1 for c in acc if c["lam"] is not None and c["lam_flag"])
            n_lam_missing = sum(1 for c in acc if c["lam"] is None)
            lam_qc = ("ok" if len(clean) == n_acc else ("all_flagged" if not clean else "partial"))
            mean = (lambda xs: round(statistics.mean(xs), 4) if xs else None)
            lam4 = mean([c["lam"] for c in clean])                     # lambda_O + lambda_R
            lam_O = mean([c["lam_ox"] / 1000.0 for c in clean if c["lam_ox"] is not None])
            lam_R = mean([c["lam_red"] / 1000.0 for c in clean if c["lam_red"] is not None])
            # --- outer-sphere + combined, per ET process (see module docstring) ---
            lo = [c["lam_o"][0] for c in clean if c["lam_o"] and c["lam_o"][0]]
            lam_o1 = mean(lo)
            lam_het = round(lam4 / 2.0 + lam_o1, 4) if (lam4 is not None and lam_o1) else None
            lam_se = None
            if lam4 is not None and lam_o1:
                a_eff = 7.19982 * pekar / lam_o1          # Born-equivalent radius, A
                lam_se = round(lam4 + lambda_o_self_exchange(lam_o1, a_eff, 2.0 * a_eff), 4)
            dvals = [disp[(gid, c["sR"])] for c in acc[:-1] if (gid, c["sR"]) in disp] \
                if pool == "anolyte" else \
                [disp[(gid, c["sO"])] for c in acc[:-1] if (gid, c["sO"]) in disp]
            rows.append(dict(
                base, status="candidate", role=pool, pool=pool,
                resting_state=rest[0], path=" ; ".join(c["couple"] for c in acc),
                path_stop=stop.get(pool, ""),
                n_accessible=n_acc,
                min_spin_gap_eV=spin_gap, spin_gap_source=spin_src,
                spin_confidence=("low" if (spin_gap is not None and spin_gap < sc.SPIN_GAP_LOW_EV)
                                 else "ok"),
                thermal_qc=("imag_modes" if any(c["thermal_qc"] == "imag_modes" for c in acc)
                            else "ok"),
                E_V=round(statistics.mean(Es), 3),
                E_anolyte_V=(round(statistics.mean(Es), 3) if pool == "anolyte" else None),
                E_catholyte_V=(round(statistics.mean(Es), 3) if pool == "catholyte" else None),
                sigma_E_V=round(max(c["sigma_E"] for c in acc), 3),
                specific_capacity_mAh_g=_cap(mw), MW=(round(mw, 1) if mw else None),
                MW_repeat_unit=(round(mw_rep, 1) if mw_rep else None),
                capacity_repeat_mAh_g=_cap(mw_rep),
                capacity_maxload_Li_mAh_g=cap_maxload.get("Li"),
                capacity_maxload_TBA_mAh_g=cap_maxload.get("TBA"),
                lambda_i_eV=lam4, lambda_i_ox_eV=lam_O, lambda_i_red_eV=lam_R,
                sigma_lambda_eV=(unc["sigma_lambda"] if lam4 is not None else None),
                lambda_qc=lam_qc, n_lambda_flagged=n_lam_flagged, n_lambda_missing=n_lam_missing,
                lambda_flags=";".join(flags_seen),
                lambda_o_pcm_eV=lam_o1,
                lambda_o_method=";".join(sorted({c["lam_o"][2] for c in clean if c["lam_o"]})),
                lambda_het_eV=lam_het,
                sigma_lambda_het_eV=(round(((unc["sigma_lambda"] / 2) ** 2
                                            + getattr(sc, "SIGMA_LAMBDA_O_EV", 0.0) ** 2) ** 0.5, 3)
                                     if lam_het is not None else None),
                lambda_se_contact_eV=lam_se,
                dG_disp_kJmol=(round(min(dvals), 1) if dvals else None),
                disp_applicable=(n_acc >= 2),
                sigma_disp_eV=(unc["sigma_disp"] if dvals else None),
                all_intact_bound=True,
                SA_score=_f(c_.get("SA_score")),
                dGsolv_proxy_eV=_f(c_.get("dGsolv_neutral_eV")),
            ))

    # --- write ---
    RESULTS.mkdir(exist_ok=True)
    cols = ["id", "name", "family", "batch", "chem_flag", "rankable", "status", "role", "pool",
            "resting_state", "path", "path_stop", "n_accessible", "E_V", "E_anolyte_V",
            "E_catholyte_V", "sigma_E_V", "specific_capacity_mAh_g", "MW",
            "MW_repeat_unit", "capacity_repeat_mAh_g", "capacity_maxload_Li_mAh_g",
            "capacity_maxload_TBA_mAh_g", "lambda_i_eV", "lambda_i_ox_eV", "lambda_i_red_eV",
            "sigma_lambda_eV", "lambda_qc", "n_lambda_flagged", "n_lambda_missing",
            "lambda_flags", "lambda_o_pcm_eV", "lambda_o_method", "lambda_het_eV",
            "sigma_lambda_het_eV", "lambda_se_contact_eV", "dG_disp_kJmol", "disp_applicable",
            "sigma_disp_eV", "all_intact_bound", "thermal_qc", "min_spin_gap_eV",
            "spin_gap_source", "spin_confidence", "SA_score", "dGsolv_proxy_eV", "reason"]
    out = RESULTS / "scorecard.csv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})

    cands = [r for r in rows if r["status"] == "candidate"]

    def _sortkey(r):
        return (r["role"], r["E_V"])

    print(f"{'id':22s} {'role':9s} {'n':>2s} {'E_an':>6s} {'E_cat':>6s} {'Cap':>5s} "
          f"{'lam':>5s} {'dGdisp':>7s} {'SA':>4s}")
    print("-" * 76)
    for r in sorted(cands, key=_sortkey):
        ea = f"{r['E_anolyte_V']:+.2f}" if r["E_anolyte_V"] is not None else "-"
        ec = f"{r['E_catholyte_V']:+.2f}" if r["E_catholyte_V"] is not None else "-"
        cap_s = f"{r['specific_capacity_mAh_g']:.0f}" if r["specific_capacity_mAh_g"] else "-"
        lam_s = f"{r['lambda_i_eV']:.2f}" if r["lambda_i_eV"] is not None else "-"
        if r.get("lambda_qc") and r["lambda_qc"] != "ok":
            lam_s += "!"
        dg_s = f"{r['dG_disp_kJmol']:.0f}" if r["dG_disp_kJmol"] is not None else "-"
        sa_s = f"{r['SA_score']:.1f}" if r["SA_score"] is not None else "-"
        spin_flag = "  <-- spin near-degenerate (low confidence)" if r.get("spin_confidence") == "low" else ""
        print(f"{r['id']:22s} {r['role']:9s} {r['n_accessible']:2d} "
              f"{ea:>6s} {ec:>6s} {cap_s:>5s} {lam_s:>5s} {dg_s:>7s} {sa_s:>4s}{spin_flag}")
    rej = [r for r in rows if r["status"] == "REJECTED"]
    inc = [r for r in rows if r["status"] == "INCOMPLETE"]
    for r in rej:
        print(f"REJECTED   {r['id']}: {r.get('reason','')}")
    for r in inc:
        print(f"INCOMPLETE {r['id']}: {r.get('reason','')}")
    print(f"\n{len(cands)} candidate pool-rows, {len(rej)} rejected, {len(inc)} incomplete. "
          f"sigma + trust per axis attached "
          f"(config/scorecard_config.py). Window {wlo}..{whi} V vs Fc.")
    print(f"wrote {out}")
    return rows


if __name__ == "__main__":
    build()
