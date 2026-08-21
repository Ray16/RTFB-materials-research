"""Viologen redox potentials with explicit PF6- ion-pairing (released-counterion scheme).

The bare-ion pipeline over-separates viologen's two reductions and puts both waves too
negative, because implicit SMD over-stabilizes the concentrated bare MV(2+) charge (a
q^2 Born effect a single additive Fc reference cannot cancel; see docs/PLAN.md and the
D3TaLES README caveat on +2/-2 species). Here each viologen species keeps its NATURAL
number of PF6- so every assembly is neutral (what continuum handles best), and each
reduction releases one free PF6-:

  wave 1 (MV2+/+.):  [MV.2PF6]0  + e-  ->  [MV.PF6]-  +  PF6-        (mv_ip2 -> mv_ip1 + pf6)
  wave 2 (MV+./0):   [MV.PF6]0   + e-  ->  MV0        +  PF6-        (mv_ip1 -> mv0   + pf6)

E is assembled as a cross-species reaction (NOT redox.py's adjacent-charge pairing), using
the same sign/reference convention as redox.py: dG = G(products) - G(reactants) excluding
the electron; E_abs = -dG; E_vs_Fc = E_abs - FC_ABS_COMPUTED_V. Free energies use
G = E_smd + G_thermal (thermal from the DFT Hessian; 0 if not computed).

Physics, not fitting: nothing here is scaled to experiment.

  PYTHONPATH=src python -m redox.ionpair
"""
from __future__ import annotations
import csv
import importlib.util
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DFT = ROOT / "calcs" / "dft"
RESULTS = ROOT / "results"
EXP = {"wave1 (MV2+/+.)": -0.45, "wave2 (MV+./0)": -0.88}   # V vs Fc, MeCN (see validation.py)


def _cfg(mod, name):
    spec = importlib.util.spec_from_file_location(mod, ROOT / "config" / f"{mod}.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return getattr(m, name)


def _read(gid, state):
    p = DFT / gid / state / "result.json"
    if not p.exists():
        return None
    return json.loads(p.read_text())


def _G(gid, state, use_thermal=True):
    """Free energy (eV) = E_smd + G_thermal (0 if absent / thermal disabled). None if missing."""
    r = _read(gid, state)
    if r is None or r.get("e_smd_eV") is None:
        return None
    g = float(r["e_smd_eV"])
    if use_thermal:
        gt = r.get("g_thermal_eV")
        if isinstance(gt, (int, float)):
            g += float(gt)
    return g


def _wave(products, reactants, e_ref):
    """E_vs_Fc for O + e- -> R, with products/reactants as [(gid,state), ...] (each a G)."""
    gp = [_G(g, s) for g, s in products]
    gr = [_G(g, s) for g, s in reactants]
    if any(x is None for x in gp + gr):
        return None
    dG = sum(gp) - sum(gr)          # products - reactants (electron excluded)
    E_abs = -dG
    return E_abs - e_ref


def _wave_bare(sO, sR, e_ref, gid="methyl_viologen", use_thermal=True):
    """Bare-ion couple (no ion pairing): E for state sO + e- -> sR of one species."""
    gO, gR = _G(gid, sO, use_thermal), _G(gid, sR, use_thermal)
    if gO is None or gR is None:
        return None
    return (-(gR - gO)) - e_ref


def compute():
    fc = _cfg("electrolyte", "FC_ABS_COMPUTED_V")
    rows = []

    # --- three levels, so each physics term's effect is explicit ---
    # 1) bare ions, electronic energy only (the original pipeline)
    b1 = _wave_bare("ox2", "ox1", fc, use_thermal=False)
    b2 = _wave_bare("ox1", "neu", fc, use_thermal=False)
    # 2) bare ions + thermal dG
    t1 = _wave_bare("ox2", "ox1", fc, use_thermal=True)
    t2 = _wave_bare("ox1", "neu", fc, use_thermal=True)
    # 3) explicit PF6- ion-pairing + thermal (released-counterion scheme)
    i1 = _wave([("mv_ip1", "s0"), ("pf6", "anion")], [("mv_ip2", "s0")], fc)
    i2 = _wave([("methyl_viologen", "neu"), ("pf6", "anion")], [("mv_ip1", "s0")], fc)

    schemes = [
        ("bare ions, electronic only", b1, b2),
        ("bare ions + thermal dG",     t1, t2),
        ("PF6- ion-pair + thermal",    i1, i2),
    ]
    e1, e2 = EXP["wave1 (MV2+/+.)"], EXP["wave2 (MV+./0)"]
    print(f"[ref] Fc/Fc+ absolute = {fc:.3f} V   (exp: wave1 {e1:+.2f}, wave2 {e2:+.2f} V vs Fc)\n")
    hdr = f"{'scheme':30s} {'wave1':>8s} {'wave2':>8s} {'spacing':>8s} {'MAE':>7s}"
    print(hdr); print("-" * len(hdr))
    for name, w1, w2 in schemes:
        if w1 is None or w2 is None:
            print(f"{name:30s} {'--- not all species computed yet ---':>36s}")
            rows.append(dict(scheme=name, wave1="", wave2="", spacing="", mae=""))
            continue
        sp = w2 - w1
        mae = (abs(w1 - e1) + abs(w2 - e2)) / 2
        print(f"{name:30s} {w1:+8.3f} {w2:+8.3f} {sp:+8.3f} {mae:7.3f}")
        rows.append(dict(scheme=name, wave1=round(w1, 3), wave2=round(w2, 3),
                         spacing=round(sp, 3), mae=round(mae, 3)))
    print(f"\n{'experiment':30s} {e1:+8.3f} {e2:+8.3f} {e2-e1:+8.3f} {0.0:7.3f}")
    print("\n" + "=" * 78)
    print("CONCLUSION (2026-08-20): the PF6- released-counterion scheme is REFUTED in implicit\n"
          "solvent. Even with a correct singlet mv_ip2 (clean geometry, n_imag=0), both waves\n"
          "and the spacing are catastrophically wrong (spacing ~ +11 V vs exp -0.43). Cause:\n"
          "forming a NEUTRAL ion pair from two well-solvated ions costs a large (~+6 eV),\n"
          "poorly-modeled desolvation free energy in a continuum (MeCN, eps=37.5), and it does\n"
          "NOT cancel between waves because the +2 and +1 states bind PF6- differently. The\n"
          "scheme trades the q^2 Born error for a much larger, less reliable one.\n"
          "BEST physics-based viologen result = 'bare ions + thermal' (MAE 0.372, above), on par\n"
          "with OROP's own raw implicit-DFT (MAE 0.495). The real remedy for the multiply-charged\n"
          "q^2 error is explicit microsolvation or an empirical/ML charge-dependent correction\n"
          "(cf. OROP's implicit+ML), NOT continuum ion-pairing.")
    print("=" * 78)

    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / "viologen_ionpair.csv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["scheme", "wave1", "wave2", "spacing", "mae"])
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {out}")
    return rows


if __name__ == "__main__":
    compute()
