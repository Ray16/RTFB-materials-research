"""Thermodynamic stability of redox intermediates from decomposition free energies.

Selecting a good flow-battery redox molecule needs BOTH a target potential AND a durable
charged/radical state. Durability is dominated by decomposition *reactions* of the reactive
intermediate, so we score it the same way we score everything else: with a Delta-G.

Reactions scored (per reactive intermediate R^q, the state that actually cycles):

  (1) DISPROPORTIONATION   2 R^q -> R^(q+1) + R^(q-1)
      dG_disp = [G(R^q+1) + G(R^q-1)] - 2 G(R^q)
      dG_disp > 0  => R^q is stable against disproportionation (radical persists).
      This is EXACT from the redox-state energies we already have (no new calc) and is
      atom- AND electron-conserving, so systematic DFT/solvation errors cancel strongly
      (much better cancellation than an absolute redox potential). Equivalent to the wave
      spacing: dG_disp = F*(E_high - E_low); K_comp = exp(+dG_disp/RT).

  (2) DIMERIZATION         2 R^q -> (R2)^(2q)          [see dimerize.py / stability_dimer]
      dG_dim = G(dimer) - 2 G(R^q).  Needs a built dimer geometry; handled separately.
      dG_dim > 0 => stable against dimerization. (THE viologen failure mode.)

Sign convention everywhere: dG_decomp > 0 == STABLE (decomposition uphill).
Free energy per species: G = e_smd_eV + g_thermal_eV (GFN2-xTB RRHO), our standard.

  PYTHONPATH=src python -m redox.stability            # disproportionation table
"""
from __future__ import annotations
import csv
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DFT = ROOT / "calcs" / "dft"
RESULTS = ROOT / "results"

EV_KJ = 96.485          # eV -> kJ/mol
KT_EV = 0.0256926       # k_B * T at 298.15 K, in eV
F_EV_PER_V = 1.0        # 1 electron: dG(eV) per volt is numerically 1


def G(gid, state):
    """Free energy (eV) = e_smd + g_thermal for one species state; None if missing."""
    p = DFT / gid / state / "result.json"
    if not p.exists():
        return None
    r = json.loads(p.read_text())
    e = r.get("e_smd_eV")
    if e is None:
        return None
    return e + (r.get("g_thermal_eV") or 0.0)


def _states_by_charge(gid):
    """Return [(state, charge)] present in calcs/dft for gid, sorted high->low charge,
    using the charge recorded in each result.json (authoritative)."""
    out = []
    d = DFT / gid
    if not d.exists():
        return out
    for sd in d.iterdir():
        rj = sd / "result.json"
        if rj.exists():
            try:
                q = json.loads(rj.read_text()).get("charge")
                if q is not None:
                    out.append((sd.name, int(q)))
            except Exception:
                pass
    return sorted(out, key=lambda x: -x[1])


def disproportionation(gid):
    """dG_disp for every interior state R^q that has BOTH neighbors R^(q+1), R^(q-1).

    Returns list of dicts. The interior state is the reactive intermediate (radical /
    semiquinone / radical cation) whose persistence we care about.
    """
    states = _states_by_charge(gid)          # high -> low charge, adjacent differ by 1e-
    by_q = {q: s for s, q in states}
    rows = []
    for s, q in states:
        if (q + 1) in by_q and (q - 1) in by_q:
            g_mid = G(gid, s)
            g_hi = G(gid, by_q[q + 1])
            g_lo = G(gid, by_q[q - 1])
            if None in (g_mid, g_hi, g_lo):
                continue
            dG = (g_hi + g_lo) - 2.0 * g_mid          # eV
            rows.append(dict(
                id=gid, intermediate=s, q=q,
                reaction=f"2 {s}(q={q:+d}) -> {by_q[q+1]}(q={q+1:+d}) + {by_q[q-1]}(q={q-1:+d})",
                dG_disp_eV=round(dG, 4),
                dG_disp_kJmol=round(dG * EV_KJ, 1),
                logK_disp=round(-dG / (KT_EV * math.log(10)), 2),   # log10 K_disp
                stable_vs_disprop=bool(dG > 0),
            ))
    return rows


def main():
    # every group in calcs/dft that has >=3 charge states can have an interior intermediate
    gids = sorted({p.parent.parent.name for p in DFT.glob("*/*/result.json")})
    all_rows = []
    for gid in gids:
        all_rows.extend(disproportionation(gid))

    if not all_rows:
        print("no interior redox states found (need 3 consecutive charge states)")
        return

    hdr = f"{'id':22s} {'interm':6s} {'q':>3s} {'dG_disp(kJ/mol)':>15s} {'log10K':>7s} {'stable?':>8s}"
    print(hdr); print("-" * len(hdr))
    for r in sorted(all_rows, key=lambda x: -x["dG_disp_kJmol"]):
        print(f"{r['id']:22s} {r['intermediate']:6s} {r['q']:+3d} "
              f"{r['dG_disp_kJmol']:15.1f} {r['logK_disp']:7.2f} "
              f"{'YES' if r['stable_vs_disprop'] else 'no':>8s}")
    print("\nSign convention: dG_disp > 0  =>  intermediate STABLE against disproportionation.")
    print("(This is exact from existing redox energies; atom+electron conserving => errors cancel.)")

    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / "stability_disproportionation.csv"
    cols = ["id", "intermediate", "q", "reaction", "dG_disp_eV", "dG_disp_kJmol",
            "logK_disp", "stable_vs_disprop"]
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for r in all_rows:
            w.writerow(r)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
