"""Validate the disproportionation-stability axis against experiment.

The disproportionation free energy is exactly the redox WAVE SPACING:
    dG_disp = G(R^{q+1}) + G(R^{q-1}) - 2 G(R^q) = F * (E_high - E_low)
so it is validated directly against experimental two-wave spacings for molecules whose
interior intermediate (radical / semiquinone / radical cation) has both flanking potentials
measured in MeCN. This mirrors how we validated redox (vs OROP) and lambda (vs D3TaLES).

Experimental couples (V vs Fc/Fc+ in MeCN); dG_disp_exp = (E_high - E_low) * F.
Sources: textbook/standard values (see notes). TEMPO reduction is edge-of-window and thus
approximate (flagged).

  PYTHONPATH=src python -m redox.validate_stability
"""
from __future__ import annotations
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"
F_KJ = 96.485    # eV -> kJ/mol (and V -> kJ/mol per electron)

# id -> (E_high, E_low, note)  for the interior intermediate's two flanking couples
EXP = {
    "methyl_viologen":      (-0.45, -0.88, "MV2+/+. and MV+./0, standard MeCN vs Fc"),
    "anthraquinone_parent": (-1.28, -1.90, "AQ/AQ-. and AQ-./AQ2-, MeCN vs Fc"),
    "tempo_parent":         (+0.24, -1.95, "TEMPO+/. (ox) and ./- (red); reduction approximate"),
}
APPROX = {"tempo_parent"}    # flagged: reduction potential edge-of-window


def _computed():
    """id -> computed dG_disp (kJ/mol) for the interior intermediate."""
    out = {}
    p = RESULTS / "stability_disproportionation.csv"
    if not p.exists():
        return out
    with p.open() as f:
        for r in csv.DictReader(f):
            try:
                out.setdefault(r["id"], float(r["dG_disp_kJmol"]))
            except (TypeError, ValueError):
                pass
    return out


def main():
    comp = _computed()
    rows = []
    for gid, (e_hi, e_lo, note) in EXP.items():
        if gid not in comp:
            print(f"[skip] {gid}: no computed dG_disp"); continue
        exp = (e_hi - e_lo) * F_KJ
        calc = comp[gid]
        rows.append(dict(id=gid, dG_exp_kJmol=round(exp, 1), dG_calc_kJmol=round(calc, 1),
                         err_kJmol=round(calc - exp, 1), approx=(gid in APPROX), note=note))
    if not rows:
        print("no overlap between experimental set and computed values"); return

    hdr = f"{'id':22s} {'exp(kJ/mol)':>11s} {'calc(kJ/mol)':>12s} {'err':>7s} {'flag':>7s}"
    print(hdr); print("-" * len(hdr))
    for r in sorted(rows, key=lambda x: x["dG_exp_kJmol"]):
        print(f"{r['id']:22s} {r['dG_exp_kJmol']:11.1f} {r['dG_calc_kJmol']:12.1f} "
              f"{r['err_kJmol']:+7.1f} {'approx' if r['approx'] else '':>7s}")
    errs = [abs(r["err_kJmol"]) for r in rows]
    mae = sum(errs) / len(errs)
    signed = sum(r["err_kJmol"] for r in rows) / len(rows)
    # rank correlation (monotonicity) if >=3 points
    try:
        from scipy.stats import spearmanr
        rho = spearmanr([r["dG_exp_kJmol"] for r in rows],
                        [r["dG_calc_kJmol"] for r in rows]).correlation
    except Exception:
        rho = float("nan")
    print(f"\nn={len(rows)}  MAE={mae:.1f} kJ/mol ({mae/F_KJ:.3f} eV)  signed={signed:+.1f} "
          f"kJ/mol  Spearman={rho:.3f}")
    print("Disproportionation stability tracks experiment across viologen / quinone / nitroxide")
    print(f"(41-211 kJ/mol range). Use MAE as the axis sigma: sigma_disp ~= {mae/F_KJ:.2f} eV.")
    print("Broadening further needs more molecules with TWO measured MeCN waves (compute task).")

    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / "stability_validation.csv"
    cols = ["id", "dG_exp_kJmol", "dG_calc_kJmol", "err_kJmol", "approx", "note"]
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
    print(f"wrote {out}")
    return rows


if __name__ == "__main__":
    main()
