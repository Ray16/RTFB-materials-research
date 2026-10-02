"""Validate the disproportionation-stability axis against experiment.

The disproportionation free energy is exactly the redox WAVE SPACING:
    dG_disp = G(R^{q+1}) + G(R^{q-1}) - 2 G(R^q) = F * (E_high - E_low)
so it is validated against measured two-wave spacings in MeCN for the interior intermediate
(radical anion / radical cation). A spacing is reference-free: it needs no Fc conversion.

Experimental spacings come ONLY from grounded sources (no hand-typed values):
  * config/benchmark_mecn.py "spacing" events (data/raw/validation/two_wave_mecn/, each with
    DOI + table/page; tier A = both waves E deg'/E1/2 of (quasi)reversible waves, tier B = e.g.
    half-peak potentials);
  * config/validation.py anchors carrying exp_dE12_V (e.g. methyl viologen, Cook et al. 2017
    Table 1, -759 - (-1179) mV).
Rows are matched to results/stability_disproportionation.csv by (id, intermediate state).
Only tier A enters sigma_disp (redox.screening.scorecard); tier B is reported alongside.

  PYTHONPATH=src python -m redox.validation.stability
"""
from __future__ import annotations

import csv

from redox.core.common import RESULTS, load_config

F_KJ = 96.485    # kJ mol^-1 V^-1 (Faraday constant / 1000), V of spacing -> kJ/mol


def _intermediate(states) -> str | None:
    """The interior state of a 3-state ladder (the species that can disproportionate)."""
    names = [s[0] for s in states]
    for cand in ("red1", "ox1"):
        if cand in names:
            return cand
    return None


def experimental_spacings() -> list[dict]:
    """[dict(id, intermediate, dE_V, tier, source)] from grounded sources only."""
    out = []
    for b in getattr(load_config("benchmark_mecn"), "BENCHMARK", []):
        if b.get("excluded"):
            continue
        inter = "ox1" if b["id"] == "methyl_viologen" else _intermediate(b["states"])
        for ev in b["events"]:
            if ev["kind"] == "spacing" and inter:
                out.append(dict(id=b["id"], intermediate=inter, dE_V=ev["exp_V"],
                                tier=ev["tier"], source=f"{ev['source'][:60]} doi:{ev['doi']}"))
    have = {(r["id"], r["intermediate"]) for r in out}
    for v in getattr(load_config("validation"), "VALIDATION", []):
        if "exp_dE12_V" in v:
            inter = _intermediate(v["states"])
            if (v["id"], inter) not in have:
                out.append(dict(id=v["id"], intermediate=inter, dE_V=v["exp_dE12_V"], tier="A",
                                source=v.get("exp_dE12_note", "config/validation.py")))
    return out


def _computed() -> dict:
    """(id, intermediate) -> computed dG_disp (kJ/mol)."""
    p = RESULTS / "stability_disproportionation.csv"
    out = {}
    if p.exists():
        with p.open() as f:
            for r in csv.DictReader(f):
                try:
                    out[(r["id"], r["intermediate"])] = float(r["dG_disp_kJmol"])
                except (TypeError, ValueError):
                    pass
    return out


def main():
    comp = _computed()
    rows = []
    for e in experimental_spacings():
        k = (e["id"], e["intermediate"])
        if k not in comp:
            print(f"[skip] {e['id']}/{e['intermediate']}: no computed dG_disp (INCOMPLETE)")
            continue
        exp = e["dE_V"] * F_KJ
        rows.append(dict(id=e["id"], intermediate=e["intermediate"],
                         dG_exp_kJmol=round(exp, 1), dG_calc_kJmol=round(comp[k], 1),
                         err_kJmol=round(comp[k] - exp, 1), tier=e["tier"], source=e["source"]))
    if not rows:
        print("no overlap between experimental set and computed values")
        return

    hdr = f"{'id':40s} {'tier':>4s} {'exp(kJ/mol)':>11s} {'calc(kJ/mol)':>12s} {'err':>7s}"
    print(hdr); print("-" * len(hdr))
    for r in sorted(rows, key=lambda x: x["dG_exp_kJmol"]):
        print(f"{r['id']:40s} {r['tier']:>4s} {r['dG_exp_kJmol']:11.1f} {r['dG_calc_kJmol']:12.1f} "
              f"{r['err_kJmol']:+7.1f}")
    for tier in ("A", "B"):
        rs = [r for r in rows if r["tier"] == tier]
        if not rs:
            continue
        mae = sum(abs(r["err_kJmol"]) for r in rs) / len(rs)
        signed = sum(r["err_kJmol"] for r in rs) / len(rs)
        rmse = (sum(r["err_kJmol"] ** 2 for r in rs) / len(rs)) ** 0.5
        try:
            from scipy.stats import spearmanr
            rho = spearmanr([r["dG_exp_kJmol"] for r in rs],
                            [r["dG_calc_kJmol"] for r in rs]).correlation if len(rs) >= 3 else float("nan")
        except Exception:
            rho = float("nan")
        print(f"tier {tier}: n={len(rs)}  MAE={mae:.1f}  RMSE={rmse:.1f}  signed={signed:+.1f} kJ/mol "
              f"(RMSE {rmse / F_KJ:.3f} eV)  Spearman={rho:.2f}")

    out = RESULTS / "stability_validation.csv"
    cols = ["id", "intermediate", "dG_exp_kJmol", "dG_calc_kJmol", "err_kJmol", "tier", "source"]
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
    print(f"wrote {out}")
    return rows


if __name__ == "__main__":
    main()
