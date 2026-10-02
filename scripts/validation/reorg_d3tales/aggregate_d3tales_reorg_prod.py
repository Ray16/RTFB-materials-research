#!/usr/bin/env python
"""Aggregate the Track-2 D3TaLES reorganization sweep (our production SMD pipeline) into a
comparison table + summary stats.

Merges three electron-couple lambda values per molecule:
  - our_smd     : our PRODUCTION pipeline (SMD-opt geoms, r2SCAN/wB97M-V, diffuse)  [prod calc/]
  - our_d3level : D3TaLES's TRUE level reproduced (lc_wpbe + their tuned omega + def2-svp, gas)
                  [d3level calc/]
  - d3tales     : D3TaLES's own reported electron_reorganization_energy

Writes results/d3tales_reorg_validation_prod/comparison.csv and prints:
  - our_d3level vs d3tales   THE MATCH TEST: at their exact level, do we reproduce them?
  - our_smd     vs d3tales   how our (better, solvated, diffuse) pipeline compares
  - our_smd     vs our_d3level  the functional/phase shift between our level and theirs

  PYTHONPATH=src python scripts/validation/reorg_d3tales/aggregate_d3tales_reorg_prod.py
"""
from __future__ import annotations
import csv
import json
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
PROD = ROOT / "results" / "d3tales_reorg_validation_prod" / "calc"
D3LEVEL = ROOT / "results" / "d3tales_reorg_validation_d3level" / "calc"
OUT = ROOT / "results" / "d3tales_reorg_validation_prod" / "comparison.csv"


def _num(x):
    try:
        v = float(x)
        return v if v == v else None   # drop NaN
    except (TypeError, ValueError):
        return None


def main():
    rows = []
    for p in sorted(PROD.glob("*.json")):
        d = json.loads(p.read_text())
        if d.get("status") not in ("ok", "partial"):
            continue
        gid = d["id"]
        our_smd = _num(d.get("our_electron"))
        d3 = _num(d.get("d3_electron"))
        # D3TaLES-level (lc_wpbe + their omega + def2-svp gas) electron lambda for this molecule
        our_d3level = None
        dlp = D3LEVEL / f"{gid}.json"
        if dlp.exists():
            try:
                dl = json.loads(dlp.read_text())
                if dl.get("status") in ("ok", "partial"):
                    our_d3level = _num(dl.get("our_electron"))
            except (OSError, ValueError):
                pass
        rows.append(dict(id=gid, smiles=d.get("smiles", ""), family=d.get("family", ""),
                         n_atoms=d.get("n_atoms", ""),
                         our_smd_eV=our_smd, our_d3level_eV=our_d3level, d3tales_eV=d3,
                         relax_neu_meV=d.get("relax_neu_meV"), relax_anion_meV=d.get("relax_anion_meV")))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    cols = ["id", "smiles", "family", "n_atoms", "our_smd_eV", "our_d3level_eV", "d3tales_eV",
            "relax_neu_meV", "relax_anion_meV"]
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)

    def _stats(pairs, label):
        # keep only sane finite pairs in (0, 2) eV (drop D3TaLES garbage / broken geoms)
        good = [(a, b) for a, b in pairs if a is not None and b is not None
                and 0 < a < 2 and 0 < b < 2]
        if not good:
            print(f"  {label}: no comparable pairs"); return
        diffs = [abs(a - b) for a, b in good]
        signed = [a - b for a, b in good]
        print(f"  {label}: n={len(good):3d}  MAD={statistics.mean(diffs):.3f} eV  "
              f"median={statistics.median(diffs):.3f}  mean_signed={statistics.mean(signed):+.3f}")

    print(f"aggregated {len(rows)} finished molecules -> {OUT}\n")
    print("Electron-couple lambda comparisons (sane 0-2 eV subset):")
    _stats([(r["our_d3level_eV"], r["d3tales_eV"]) for r in rows], "our_D3TaLES-level (lc_wpbe)  vs  D3TaLES   [THE MATCH TEST]")
    _stats([(r["our_smd_eV"], r["d3tales_eV"]) for r in rows], "our_SMD (production)         vs  D3TaLES")
    _stats([(r["our_smd_eV"], r["our_d3level_eV"]) for r in rows], "our_SMD  vs  our_D3TaLES-level   [functional/phase shift]")


if __name__ == "__main__":
    main()
