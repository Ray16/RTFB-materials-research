#!/usr/bin/env python
"""Aggregate the Track-2 D3TaLES reorganization sweep (our production SMD pipeline) into a
comparison table + summary stats.

Merges three electron-couple lambda values per molecule:
  - our_smd   : Track-2, this sweep  (SMD-opt geometries, r2SCAN/wB97M-V)   [prod calc/]
  - our_b3lyp : Track-1               (B3LYP/6-31G* gas-opt, D3TaLES protocol) [validation calc/]
  - d3tales   : D3TaLES's own reported electron_reorganization_energy

Writes results/d3tales_reorg_validation_prod/comparison.csv and prints:
  - our_smd vs d3tales           (how our full pipeline compares to their reported value)
  - our_smd vs our_b3lyp         (the SOLVENT-GEOMETRY shift: same molecule, SMD-opt vs gas-opt)

  PYTHONPATH=src python scripts/aggregate_d3tales_reorg_prod.py
"""
from __future__ import annotations
import csv
import json
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROD = ROOT / "results" / "d3tales_reorg_validation_prod" / "calc"
T1 = ROOT / "results" / "d3tales_reorg_validation" / "calc"
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
        # Track-1 (B3LYP gas) electron lambda for the same molecule, if computed
        our_b3lyp = None
        t1p = T1 / f"{gid}.json"
        if t1p.exists():
            try:
                our_b3lyp = _num(json.loads(t1p.read_text()).get("our_electron"))
            except (OSError, ValueError):
                pass
        rows.append(dict(id=gid, smiles=d.get("smiles", ""), family=d.get("family", ""),
                         n_atoms=d.get("n_atoms", ""),
                         our_smd_eV=our_smd, our_b3lyp_gas_eV=our_b3lyp, d3tales_eV=d3,
                         relax_neu_meV=d.get("relax_neu_meV"), relax_anion_meV=d.get("relax_anion_meV")))

    OUT.parent.mkdir(parents=True, exist_ok=True)
    cols = ["id", "smiles", "family", "n_atoms", "our_smd_eV", "our_b3lyp_gas_eV", "d3tales_eV",
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
    _stats([(r["our_smd_eV"], r["d3tales_eV"]) for r in rows], "our_SMD  vs  D3TaLES(reported)")
    _stats([(r["our_smd_eV"], r["our_b3lyp_gas_eV"]) for r in rows], "our_SMD  vs  our_B3LYP(gas)  [solvent-geometry shift]")
    _stats([(r["our_b3lyp_gas_eV"], r["d3tales_eV"]) for r in rows], "our_B3LYP vs  D3TaLES(reported)")


if __name__ == "__main__":
    main()
