#!/usr/bin/env python
"""How can the D3TaLES dump (35,729 molecules, implicit-MeCN DFT) serve THIS project?

Answers three concrete questions and writes small, committable artifacts:

  1) COVERAGE — how many D3TaLES molecules match each of our 6 redox families
     (pyridinium, viologen/bipyridinium, quinone, phenothiazine, nitroxide/TEMPO,
     metallocene)? -> results/d3tales_family_coverage.csv
  2) CANDIDATE POOL — for each family, the most synthetically-accessible members
     (low sa_score) with their computed ox/red potentials -> results/d3tales_candidates_<family>.csv
  3) DICATION CROSS-CHECK — pull the two-N+ (viologen-like) subset; these are exactly
     the species D3TaLES's own README flags as implicit-solvent-unreliable (+2/-2),
     independently corroborating our viologen diagnosis. -> results/d3tales_dications.csv

Potentials are ABSOLUTE (eV, vs the free electron); subtract a level-matched Fc
absolute to get V vs Fc/Fc+. We only RANK within D3TaLES here (reference cancels), so
no conversion is applied; conversion is a separate step when comparing to our numbers.

  python scripts/analyze_d3tales.py
"""
from __future__ import annotations
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CSV = ROOT / "data" / "raw" / "validation" / "D3TaLES" / "d3tales_public.csv"
OUT = ROOT / "results"

# SMARTS for our 6 families (kept deliberately broad; we RANK/COUNT, not gate on these)
FAMILIES = {
    "pyridinium":    ["[n+]1ccccc1", "[n+]1ccccc1"],           # aromatic N+ 6-ring
    "viologen":      ["c1cc[n+]cc1-c1cc[n+]cc1", "[n+]1ccc(cc1)-c1cc[n+]cc1"],  # 4,4'-bipyridinium
    "quinone":       ["O=C1C=CC(=O)C=C1", "O=C1c2ccccc2C(=O)c2ccccc21"],        # p-quinone / AQ
    "phenothiazine": ["c1ccc2c(c1)Nc1ccccc1S2", "c1ccc2c(c1)[#7]c1ccccc1[#16]2"],
    "nitroxide":     ["[#7]([#6])([#6])[O]", "[O][N]"],        # aminoxyl / TEMPO-like
    "metallocene":   ["[Fe]", "[Ru]"],
}


def main():
    import pandas as pd
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")

    if not CSV.exists():
        sys.exit(f"missing {CSV} — download via docs/DATASETS.md (D3TaLES section)")
    df = pd.read_csv(CSV, low_memory=False)
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"loaded {len(df)} molecules, {df.shape[1]} cols")

    # pre-parse SMILES once
    mols = [Chem.MolFromSmiles(s) if isinstance(s, str) else None for s in df["smiles"]]
    patt = {fam: [Chem.MolFromSmarts(p) for p in pats if Chem.MolFromSmarts(p)]
            for fam, pats in FAMILIES.items()}

    def match(mol, pats):
        return mol is not None and any(mol.HasSubstructMatch(p) for p in pats)

    # 1) coverage
    cov_rows = []
    fam_masks = {}
    for fam, pats in patt.items():
        mask = [match(m, pats) for m in mols]
        fam_masks[fam] = mask
        n = sum(mask)
        cov_rows.append(dict(family=fam, n_matches=n,
                             pct=round(100 * n / len(df), 2)))
        print(f"  {fam:14s} {n:6d}  ({100*n/len(df):.1f}%)")
    pd.DataFrame(cov_rows).to_csv(OUT / "d3tales_family_coverage.csv", index=False)
    print(f"  -> results/d3tales_family_coverage.csv")

    # 2) per-family candidate shortlists (most synthesizable, with potentials)
    keep = ["_id", "smiles", "source_group", "groundState_charge", "molecular_weight",
            "sa_score", "solv_oxidation_potential", "solv_reduction_potential",
            "hole_reorganization_energy", "electron_reorganization_energy"]
    keep = [c for c in keep if c in df.columns]
    for fam, mask in fam_masks.items():
        sub = df[mask][keep].copy()
        if "sa_score" in sub:
            sub = sub.sort_values("sa_score")
        sub.head(50).to_csv(OUT / f"d3tales_candidates_{fam}.csv", index=False)

    # 3) dication cross-check (two aromatic N+) — README-flagged implicit-solvent-hard
    di = df[df["smiles"].astype(str).str.count(r"\[n\+\]") >= 2][keep].copy()
    di.to_csv(OUT / "d3tales_dications.csv", index=False)
    print(f"  dications (>=2 [n+]): {len(di)}  -> results/d3tales_dications.csv")

    # potential-scale sanity
    for c in ["solv_oxidation_potential", "solv_reduction_potential"]:
        s = pd.to_numeric(df[c], errors="coerce").dropna()
        print(f"  {c}: n={len(s)} range {s.min():.2f}..{s.max():.2f} eV (ABSOLUTE; "
              f"subtract Fc_abs for vs-Fc)")


if __name__ == "__main__":
    main()
