#!/usr/bin/env python
"""Mine the D3TaLES dump for NEW candidates with low reorganization energy AND low synthetic
accessibility (SA) — the two "free" screening axes we can read straight from D3TaLES.

Design choices:
  - Per family we rank on the RELEVANT reorg column: electron_reorganization_energy for the
    reductive (n-type) families (quinone, pyridinium, viologen/dication) and
    hole_reorganization_energy for the oxidative (p-type) families (phenothiazine, nitroxide).
    D3TaLES's electron column is noisy (FINDINGS #10: negatives, >3 eV outliers), so reductive
    hits are flagged reorg_confidence='low' and must be re-checked with our pipeline.
  - Sane-value filter: relevant reorg in (0.05, 0.60) eV, sa_score in (1.0, 3.5).
  - NO scalarized figure-of-merit (project preference): we mark the Pareto-nondominated set on
    (reorg, SA) both-minimized and sort by reorg (the primary axis); read both raw axes.
  - Molecules already in library/manifest.csv are excluded (these are the OTHER, missing ones).
  - D3TaLES potentials are ABSOLUTE (their +4.42 V SHE add); relabel_potentials tags them and
    adds an approximate vs-Fc estimate, so a raw value never sits unlabeled next to our numbers.

  python scripts/mine_d3tales_candidates.py            # writes results/d3tales_mined_*.csv
"""
from __future__ import annotations
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from redox.d3tales_ingest import relabel_potentials  # noqa: E402

CSV = ROOT / "data" / "raw" / "validation" / "D3TaLES" / "d3tales_public.csv"
OUT = ROOT / "results"

REORG_LO, REORG_HI = 0.05, 0.60      # eV, sane inner-sphere window
SA_LO, SA_HI = 1.0, 3.5              # Ertl SA (1 easy .. 10 hard); <3.5 = readily accessible

# family -> (SMARTS list, redox direction). direction picks the relevant reorg column.
# SMARTS tightened so each matches the actual redox core, not look-alikes:
#   pyridinium   : aromatic N+ 6-ring bearing an exocyclic C (N-alkyl/aryl pyridinium),
#                  explicitly NOT a pyridine-N-oxide ([n+]-[O-]).
#   nitroxide    : genuine aminoxyl radical R2N-O* (trivalent neutral N, two C, terminal
#                  single-bonded O), which excludes oximes (C=N-O) and nitro ([N+](=O)[O-]).
FAMILIES = {
    "quinone":       (["O=C1C=CC(=O)C=C1", "O=C1c2ccccc2C(=O)c2ccccc21"], "reductive"),
    "pyridinium":    (["[n+;!$([n+][O-])]1(ccccc1)[#6]"], "reductive"),
    "phenothiazine": (["c1ccc2c(c1)[#7]c1ccccc1[#16]2"], "oxidative"),
    "nitroxide":     (["[NX3;!+;!$([N]=*)]([#6])([#6])[OX1;!-;!$([O]=*)]"], "oxidative"),
}


def _existing_canon():
    from rdkit import Chem
    seen = set()
    mf = ROOT / "library" / "manifest.csv"
    if mf.exists():
        for r in csv.DictReader(mf.open()):
            m = Chem.MolFromSmiles(r["smiles"]) if r.get("smiles") else None
            if m:
                seen.add(Chem.MolToSmiles(m))
    return seen


def _pareto_min(df, ax, ay):
    """Boolean mask of rows nondominated when MINIMIZING both ax and ay (no ties handling
    needed for a flag). O(n^2) is fine for a few-thousand-row family slice."""
    import numpy as np
    a = df[[ax, ay]].to_numpy(float)
    n = len(a)
    keep = np.ones(n, dtype=bool)
    for i in range(n):
        if not keep[i]:
            continue
        # dominated if some j is <= on both and < on at least one
        dom = ((a[:, 0] <= a[i, 0]) & (a[:, 1] <= a[i, 1]) &
               ((a[:, 0] < a[i, 0]) | (a[:, 1] < a[i, 1])))
        if dom.any():
            keep[i] = False
    return keep


def main():
    import pandas as pd
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    if not CSV.exists():
        sys.exit(f"missing {CSV} — download via docs/DATASETS.md (D3TaLES section)")

    df = pd.read_csv(CSV, low_memory=False)
    df["_mol"] = [Chem.MolFromSmiles(s) if isinstance(s, str) else None for s in df["smiles"]]
    df = df[df["_mol"].notna()].copy()
    df["_canon"] = [Chem.MolToSmiles(m) for m in df["_mol"]]
    df["sa_score"] = pd.to_numeric(df["sa_score"], errors="coerce")
    df["hole_reorganization_energy"] = pd.to_numeric(df["hole_reorganization_energy"], errors="coerce")
    df["electron_reorganization_energy"] = pd.to_numeric(df["electron_reorganization_energy"], errors="coerce")

    existing = _existing_canon()
    print(f"loaded {len(df)} parseable molecules; excluding {len(existing)} already in library")

    def viologen_mask():   # 4,4'-bipyridinium == the >=2 aromatic-N+ dication pool
        return df["smiles"].astype(str).str.count(r"\[n\+\]") >= 2

    patt = {f: [Chem.MolFromSmarts(p) for p in pats] for f, (pats, _) in FAMILIES.items()}
    keep_cols = ["smiles", "source_group", "groundState_charge", "molecular_weight",
                 "sa_score", "reorg_eV", "reorg_type", "reorg_confidence",
                 "solv_oxidation_potential", "solv_reduction_potential", "pareto_optimal"]

    all_hits = []
    for fam, (pats, direction) in list(FAMILIES.items()) + [("viologen", (None, "reductive"))]:
        if fam == "viologen":
            mask = viologen_mask()
        else:
            ps = patt[fam]
            mask = df["_mol"].apply(lambda m: any(m.HasSubstructMatch(p) for p in ps))
        col = "electron_reorganization_energy" if direction == "reductive" else "hole_reorganization_energy"
        conf = "low" if direction == "reductive" else "ok"   # electron column is noisy
        sub = df[mask].copy()
        sub["reorg_eV"] = sub[col]
        sub["reorg_type"] = "electron" if direction == "reductive" else "hole"
        sub["reorg_confidence"] = conf
        # exclude what we already have; apply sane windows
        sub = sub[~sub["_canon"].isin(existing)]
        sub = sub[sub["reorg_eV"].between(REORG_LO, REORG_HI) &
                  sub["sa_score"].between(SA_LO, SA_HI)]
        sub = sub.drop_duplicates("_canon")
        if sub.empty:
            print(f"  {fam:14s} 0 hits"); continue
        sub["pareto_optimal"] = _pareto_min(sub, "reorg_eV", "sa_score")
        sub = relabel_potentials(sub)
        sub = sub.sort_values("reorg_eV")
        outcols = [c for c in keep_cols if c in sub.columns] + \
                  [c for c in sub.columns if c.endswith("_vs_Fc_est_V") or c.endswith("_abs_V")]
        outcols = list(dict.fromkeys(outcols))
        sub[outcols].to_csv(OUT / f"d3tales_mined_{fam}.csv", index=False)
        npf = int(sub["pareto_optimal"].sum())
        print(f"  {fam:14s} {len(sub):4d} hits ({npf} Pareto)  reorg {sub['reorg_eV'].min():.3f}"
              f"..{sub['reorg_eV'].max():.3f} eV [{conf}]  -> results/d3tales_mined_{fam}.csv")
        top = sub.head(5)
        for _, r in top.iterrows():
            print(f"      reorg={r['reorg_eV']:.3f} SA={r['sa_score']:.2f} MW={r['molecular_weight']:.0f}"
                  f"  {'PARETO ' if r['pareto_optimal'] else ''}{r['smiles'][:60]}")
        sub["family"] = fam
        all_hits.append(sub[outcols + ["family", "_canon"]])

    if all_hits:
        allh = pd.concat(all_hits).drop_duplicates("_canon").drop(columns="_canon")
        allh = allh.sort_values(["reorg_confidence", "reorg_eV"])  # 'ok' (hole) first, then low
        allh.to_csv(OUT / "d3tales_mined_all.csv", index=False)
        print(f"\ntotal unique mined candidates: {len(allh)}  -> results/d3tales_mined_all.csv")
        print("NOTE: reductive (electron-reorg) hits are reorg_confidence='low' (D3TaLES electron "
              "column is noisy) — re-verify shortlisted ones with our SMD pipeline before trusting.")


if __name__ == "__main__":
    main()
