#!/usr/bin/env python
"""Aggregate the D3TaLES reorg validation batch: our (matched-protocol) lambda_i vs D3TaLES,
across all quinone + imide molecules. Produces a comparison CSV + summary stats, and tests
whether disagreements cluster on conformationally FLEXIBLE molecules (rotatable bonds).

  python scripts/aggregate_d3tales_reorg.py
"""
from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results" / "d3tales_reorg_validation"
CALC = BASE / "calc"


def main():
    rows = []
    for f in sorted(CALC.glob("*.json")):
        try:
            rows.append(json.loads(f.read_text()))
        except Exception:
            pass
    df = pd.DataFrame(rows)
    n_all = len(df)
    ok = df[df["status"].isin(["ok", "partial"])].copy()
    fail = df[df["status"] == "fail"]
    print(f"molecules: {n_all}  ok/partial: {len(ok)}  failed: {len(fail)}")
    if len(fail):
        from collections import Counter
        errs = Counter(e.split(":")[0] for e in fail["error"].fillna(""))
        print("  failure reasons:", dict(errs))

    def block(col_our, col_d3, label):
        s = ok.dropna(subset=[col_our, col_d3]).copy()
        s = s[pd.to_numeric(s[col_d3], errors="coerce").notna()]
        s[col_our] = s[col_our].astype(float); s[col_d3] = s[col_d3].astype(float)
        d = (s[col_our] - s[col_d3])
        ad = d.abs()
        print(f"\n=== {label}  (n={len(s)}) ===")
        print(f"  MAD={ad.mean():.3f}  median|diff|={ad.median():.3f}  RMSD={np.sqrt((d**2).mean()):.3f}")
        print(f"  within 0.10 eV: {(ad<=0.10).mean()*100:.0f}%   within 0.20 eV: {(ad<=0.20).mean()*100:.0f}%")
        r = np.corrcoef(s[col_our], s[col_d3])[0, 1]
        print(f"  Pearson r={r:.3f}")
        # flexibility split: does disagreement grow with rotatable bonds?
        rigid = s[s["n_rot"] <= 1]; flex = s[s["n_rot"] >= 3]
        print(f"  rigid (n_rot<=1, n={len(rigid)}):  MAD={ (rigid[col_our]-rigid[col_d3]).abs().mean():.3f}")
        print(f"  flexible (n_rot>=3, n={len(flex)}): MAD={ (flex[col_our]-flex[col_d3]).abs().mean():.3f}")
        # biggest disagreements
        s["absdiff"] = ad
        top = s.sort_values("absdiff", ascending=False).head(8)
        print("  largest disagreements:")
        for _, t in top.iterrows():
            print(f"    {t['id']} {t['family']:7s} n_rot={int(t['n_rot'])}  "
                  f"ours={t[col_our]:.3f} d3={t[col_d3]:.3f} |d|={t['absdiff']:.3f}")
        return s

    e = block("our_electron", "d3_electron", "ELECTRON reorg (neutral<->anion)")
    h = block("our_hole", "d3_hole", "HOLE reorg (neutral<->cation)")

    out = BASE / "comparison.csv"
    ok.to_csv(out, index=False)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
