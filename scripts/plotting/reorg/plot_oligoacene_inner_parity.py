#!/usr/bin/env python
"""Parity figure: our inner-sphere hole lambda vs literature, 4-point Nelsen, LEVEL-MATCHED.

Each point is compared to its OWN same-basis literature reference, and the FULL level of theory
(functional + basis) is printed under the molecule name (and color-coded), because the anchors
come from two sources at two bases:
  - naphthalene/anthracene/tetracene/pentacene : UB3LYP/6-311G**   vs Deng & Goddard 2004 (Table 2)
  - tetracene/pentacene/rubrene                : UB3LYP/6-31G(d,p)  vs da Silva Filho, Kim & Bredas 2005
(tetracene & pentacene appear at BOTH bases -> shows the value barely shifts with basis.)

Independent validation of the reorg pipeline (different functional/basis/family than D3TaLES,
bound cations only -> no anion pathology). PNG only (no PDF unless asked).

  python scripts/plotting/reorg/plot_oligoacene_inner_parity.py
"""
from __future__ import annotations
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
import sys
sys.path.insert(0, str(ROOT / "scripts" / "plotting"))
import plot_style as ps  # noqa: E402

OUTDIR = ROOT / "results" / "figures" / "reorg"

# one block per (level, reference). `side` steers labels so doubled points don't collide.
BLOCKS = [
    dict(ref="Deng & Goddard 2004", functional="UB3LYP", basis="6-311G**",
         source="results/reorg_anchors/oligoacene_denggoddard_6311.csv",
         calc="results/reorg_anchors/inner_acene/b3lyp_6311gdp/calc",
         color="#2B6CB0", side="right",
         ids=["naphthalene", "anthracene", "tetracene", "pentacene"]),
    dict(ref="da Silva Filho 2005", functional="UB3LYP", basis="6-31G(d,p)",
         source="results/reorg_anchors/oligoacene_inner.csv",
         calc="results/reorg_anchors/inner_acene/calc",
         color="#DD6B20", side="left",
         ids=["tetracene", "pentacene", "rubrene"]),
]


def main():
    ps.apply_style()
    fig, ax = plt.subplots(figsize=(6.8, 6.8))
    lo, hi = 0.06, 0.20
    ax.plot([lo, hi], [lo, hi], color="#333333", lw=1.6, zorder=2, label="y = x")

    all_d = []
    texts, px, py = [], [], []
    for blk in BLOCKS:
        level = f"{blk['functional']}/{blk['basis']}"
        lit = {r["id"]: float(r["lit_hole_eV"]) for r in csv.DictReader(open(ROOT / blk["source"]))}
        calc = ROOT / blk["calc"]
        xs, ys, labs = [], [], []
        for gid in blk["ids"]:
            j = calc / f"{gid}.json"
            if not j.exists() or gid not in lit:
                continue
            our = json.loads(j.read_text()).get("our_hole")
            if our is None:
                continue
            xs.append(lit[gid]); ys.append(float(our)); labs.append(gid)
            all_d.append(abs(float(our) - lit[gid]))
        ax.scatter(xs, ys, s=140, c=blk["color"], edgecolors="white", linewidths=1.2, zorder=4,
                   label=f"{level}  ({blk['ref']})")
        # one 2-line label (name + level) per point; adjust_text de-overlaps them below
        for x, y, lab in zip(xs, ys, labs):
            texts.append(ax.text(x, y, f"{lab}\n{level}", fontsize=10, color=blk["color"],
                                 ha="center", va="center", zorder=6, linespacing=1.2))
            px.append(x); py.append(y)

    from adjustText import adjust_text
    adjust_text(texts, x=px, y=py, ax=ax, expand=(1.3, 1.6), force_text=(0.5, 0.8),
                arrowprops=dict(arrowstyle="-", color="#8A8A8A", lw=0.6))

    d = np.array(all_d); mad, mx = d.mean(), d.max()
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
    ax.set_xlabel("literature $\\lambda_{hole}$  (eV)")
    ax.set_ylabel("this work  $\\lambda_{hole}$  (eV)")
    ax.set_title("Oligoacene inner-sphere hole $\\lambda$\ngas, 4-point (level-matched per point)")
    ax.grid(True, color=ps.C["grid"], lw=0.6, alpha=0.30); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.text(0.04, 0.96, f"n = {len(d)}\nMAD = {mad*1000:.1f} meV\nmax = {mx*1000:.1f} meV",
            transform=ax.transAxes, va="top", ha="left", fontsize=14,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#CBD5E0"))
    ax.legend(loc="lower right", frameon=False, fontsize=11)

    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTDIR / "oligoacene_inner_parity.png", dpi=300, bbox_inches="tight")   # PNG only
    print(f"wrote {OUTDIR/'oligoacene_inner_parity.png'}  (n={len(d)}, MAD={mad*1000:.1f} meV, max={mx*1000:.1f} meV)")


if __name__ == "__main__":
    main()
