#!/usr/bin/env python
"""Parity figure: our inner-sphere reorganization energy vs literature for the perfluoroacenes,
BOTH couples (hole + electron), level-matched to Ruiz Delgado et al. JACS 2009 Table 5 (AP)
at B3LYP/6-31G**. Companion to the oligoacene hole parity figure. PNG only.

  python scripts/plotting/reorg/plot_perfluoroacene_parity.py
"""
from __future__ import annotations
import csv, json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
import sys; sys.path.insert(0, str(ROOT / "scripts" / "plotting"))
import plot_style as ps  # noqa: E402

OUT = ROOT / "results" / "figures" / "reorg" / "perfluoroacene_parity.png"
SRC = ROOT / "results" / "reorg_anchors" / "electron_acceptors.csv"
# NOTE: generic worker writes every tag under results/reorg_anchors/electron/<tag>/calc
CALC = {"electron": ROOT / "results/reorg_anchors/electron/b3lyp_631gdp_electron/calc",
        "hole":     ROOT / "results/reorg_anchors/electron/b3lyp_631gdp_hole/calc"}
COLOR = {"electron": "#2B6CB0", "hole": "#DD6B20"}
LABEL = {"perfluorotetracene": "PFT", "perfluoropentacene": "PFP"}


def main():
    ps.apply_style()
    lit = {r["id"]: r for r in csv.DictReader(open(SRC))}
    fig, ax = plt.subplots(figsize=(6.6, 6.6))
    lo, hi = 0.20, 0.31
    ax.plot([lo, hi], [lo, hi], color="#333333", lw=1.6, zorder=2, label="y = x")
    LEVEL = "UB3LYP/6-31G(d,p)"
    alld = []
    texts, px, py = [], [], []
    for couple in ("electron", "hole"):
        xs, ys, labs = [], [], []
        for gid, r in lit.items():
            j = CALC[couple] / f"{gid}.json"
            v = r.get(f"lit_{couple}_eV")
            if not j.exists() or not v:
                continue
            our = json.loads(j.read_text()).get("our_lambda")
            if our is None:
                continue
            xs.append(float(v)); ys.append(float(our)); labs.append(LABEL.get(gid, gid))
            alld.append(abs(float(our) - float(v)))
        ax.scatter(xs, ys, s=150, c=COLOR[couple], edgecolors="white", linewidths=1.2,
                   zorder=4, label=f"{couple} λ")
        # one 2-line label (molecule+couple, then level) per point; adjust_text de-overlaps
        for x, y, lab in zip(xs, ys, labs):
            texts.append(ax.text(x, y, f"{lab} {couple}\n{LEVEL}", fontsize=10.5,
                                 color=COLOR[couple], ha="center", va="center",
                                 zorder=6, linespacing=1.2))
            px.append(x); py.append(y)

    from adjustText import adjust_text
    adjust_text(texts, x=px, y=py, ax=ax, expand=(1.3, 1.6), force_text=(0.5, 0.9),
                arrowprops=dict(arrowstyle="-", color="#8A8A8A", lw=0.6))
    d = np.array(alld); mad = d.mean() if len(d) else float("nan")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
    ax.set_xlabel("literature $\\lambda$  (eV)")
    ax.set_ylabel("this work  $\\lambda$  (eV)")
    ax.set_title("Perfluoroacene inner-sphere $\\lambda$\nB3LYP/6-31G**, gas, 4-point  (vs Delgado 2009, Table 5)")
    ax.grid(True, color=ps.C["grid"], lw=0.6, alpha=0.30); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.text(0.04, 0.96, f"n = {len(d)}\nMAD = {mad*1000:.1f} meV",
            transform=ax.transAxes, va="top", ha="left", fontsize=14,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#CBD5E0"))
    ax.legend(loc="lower right", frameon=False, fontsize=12)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=300, bbox_inches="tight")   # PNG only
    print(f"wrote {OUT}  (n={len(d)}, MAD={mad*1000:.1f} meV)")


if __name__ == "__main__":
    main()
