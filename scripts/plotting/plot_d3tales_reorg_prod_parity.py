#!/usr/bin/env python
"""Track-2 D3TaLES reorganization sweep: 3-panel parity showing WHY the gap vs D3TaLES is not a
solvent effect.

  (a) our SMD-opt   vs D3TaLES(reported)
  (b) our B3LYP-gas  vs D3TaLES(reported)     <- same gap as (a) with NO solvent anywhere
  (c) our SMD-opt   vs our B3LYP-gas          <- solvent-geometry shift: tiny

Reads results/d3tales_reorg_validation_prod/comparison.csv (from aggregate_d3tales_reorg_prod.py).
  PYTHONPATH=src python scripts/plotting/plot_d3tales_reorg_prod_parity.py
"""
from __future__ import annotations
import csv
import statistics
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_style import apply_style, FAM_COLOR  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
CSV = ROOT / "results" / "d3tales_reorg_validation_prod" / "comparison.csv"
OUT = ROOT / "results" / "figures" / "reorg" / "reorg_prod_parity_3panel.png"

LO, HI = 0.0, 2.0   # sane window (drop D3TaLES garbage / broken geoms)


def _f(x):
    try:
        v = float(x); return v if v == v else None
    except (TypeError, ValueError):
        return None


def _famcolor(fam):
    for k, c in FAM_COLOR.items():
        if k in (fam or ""):
            return c
    return "#888888"


def main():
    rows = list(csv.DictReader(CSV.open()))
    apply_style()
    panels = [
        ("our_d3level_eV", "d3tales_eV", "our LC-$\\omega$HPBE (their level)  vs  D3TaLES", "D3TaLES reported $\\lambda$ (eV)", "our LC-$\\omega$HPBE $\\lambda$ (eV)"),
        ("our_smd_eV", "d3tales_eV", "our production (SMD)  vs  D3TaLES", "D3TaLES reported $\\lambda$ (eV)", "our production $\\lambda$ (eV)"),
        ("our_smd_eV", "our_d3level_eV", "our production  vs  our LC-$\\omega$HPBE", "our LC-$\\omega$HPBE $\\lambda$ (eV)", "our production $\\lambda$ (eV)"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(19, 6.4))
    for ax, (xk, yk, title, xl, yl) in zip(axes, panels):
        # QC: drop molecules whose inner-sphere lambda is untrustworthy (unbound vertical anion /
        # unrescuable outlier — see finalize_reorg_qc.py). Then keep a sane window (drops D3TaLES's
        # own broken-geometry garbage on the reference axis).
        rel = [r for r in rows if str(r.get("reliable", "True")) == "True"]
        n_excl = len(rows) - len(rel)
        pts = [(_f(r[xk]), _f(r[yk]), _famcolor(r["family"])) for r in rel]
        pts = [(x, y, c) for x, y, c in pts if x is not None and y is not None and LO < x < HI and LO < y < HI]
        xs = [p[0] for p in pts]; ys = [p[1] for p in pts]; cs = [p[2] for p in pts]
        ax.scatter(xs, ys, c=cs, s=42, alpha=0.75, edgecolors="white", linewidths=0.4, zorder=3)
        lo, hi = 0.0, 1.4
        ax.plot([lo, hi], [lo, hi], "-", color="#444444", lw=1.4, zorder=2, label="y = x")
        diffs = [abs(a - b) for a, b in zip(xs, ys)]
        signed = [a - b for a, b in zip(xs, ys)]
        mad = statistics.mean(diffs); med = statistics.median(diffs); bias = statistics.mean(signed)
        ax.text(0.04, 0.96, f"n = {len(xs)}\nMAD = {mad:.3f} eV\nmedian = {med:.3f}\nbias = {bias:+.3f}"
                f"\nexcl = {n_excl} (unbound $\\lambda$)",
                transform=ax.transAxes, va="top", ha="left", fontsize=14,
                bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#bbbbbb", alpha=0.9))
        ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
        ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(title, fontsize=18, fontweight="bold")
        ax.grid(True, alpha=0.3)
    # family legend
    handles = [plt.Line2D([0], [0], marker="o", ls="", mfc=c, mec="white", ms=11, label=k)
               for k, c in FAM_COLOR.items()]
    axes[2].legend(handles=handles, loc="lower right", fontsize=13, title="family", framealpha=0.9)
    fig.suptitle("Inner-sphere $\\lambda$ (electron): matching D3TaLES's IP-tuned LC-$\\omega$HPBE reproduces their $\\lambda$; "
                 "our production level differs by functional+basis (not solvent)",
                 fontsize=17, fontweight="bold", y=1.02)
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  -> {OUT}")


if __name__ == "__main__":
    main()
