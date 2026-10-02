#!/usr/bin/env python
"""Parity figure for the sourced MeCN benchmark (results/validation/benchmark_mecn.csv,
written by scripts/validation/benchmark/score_benchmark_mecn.py): computed vs experimental
E (V vs Fc/Fc+), tier A filled / tier B open, coloured by family; n, MAE and bias per tier in
the legend. Excluded waves are not plotted. PNG only.

  python scripts/plotting/validation/plot_benchmark_mecn.py
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_style import C, FAMILY_COLOR, apply_style  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np               # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "results" / "validation" / "benchmark_mecn.csv"
OUT = ROOT / "results" / "figures" / "validation" / "benchmark_mecn_parity.png"
FAM = dict(FAMILY_COLOR, **{"imide (n-type)": "#009E73"})


def main():
    rows = [r for r in csv.DictReader(open(SRC))
            if r["quantity"].startswith("E[") and r["err_V"] and not r["excluded_molecule"]
            and r["tier"] in ("A", "B")]
    apply_style()
    fig, ax = plt.subplots(figsize=(9.5, 9.0))
    lo, hi = -2.3, 0.6
    ax.fill_between([lo, hi], [lo - 0.2, hi - 0.2], [lo + 0.2, hi + 0.2], color="#B8C0C8",
                    alpha=0.28, lw=0, zorder=1)
    ax.plot([lo, hi], [lo, hi], color="#333333", lw=1.6, zorder=2)
    for tier, filled in (("A", True), ("B", False)):
        sel = [r for r in rows if r["tier"] == tier]
        if not sel:
            continue
        x = np.array([float(r["exp_V"]) for r in sel]); y = np.array([float(r["calc_V"]) for r in sel])
        cols = [FAM.get(r["family"], C["measured"]) for r in sel]
        d = y - x
        ax.scatter(x, y, s=95, zorder=3, facecolors=cols if filled else "none",
                   edgecolors=cols, linewidths=2.0,
                   label=f"tier {tier}: n={len(sel)}, MAE {np.mean(np.abs(d)):.2f} V, "
                         f"bias {np.mean(d):+.2f} V")
    for fam in sorted({r["family"] for r in rows if r["family"]}):
        ax.scatter([], [], s=95, color=FAM.get(fam, C["measured"]), label=fam)
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
    ax.set_xlabel(r"experimental $E$ (V vs Fc/Fc$^+$, MeCN)")
    ax.set_ylabel(r"computed $E$ (V vs Fc/Fc$^+$)")
    ax.grid(True, color=C["grid"], lw=0.6, alpha=0.3)
    ax.legend(loc="upper left", frameon=True, framealpha=0.95)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, bbox_inches="tight"); plt.close(fig)
    print(f"  -> {OUT}  ({len(rows)} points)")


if __name__ == "__main__":
    main()
