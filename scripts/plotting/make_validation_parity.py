#!/usr/bin/env python
"""Fig 1b — validation parity: computed vs measured redox potential (vs Fc), 44 tier-A couples.
Reads archive benchmark_mecn.csv. Colored by family; y=x identity + stratified-bias guide;
MAE / bias annotated. Publication style (18 pt, no overlaps, PNG only)."""
from __future__ import annotations
import csv, sys
sys.path.insert(0, "/nfs/lambda_stor_01/homes/rzhu/0_redox_refactor/scripts/plotting")
import numpy as np, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from plot_style import apply_style, C

BM = "/nfs/lambda_stor_01/homes/rzhu/archive/0_redox/results/validation/benchmark_mecn.csv"
OUT = "/nfs/lambda_stor_01/homes/rzhu/0_redox_refactor/results/figures/validation/parity_mecn.png"

def num(x):
    try: return float(x)
    except: return None

rows = [r for r in csv.DictReader(open(BM))
        if r.get("tier","").strip().upper().startswith("A")
        and num(r["exp_V"]) is not None and num(r["calc_V"]) is not None]

FAM = {"quinone": "#0072B2", "imide": "#009E73", "pyridine": "#CC79A7",
       "nitroxide": "#E69F00", "other": "#999999"}
def fk(f):
    f = (f or "").lower()
    for k in ("quinone","imide","pyridine","nitroxide"):
        if k in f: return k
    return "other"

xs = np.array([num(r["exp_V"]) for r in rows])
ys = np.array([num(r["calc_V"]) for r in rows])
cols = [FAM[fk(r["family"])] for r in rows]
d = ys - xs
mae = float(np.mean(np.abs(d))); bias = float(np.mean(d))

apply_style()
fig, ax = plt.subplots(figsize=(8.2, 8.0))
lo, hi = min(xs.min(), ys.min()) - 0.2, max(xs.max(), ys.max()) + 0.2
ax.fill_between([lo, hi], [lo-0.226, hi-0.226], [lo+0.226, hi+0.226],
                color="#B8C0C8", alpha=0.25, lw=0, zorder=1, label="±MAE")
ax.plot([lo, hi], [lo, hi], color="#333333", lw=1.6, zorder=2)
ax.plot([lo, hi], [lo+bias, hi+bias], color=C["warn"], lw=1.4, ls="--", zorder=2)
ax.scatter(xs, ys, c=cols, s=70, edgecolors="black", linewidths=0.6, alpha=0.9, zorder=3)
ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
ax.set_xlabel("Measured $E$ vs Fc/Fc$^+$ (V)")
ax.set_ylabel("Computed $E$ vs Fc/Fc$^+$ (V)")
ax.set_title("Potential validation — 44 tier-A couples")
ax.text(0.04, 0.96, f"MAE {mae:.3f} V\nbias {bias:+.3f} V\nn = {len(rows)}",
        transform=ax.transAxes, va="top", ha="left", fontsize=16,
        bbox=dict(boxstyle="round,pad=0.4", fc="#F7F7F7", ec="#CCCCCC"))
handles = [Line2D([0],[0], marker="o", ls="", mfc=FAM[k], mec="black", ms=11,
                  label=lab) for k, lab in
           [("quinone","quinone"),("imide","imide"),("pyridine","pyridine/viologen"),
            ("nitroxide","nitroxide")]]
handles.append(Line2D([0],[0], color=C["warn"], ls="--", lw=1.4, label="systematic bias"))
ax.legend(handles=handles, loc="lower right", fontsize=13, framealpha=0.95)
fig.savefig(OUT, bbox_inches="tight")
print("wrote", OUT, f"| MAE {mae:.3f} bias {bias:+.3f} n {len(rows)}")
