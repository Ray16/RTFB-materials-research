#!/usr/bin/env python
"""Publication-ready scatter plots comparing our (matched-protocol, B3LYP/6-31G* + RI-J)
reorganization energies against D3TaLES, for the quinone+imide validation set.

Fig 1 (parity): x=D3TaLES lambda, y=our lambda; y=x is perfect agreement. Panels: electron, hole.
Fig 2 (landscape): x=electron reorg, y=hole reorg; our points vs D3TaLES points.

Styling via src/redox/plotstyle.py (18 pt text, no overlaps, top-journal theme).
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from redox.plotstyle import set_pub_style, FAM_COLOR, OKABE_ITO
set_pub_style()

BASE = Path(__file__).resolve().parents[1] / "results" / "d3tales_reorg_validation"
OURS_C, D3_C = OKABE_ITO["blue"], OKABE_ITO["vermillion"]
WIN = 2.0  # eV display window; larger artifacts noted off-scale

df = pd.read_csv(BASE / "comparison.csv")
for c in ["our_electron", "d3_electron", "our_hole", "d3_hole"]:
    df[c] = pd.to_numeric(df[c], errors="coerce")


# ---------- Fig 1: parity (ours vs D3TaLES) ----------
fig, axes = plt.subplots(1, 2, figsize=(15, 7.6))
handles = None
for ax, (co, cd, title) in zip(axes, [("our_electron", "d3_electron", "Electron reorg"),
                                      ("our_hole", "d3_hole", "Hole reorg")]):
    s = df.dropna(subset=[co, cd])
    off = int(((s[co] > WIN) | (s[cd] > WIN)).sum())
    ax.axhspan(1.5, WIN, color="0.92", zorder=0)
    ax.axvspan(1.5, WIN, color="0.92", zorder=0)
    ax.plot([0, WIN], [0, WIN], "--", color="0.35", lw=1.6, zorder=1, label="perfect agreement")
    for fam, g in s.groupby("family"):
        ax.scatter(g[cd].clip(upper=WIN), g[co].clip(upper=WIN), s=42, alpha=0.75,
                   color=FAM_COLOR.get(fam, "0.5"), edgecolor="white", linewidth=0.4,
                   label=fam, zorder=3)
    ax.set_xlim(0, WIN); ax.set_ylim(0, WIN); ax.set_aspect("equal")
    ax.set_xticks([0, 0.5, 1.0, 1.5, 2.0]); ax.set_yticks([0, 0.5, 1.0, 1.5, 2.0])
    ax.set_xlabel("D3TaLES  $\\lambda$  (eV)")
    ax.set_ylabel("our  $\\lambda$  (eV)")
    ax.set_title(title, pad=10)
    ad = (s[co] - s[cd]).abs()
    # stats box in the empty upper-left triangle
    ax.text(0.04, 0.96, f"n = {len(s)}\nMAD = {ad.mean():.2f} eV\n{off} off-scale "
            f"($\\rightarrow${s[cd].max():.1f})", transform=ax.transAxes, va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="0.7", lw=1))
    if handles is None:
        handles, labels = ax.get_legend_handles_labels()
# 'outside' placement reserves its own band so it never overlaps the axis labels
fig.legend(handles, labels, loc="outside lower center", ncol=len(labels))
fig.suptitle("Reorganization energy: ours vs D3TaLES  (matched protocol, B3LYP/6-31G*)")
f1 = BASE / "reorg_parity_ours_vs_d3tales.png"
fig.savefig(f1); plt.close(fig)


# ---------- Fig 2: electron-vs-hole landscape ----------
fig, ax = plt.subplots(figsize=(8.8, 8.4))
o = df.dropna(subset=["our_electron", "our_hole"])
d = df.dropna(subset=["d3_electron", "d3_hole"])
ax.axhspan(1.5, WIN, color="0.94", zorder=0); ax.axvspan(1.5, WIN, color="0.94", zorder=0)
ax.scatter(d["d3_electron"].clip(upper=WIN), d["d3_hole"].clip(upper=WIN), s=48, alpha=0.6,
           color=D3_C, marker="s", edgecolor="white", linewidth=0.4, label=f"D3TaLES (n={len(d)})")
ax.scatter(o["our_electron"].clip(upper=WIN), o["our_hole"].clip(upper=WIN), s=48, alpha=0.75,
           color=OURS_C, marker="o", edgecolor="white", linewidth=0.4, label=f"ours (n={len(o)})")
ax.set_xlim(0, WIN); ax.set_ylim(0, WIN); ax.set_aspect("equal")
ax.set_xticks([0, 0.5, 1.0, 1.5, 2.0]); ax.set_yticks([0, 0.5, 1.0, 1.5, 2.0])
ax.set_xlabel("electron reorg  $\\lambda_e$  (eV)")
ax.set_ylabel("hole reorg  $\\lambda_h$  (eV)")
ax.set_title("Reorganization-energy landscape")
ax.legend(loc="upper right")
ax.text(0.04, 0.96, "shaded: implausible\n(> 1.5 eV)", transform=ax.transAxes, va="top",
        ha="left", color="0.4")
f2 = BASE / "reorg_landscape_e_vs_h.png"
fig.savefig(f2); plt.close(fig)

print(f"wrote {f1}")
print(f"wrote {f2}")
