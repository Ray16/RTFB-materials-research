#!/usr/bin/env python
"""lambda vs 2e-capacity Pareto for the RFB grafting decision (the two payoff axes:
fast kinetics = low lambda_i; high gravimetric capacity). New merrifield candidates + the
first-batch grafted candidates for context. No scalarized figure-of-merit (raw axes, marked
non-dominated front). Confidence encoded: filled = QC-clean lambda, open = preliminary
(conformer-flagged, being refined). Data read live from results/reorganization.csv +
results/scorecard.csv so it stays in sync.
"""
from __future__ import annotations
import csv
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_style import apply_style  # noqa
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
R = ROOT / "results"
FIG = R / "figures" / "candidates"; FIG.mkdir(parents=True, exist_ok=True)

FAMCOL = {"quinone": "#0072B2", "viologen": "#CC79A7", "imide": "#009E73", "phenol": "#762A83"}

# --- merrifield: primary-couple lambda (+clean/flag) from reorganization.csv; cap/SA verified
PRIM = {"aq_benzyloxy": "neu->red1", "aq_benzylamino": "neu->red1", "nq_benzyloxy": "neu->red1",
        "bisviologen": "ox2->ox1", "dtbc_phenol": "ox->neu"}
MERRI = {  # id: (display, family, capacity mAh/g, n_e)
    "aq_benzyloxy":  ("benzyloxy-AQ", "quinone", 163.2, 2),
    "aq_benzylamino":("benzylamino-AQ", "quinone", 163.7, 2),
    "nq_benzyloxy":  ("benzyloxy-NQ", "quinone", 192.6, 2),
    "bisviologen":   ("bis-viologen", "viologen", 146.3, 2),
    "dtbc_phenol":   ("DTB-catechol (1e)", "phenol", 82.1, 1),
}
lam_flag = {}
for r in csv.DictReader(open(R / "reorganization.csv")):
    if r["id"] in PRIM and r["couple"] == PRIM[r["id"]]:
        try: lam_flag[r["id"]] = (float(r["lambda_i_eV"]), bool(r.get("flag")))
        except: pass

pts = []  # (lam, cap, color, name, cohort, clean)
for i, (nm, fam, cap, ne) in MERRI.items():
    if i in lam_flag:
        lam, flag = lam_flag[i]
        pts.append((lam, cap, FAMCOL[fam], nm, "merrifield", not flag))

# --- batch-1 grafted candidates from scorecard (QC-clean lambda already filtered)
FAM1 = {"pyridine-multi-e": "viologen", "imide (n-type)": "imide", "quinone (n-type)": "quinone"}
N1 = {"viologen": "Me-viologen", "pmdi": "PMDI", "ndi_ammonium": "ammonium-NDI",
      "mophquinone": "MeO-Ph-quinone", "dmophquinone": "(MeO)2-Ph-quinone"}
for r in csv.DictReader(open(R / "scorecard.csv")):
    if r.get("status") != "candidate":
        continue
    try:
        lam = float(r["lambda_i_eV"]); cap = float(r["specific_capacity_mAh_g"])
    except: continue
    fam = FAM1.get(r.get("family", ""), "quinone")
    pts.append((lam, cap, FAMCOL[fam], N1.get(r["id"], r["id"]), "batch1", True))

# --- non-dominated front (minimize lambda, maximize capacity) over ALL placed points
def nondominated(P):
    keep = []
    for a in P:
        dom = any((b[0] <= a[0] and b[1] >= a[1]) and (b[0] < a[0] or b[1] > a[1]) for b in P)
        if not dom: keep.append(a)
    return sorted(keep, key=lambda p: p[0])
front = nondominated(pts)

apply_style()
fig, ax = plt.subplots(figsize=(13, 9))
# desirable corner (low lambda, high capacity)
ax.axvspan(0.30, 0.45, color="#1B7837", alpha=0.05, zorder=0)
# front line
fx = [p[0] for p in front]; fy = [p[1] for p in front]
ax.plot(fx, fy, "--", color="#555", lw=1.8, zorder=1, label="non-dominated front")
texts = []
for lam, cap, col, nm, cohort, clean in pts:
    mk = "o" if cohort == "merrifield" else "s"
    if clean:
        ax.scatter([lam], [cap], s=230, marker=mk, color=col, edgecolor="white", linewidth=1.5, zorder=3)
    else:
        ax.scatter([lam], [cap], s=230, marker=mk, facecolor="none", edgecolor=col, linewidth=2.6, zorder=3)
    texts.append(ax.text(lam, cap, nm, fontsize=15, color="#222", zorder=4))
ax.set_xlabel("inner-sphere $\\lambda_i$ (eV)  —  lower = faster kinetics")
ax.set_ylabel("specific capacity (mAh g$^{-1}$)  —  higher = better")
ax.set_title("RFB payoff: reorganization energy vs. 2e capacity", fontweight="bold")
ax.set_xlim(0.30, 0.95); ax.set_ylim(60, 216)
from adjustText import adjust_text
adjust_text(texts, x=[p[0] for p in pts], y=[p[1] for p in pts], ax=ax,
            expand=(1.4, 1.8), force_text=(0.4, 0.6), force_static=(0.6, 0.9),
            arrowprops=dict(arrowstyle="-", color="#888", lw=0.9))
# legend: cohort + confidence + families
import matplotlib.lines as ml
leg1 = [ml.Line2D([], [], marker="o", ls="", color="#444", ms=12, label="merrifield (new)"),
        ml.Line2D([], [], marker="s", ls="", color="#444", ms=11, label="first batch"),
        ml.Line2D([], [], marker="o", ls="", mfc="none", mec="#444", mew=2.4, ms=12, label="open = preliminary $\\lambda$ (conformer)"),
        ml.Line2D([], [], ls="--", color="#555", label="non-dominated front")]
fam_h = [ml.Line2D([], [], marker="o", ls="", color=FAMCOL[f], ms=11, label=f) for f in ["quinone","viologen","imide","phenol"]]
l1 = ax.legend(handles=leg1, loc="lower center", fontsize=14, title="cohort / confidence", title_fontsize=15)
ax.add_artist(l1)
ax.legend(handles=fam_h, loc="upper right", fontsize=14, title="redox family", title_fontsize=15)
fig.text(0.5, -0.02,
         "Absolute $\\lambda_i$ systematic uncertainty ~±0.2 eV > point spread, so ranking is relative.\n"
         "Open points: grafted $\\lambda$ inflated by tether rotation (conformer-matching in progress). 1e phenol shown for context.",
         ha="center", va="top", fontsize=14, style="italic", color="#555")
fig.tight_layout()
out = FIG / "pareto_lambda_capacity.png"
fig.savefig(out, bbox_inches="tight")
print("wrote", out, "| front:", [p[3] for p in front])
