#!/usr/bin/env python
"""Total reorganization energy = inner-sphere (Nelsen 4-point, gas) + outer-sphere (Born
continuum, MeCN) for the merrifield candidates. Shows WHERE the outer-sphere term sits (it
dominates in MeCN). Inner lambda read from reorganization.csv (primary couple); outer (Born)
from results/lambda_outer_merrifield_born.json. Inner clean = solid, conformer-flagged =
hatched (preliminary). Outer is a Born estimate (PCM ~1.2-1.3x higher)."""
from __future__ import annotations
import csv, json
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_style import apply_style  # noqa
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]; R = ROOT / "results"
FIG = R / "figures" / "candidates"; FIG.mkdir(parents=True, exist_ok=True)

PRIM = {"aq_benzyloxy":"neu->red1","aq_benzylamino":"neu->red1","nq_benzyloxy":"neu->red1",
        "bisviologen":"ox2->ox1","dtbc_phenol":"ox->neu"}
DISP = {"aq_benzyloxy":"benzyloxy-AQ","aq_benzylamino":"benzylamino-AQ","nq_benzyloxy":"benzyloxy-NQ",
        "bisviologen":"bis-viologen","dtbc_phenol":"DTB-catechol (1e)"}
FAM = {"aq_benzyloxy":"#0072B2","aq_benzylamino":"#0072B2","nq_benzyloxy":"#0072B2",
       "bisviologen":"#CC79A7","dtbc_phenol":"#762A83"}

li = {}
for r in csv.DictReader(open(R/"reorganization.csv")):
    if r["id"] in PRIM and r["couple"]==PRIM[r["id"]]:
        try: li[r["id"]]=(float(r["lambda_i_eV"]), bool(r.get("flag")))
        except: pass
lo = json.load(open(R/"lambda_outer_merrifield_born.json"))

ids = sorted(PRIM, key=lambda i: li[i][0]+lo[i]["lam_o_born_eV"])
apply_style()
fig, ax = plt.subplots(figsize=(13, 7.5))
y = np.arange(len(ids))
for k,i in enumerate(ids):
    lam_i, flag = li[i]; lam_o = lo[i]["lam_o_born_eV"]; col = FAM[i]
    ax.barh(k, lam_i, color=col, edgecolor="white",
            hatch=("//" if flag else None), alpha=0.95, height=0.6, zorder=3)
    ax.barh(k, lam_o, left=lam_i, color=col, alpha=0.32, edgecolor="white", height=0.6, zorder=3)
    ax.text(lam_i/2, k, f"{lam_i:.2f}", va="center", ha="center", fontsize=15, color="white", fontweight="bold")
    ax.text(lam_i+lam_o/2, k, f"{lam_o:.2f}", va="center", ha="center", fontsize=15, color="#333")
    ax.text(lam_i+lam_o+0.02, k, f"$\\lambda_{{tot}}$ = {lam_i+lam_o:.2f}", va="center", fontsize=15, fontweight="bold", color="#222")
ax.set_yticks(y); ax.set_yticklabels([DISP[i] for i in ids], fontsize=15)
ax.set_xlabel("reorganization energy $\\lambda$ (eV)")
ax.set_xlim(0, 1.85)
ax.set_title("Total $\\lambda$ = inner-sphere + outer-sphere (MeCN)", fontweight="bold")
import matplotlib.patches as mp
h=[mp.Patch(facecolor="#555", label="inner-sphere $\\lambda_i$ (clean)"),
   mp.Patch(facecolor="#555", hatch="//", label="inner-sphere $\\lambda_i$ (preliminary, conformer)"),
   mp.Patch(facecolor="#555", alpha=0.32, label="outer-sphere $\\lambda_o$ (Born, MeCN)")]
ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.45, -0.16), ncol=2, fontsize=15, frameon=False)
fig.text(0.5,-0.05,
   "Outer-sphere = single-sphere Born (SASA radius); PCM continuum runs ~1.2-1.3x higher.\n"
   "Outer-sphere dominates in MeCN and is smoother across candidates, so it shifts $\\lambda_{tot}$ up\n"
   "roughly uniformly and preserves the inner-sphere ranking.",
   ha="center", fontsize=15, style="italic", color="#555")
fig.tight_layout()
out=FIG/"lambda_decomposition.png"; fig.savefig(out, bbox_inches="tight")
print("wrote", out)
