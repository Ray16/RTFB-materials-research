#!/usr/bin/env python
"""Fig 3 — dihydrophenazine (DHP) catholyte design rules + per-couple validation.
One panel: dumbbell per DHP variant (E1 -> E2 couples vs Fc), grouped by design lever, with the
measured + computed BMEPZ parent as the validation reference. All numbers from the scored pool +
the BMEPZ anchor run. Publication style (18 pt, no overlaps, PNG only)."""
from __future__ import annotations
import csv, sys
sys.path.insert(0, "/nfs/lambda_stor_01/homes/rzhu/0_redox_refactor/scripts/plotting")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from plot_style import apply_style, C

CSV = "/nfs/lambda_stor_01/homes/rzhu/0_redox_refactor/results/validation/full_pool_scored_legacyG_MeCN_TBAPF6.csv"
OUT = "/nfs/lambda_stor_01/homes/rzhu/0_redox_refactor/results/figures/candidates/dhp_design_rules.png"
rows = {r["id"]: r for r in csv.DictReader(open(CSV))}
def couples(gid): return [float(x) for x in rows[gid]["E_couples"].split(";")]

# (label, id-or-None, E1, E2, group)  group: ncap / ring / ref
items = [
    ("BMEPZ parent (measured)",        None, -0.29, 0.50, "ref"),
    ("BMEPZ parent (computed)",        None, -0.062, 0.869, "ref"),
    ("N'-methyl",      "exp_dmpz_mb_methyl",          None, None, "ncap"),
    ("N'-ethyl",       "exp_dmpz_mb_ethyl",           None, None, "ncap"),
    ("N'-(2-methoxyethyl)", "hyp_dmpz_mb_methoxyethyl", None, None, "ncap"),
    ("N'-(2,2,2-trifluoroethyl) [EWG]", "hyp_dmpz_mb_trifluoroethyl", None, None, "ncap"),
    ("2,7-dimethyl [ring donor]",  "hyp_dmpz_27dimethyl_ethyl",  None, None, "ring"),
    ("2,7-dimethoxy [ring donor]", "hyp_dmpz_27dimethoxy_ethyl", None, None, "ring"),
]
data = []
for lab, gid, e1, e2, grp in items:
    if gid is not None:
        c = couples(gid); e1, e2 = c[0], c[1]
    data.append((lab, e1, e2, grp))
# order: refs on top, then grafted by mean potential
refs = [d for d in data if d[3] == "ref"]
graf = sorted([d for d in data if d[3] != "ref"], key=lambda d: (d[1]+d[2])/2)
ordered = refs + graf

GC = {"ncap": "#E69F00", "ring": "#0072B2", "ref": "#4D4D4D"}
apply_style()
fig, ax = plt.subplots(figsize=(11.5, 7.2))
ys = list(range(len(ordered)))[::-1]
for y, (lab, e1, e2, grp) in zip(ys, ordered):
    style = dict(color=GC[grp], lw=3 if grp != "ref" else 2.4, alpha=0.9,
                 ls="-" if grp != "ref" else (0, (4, 2)))
    ax.plot([e1, e2], [y, y], **style, zorder=2)
    ax.scatter([e1, e2], [y, y], c=GC[grp], s=120, edgecolors="black",
               linewidths=1.0, zorder=3)
    ax.annotate(f"{e1:+.2f}", (e1, y), textcoords="offset points", xytext=(-6, 9),
                ha="right", fontsize=12, color="#333333")
    ax.annotate(f"{e2:+.2f}", (e2, y), textcoords="offset points", xytext=(6, 9),
                ha="left", fontsize=12, color="#333333")
ax.set_yticks(ys); ax.set_yticklabels([d[0] for d in ordered], fontsize=14)
ax.axvline(0, color="#999999", lw=1.0, ls=":", zorder=1)
ax.set_xlabel("Oxidation potential couples $E_1,\\,E_2$ vs Fc/Fc$^+$ (V)")
ax.set_title("Dihydrophenazine catholyte design rules + BMEPZ validation")
ax.set_xlim(-0.55, 1.25)
# group legend — placed OUTSIDE the axes (upper right) so it never overlaps a dumbbell
handles = [Line2D([0],[0], color=GC["ncap"], lw=3, marker="o", mec="black", ms=10, label="N'-cap series"),
           Line2D([0],[0], color=GC["ring"], lw=3, marker="o", mec="black", ms=10, label="ring-substituent series"),
           Line2D([0],[0], color=GC["ref"], lw=2.4, ls=(0,(4,2)), marker="o", mec="black", ms=10, label="BMEPZ parent (ref)")]
leg = ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0),
                fontsize=13, framealpha=0.95)
fig.savefig(OUT, bbox_inches="tight", bbox_extra_artists=[leg])
print("wrote", OUT)
