#!/usr/bin/env python
"""Full-pool candidate landscape: redox potential (vs Fc) x usable capacity (TBAPF6/MeCN).

One canonical panel over all 88 computed candidates on a uniform legacy-G (xTB-RRHO)
thermal basis. Anolytes (E<0) and catholytes (E>0) with their sigma-aware Pareto fronts
(raw axes: most-negative/positive potential x highest usable capacity; no scalarized FoM).
Front members get numbered markers + a real-name key (no coined ids, no overlaps).
1e controls are drawn hollow and excluded from the fronts.
"""
from __future__ import annotations
import csv, sys, os, re
sys.path.insert(0, "/nfs/lambda_stor_01/homes/rzhu/0_redox_refactor/scripts/plotting")
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from plot_style import apply_style, FAM_COLOR, C

CSV = sys.argv[1]
OUT = sys.argv[2]

# ---- real, compact labels. The 4-methylbenzyl graft handle is the UNIVERSAL Merrifield
#      tether (same on every candidate) -> stated in the caption, dropped from per-point labels. ----
_CORE = [("anthracene-9,10-dione", "AQ"), ("anthraquinone", "AQ"),
         ("1,4-naphthoquinone", "1,4-NQ"), ("1,2-naphthoquinone", "1,2-NQ"), ("naphthoquinone", "NQ"),
         ("1,4-benzoquinone", "1,4-BQ"), ("1,2-benzoquinone", "1,2-BQ"), ("benzoquinone", "BQ"),
         ("5,10-dihydrophenazine", "DHP"), ("dihydrophenazine", "DHP"), ("phenothiazine", "PTZ")]
_C = {"pmdi": "PMDI", "ndi": "NDI", "viologen": "viologen"}
_NGRP = {"ethyl": "N'-Et", "npropyl": "N'-nPr", "methyl": "N'-Me"}
_ARYL = {"H": "", "Me": "4-Me", "OMe": "4-OMe", "F": "4-F", "CN": "4-CN"}

def label_for(name: str) -> str:
    if " / " in name:  # coined design-series "core / tether / ngroup / aryl"
        p = [x.strip() for x in name.split("/")]
        bits = [_C.get(p[0], p[0])]
        if len(p) > 2 and p[2] in _NGRP: bits.append(_NGRP[p[2]])
        if len(p) > 3 and _ARYL.get(p[3]): bits.append(_ARYL[p[3]])
        if len(p) > 1 and p[1] == "ether_C2": bits.append("O-link")
        return " ".join(bits)
    s = name
    for lo, ab in _CORE:
        if lo in s.lower():
            i = s.lower().index(lo)
            pre = s[:i]
            if pre and pre[-1].isalnum(): pre += " "  # keep a space before the abbrev
            s = pre + ab + s[i + len(lo):]
            break
    s = re.sub(r"N-\(4-methylbenzyl\)-?", "", s)
    s = re.sub(r"\(?4-methylbenzyloxy\)?", "O-graft ", s)
    s = re.sub(r"benzyl\(ethyl\)amino", "N-ethylamino", s)
    s = re.sub(r"benzylamino", "amino", s)
    s = s.replace("(amino)", "amino")
    s = re.sub(r"\s+", " ", s).replace("graft -", "graft ").strip(" -,")
    return s

def fam_key(family: str) -> str:
    f = family.lower()
    if "quinone" in f: return "quinone"
    if "imide" in f: return "imide"
    if "pyridine" in f or "viologen" in f: return "viologen"
    if "phenazine" in f: return "phenazine"
    if "phenothiazine" in f: return "phenothiazine"
    return "other"

COL = dict(FAM_COLOR)
COL["phenazine"] = "#E69F00"      # orange  (2e catholyte)
COL["phenothiazine"] = "#762A83"  # purple  (1e/2e p-type)
COL["other"] = "#999999"

rows = list(csv.DictReader(open(CSV)))
for r in rows:
    r["E"] = float(r["E_mean"]); r["cap"] = float(r["cap_active_mAh_g"])
    r["ctrl"] = r["id"].endswith("_control") or r["n"] == "1"
    r["fk"] = fam_key(r["family"]); r["lab"] = label_for(r["name"])

def pareto(cands, want_low_E):
    front = []
    for a in cands:
        dom = False
        for b in cands:
            if b is a: continue
            cge = b["cap"] >= a["cap"]
            pge = (b["E"] <= a["E"]) if want_low_E else (b["E"] >= a["E"])
            strict = (b["cap"] > a["cap"]) or ((b["E"] < a["E"]) if want_low_E else (b["E"] > a["E"]))
            if cge and pge and strict: dom = True; break
        if not dom: front.append(a)
    return front

ano = [r for r in rows if r["pool"] == "anolyte" and not r["ctrl"]]
cat = [r for r in rows if r["pool"] == "catholyte" and not r["ctrl"]]
ctrls = [r for r in rows if r["ctrl"]]
af = sorted(pareto(ano, True), key=lambda r: r["E"])
cf = sorted(pareto(cat, False), key=lambda r: -r["E"])

apply_style()
fig, (ax, axk) = plt.subplots(1, 2, figsize=(16.5, 8.4),
                              gridspec_kw=dict(width_ratios=[2.45, 1.0], wspace=0.04))
axk.axis("off")

# background: all non-front candidates, faint, colored by family
front_ids = {r["id"] for r in af + cf}
for r in rows:
    if r["id"] in front_ids or r["ctrl"]: continue
    ax.scatter(r["E"], r["cap"], s=46, c=COL[r["fk"]], alpha=0.33, edgecolors="none", zorder=2)
# controls: hollow grey
for r in ctrls:
    ax.scatter(r["E"], r["cap"], s=60, facecolors="none", edgecolors="#777777",
               linewidths=1.4, zorder=3)

# Pareto fronts as connected step lines
def draw_front(front, color):
    xs = [r["E"] for r in front]; ys = [r["cap"] for r in front]
    ax.plot(xs, ys, "-", color=color, lw=2.2, alpha=0.85, zorder=4)

draw_front(af, "#2166AC"); draw_front(cf, "#B2182B")

# numbered markers for front members
key_lines = []
num = 0
for grp, front, mk in [("A", af, "o"), ("C", cf, "s")]:
    for r in front:
        num += 1
        ax.scatter(r["E"], r["cap"], s=150, c=COL[r["fk"]], edgecolors="black",
                   linewidths=1.3, marker=mk, zorder=6)
        ax.annotate(str(num), (r["E"], r["cap"]), ha="center", va="center",
                    fontsize=11, fontweight="bold", color="white", zorder=7)
        key_lines.append(f"{num}. {r['lab']}  ({r['E']:+.2f} V, {r['cap']:.0f})")

# axes
ax.axvline(0, color="#555555", lw=1.0, ls="--", alpha=0.7, zorder=1)
ax.set_xlabel("Redox potential  $E$  vs Fc/Fc$^+$  (V)")
ax.set_ylabel("Usable capacity  (mAh g$^{-1}$,  TBAPF$_6$ / MeCN)")
ax.set_title(f"Grafted redox-motif landscape — {len(rows)} candidates, uniform xTB-RRHO thermal basis")
ax.margins(y=0.08)  # headroom so the top front marker + labels never collide
ax.text(0.015, 0.97, "anolytes", transform=ax.transAxes, ha="left", va="top",
        fontsize=15, color="#2166AC", fontweight="bold")
ax.text(0.985, 0.55, "catholytes", transform=ax.transAxes, ha="right", va="top",
        fontsize=15, color="#B2182B", fontweight="bold")

# family legend (color) + marker legend
fam_handles = [Line2D([0],[0], marker="o", ls="", mfc=COL[k], mec="none", ms=11,
                      label=lab) for k, lab in
               [("quinone","quinone (AQ/NQ/BQ)"),("imide","imide (PMDI/NDI)"),
                ("viologen","viologen"),("phenazine","dihydrophenazine (DHP)"),
                ("phenothiazine","phenothiazine (PTZ)")]]
fam_handles += [Line2D([0],[0], marker="o", ls="", mfc="none", mec="#777777", ms=11,
                       label="1e control (excluded)")]
leg1 = ax.legend(handles=fam_handles, loc="lower center", ncol=2, fontsize=13,
                 title="family", title_fontsize=13, framealpha=0.95)
ax.add_artist(leg1)

# numbered key in its own right-hand panel (never clipped)
key_txt = ("Pareto front\n(○ anolyte,  □ catholyte)\n"
           "graft = 4-methylbenzyl (Merrifield proxy)\n\n" + "\n".join(key_lines))
axk.text(0.0, 1.0, key_txt, transform=axk.transAxes, ha="left", va="top", fontsize=12.5,
         family="monospace", linespacing=1.45,
         bbox=dict(boxstyle="round,pad=0.6", fc="#F7F7F7", ec="#CCCCCC"))

fig.savefig(OUT, bbox_inches="tight")
print("wrote", OUT)
print("\nANOLYTE FRONT:", len(af), " CATHOLYTE FRONT:", len(cf))
for i, l in enumerate(key_lines, 1): print(" ", l)
