#!/usr/bin/env python
"""Manuscript Figure 5 — sigma-aware Pareto fronts of the 91-candidate pool.

Two SEPARATE one-panel figures (JCIM house rule + STYLE_CONTRACT.md: one panel per file):

  fig5a_anolyte_pareto.png    85 n-type anolytes
  fig5b_catholyte_pareto.png   6 p-type catholytes

Axes are the two raw performance quantities, kept SEPARATE (no scalarized / weighted-sum
figure of merit):
    x = redox potential  E  vs Fc/Fc+   (V)
    y = usable capacity                  (mAh/g, TBAPF6/MeCN)

One point per candidate, coloured by redox family (colours per STYLE_CONTRACT §2, identical
across both panels). The Pareto-optimal front within each pool is marked (step line + outlined,
numbered markers keyed to real chemical names in a right-margin key).

"sigma-aware": front members carry a horizontal whisker of +/- 0.25 V (MAE 0.247), the
validation-benchmark potential uncertainty; front membership that would flip under a potential
shift smaller than that is not meaningful, shown without collapsing the two axes into a score.

STYLE_CONTRACT compliance: DMPZ (never DHP); PNG only (no pdf); every text element >= 18 pt;
zero overlapping elements (key + family legend live in a reserved right margin, clear of all
data and of each other); verbatim axis labels; frozen-91 pool enforced.

Data: results/validation/full_pool_scored_legacyG_MeCN_TBAPF6.csv (legacy xTB-RRHO tier).
The live expand_catholyte campaign may append extra catholyte rows; they are dropped (and
reported) so the figure matches the frozen 91-candidate manuscript pool.

CPU / matplotlib only. No GPU, no DFT, no number is altered.

Usage:
    python scripts/plotting/manuscript/fig5_pareto.py [input_csv] [output_dir]
"""
from __future__ import annotations
import csv
import re
import sys
from pathlib import Path

ROOT = Path("/nfs/lambda_stor_01/homes/rzhu/0_redox_refactor")
sys.path.insert(0, str(ROOT / "scripts" / "plotting"))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from plot_style import apply_style, BASE

CSV = Path(sys.argv[1]) if len(sys.argv) > 1 else \
    ROOT / "results/validation/full_pool_scored_legacyG_MeCN_TBAPF6.csv"
OUTDIR = Path(sys.argv[2]) if len(sys.argv) > 2 else ROOT / "manuscript/figures"

# Benchmark potential uncertainty (MAE over 39 tier-A redox endpoints, archive benchmark.json).
SIGMA_E = 0.247  # V
WHISKER_LABEL = "whisker ±0.25 V (MAE 0.247)"

# Live-campaign rows appended beyond the frozen 91-candidate pool (dropped, and reported).
# Frozen catholytes = 5 DMPZ phenazines + 1 phenothiazine (manuscript CONTEXT.md).
EXTRA_IDS = {"hyp_dmpz_27dimethoxy_ethyl"}

# Family -> colour (STYLE_CONTRACT §2; identical across both panels). Keyed by entity.
COL = {
    "quinone": "#0072B2",        # blue           (n-type anolyte)
    "imide": "#009E73",          # green          (n-type anolyte)
    "viologen": "#CC79A7",       # reddish-purple (pyridine-multi-e anolyte)
    "phenazine": "#E69F00",      # orange         (p-type catholyte, DMPZ)
    "phenothiazine": "#762A83",  # purple         (p-type catholyte, PTZ)
    "other": "#999999",
}
ANO_FRONT = "#2166AC"
CAT_FRONT = "#B2182B"

# Verbatim axis labels (STYLE_CONTRACT §3).
XLABEL = "Redox potential, $E$ vs Fc/Fc$^+$ (V)"
YLABEL = "Usable capacity (mAh g$^{-1}$, TBAPF$_6$/MeCN)"

# ---- compact real-chemistry labels (adapted from make_full_pool_figure.py) ----
_CORE = [("anthracene-9,10-dione", "AQ"), ("anthraquinone", "AQ"),
         ("1,4-naphthoquinone", "1,4-NQ"), ("1,2-naphthoquinone", "1,2-NQ"), ("naphthoquinone", "NQ"),
         ("1,4-benzoquinone", "1,4-BQ"), ("1,2-benzoquinone", "1,2-BQ"), ("benzoquinone", "BQ"),
         ("5,10-dihydrophenazine", "DMPZ"), ("dihydrophenazine", "DMPZ"), ("phenothiazine", "PTZ")]
_C = {"pmdi": "PMDI", "ndi": "NDI", "viologen": "viologen"}
_NGRP = {"ethyl": "N'-Et", "npropyl": "N'-nPr", "methyl": "N'-Me"}
_ARYL = {"H": "", "Me": "4-Me", "OMe": "4-OMe", "F": "4-F", "CN": "4-CN"}


def label_for(name: str) -> str:
    if " / " in name:  # coined design-series "core / tether / ngroup / aryl"
        p = [x.strip() for x in name.split("/")]
        bits = [_C.get(p[0], p[0])]
        if len(p) > 2 and p[2] in _NGRP:
            bits.append(_NGRP[p[2]])
        if len(p) > 3 and _ARYL.get(p[3]):
            bits.append(_ARYL[p[3]])
        if len(p) > 1 and p[1] == "ether_C2":
            bits.append("O-link")
        return " ".join(bits)
    s = name
    for lo, ab in _CORE:
        if lo in s.lower():
            i = s.lower().index(lo)
            pre = s[:i]
            if pre and pre[-1].isalnum():
                pre += " "
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
    if "quinone" in f:
        return "quinone"
    if "phenazine" in f:
        return "phenazine"
    if "phenothiazine" in f:
        return "phenothiazine"
    if "imide" in f:  # imide + mono-imide
        return "imide"
    if "pyridine" in f or "viologen" in f:
        return "viologen"
    return "other"


def _mm(v: float) -> str:
    """Format a potential with explicit + and a TRUE minus sign (STYLE_CONTRACT §3)."""
    return f"{v:+.2f}".replace("-", "−")


def pareto(cands, want_low_E):
    """Raw-axis Pareto front: maximise capacity AND (minimise E for anolytes /
    maximise E for catholytes). No scalarization."""
    front = []
    for a in cands:
        dominated = False
        for b in cands:
            if b is a:
                continue
            cge = b["cap"] >= a["cap"]
            pge = (b["E"] <= a["E"]) if want_low_E else (b["E"] >= a["E"])
            strict = (b["cap"] > a["cap"]) or \
                     ((b["E"] < a["E"]) if want_low_E else (b["E"] > a["E"]))
            if cge and pge and strict:
                dominated = True
                break
        if not dominated:
            front.append(a)
    return front


def load_frozen():
    rows = list(csv.DictReader(open(CSV)))
    dropped = [r for r in rows if r["id"] in EXTRA_IDS]
    rows = [r for r in rows if r["id"] not in EXTRA_IDS]
    for r in rows:
        r["E"] = float(r["E_mean"])
        r["cap"] = float(r["cap_mAh_g"])
        r["ctrl"] = r["id"].endswith("_control") or r["n"] == "1"
        r["fk"] = fam_key(r["family"])
        r["lab"] = label_for(r["name"])
    return rows, dropped


def _save(fig, stem):
    p = OUTDIR / f"{stem}.png"          # PNG ONLY (STYLE_CONTRACT §6)
    fig.savefig(p, dpi=300)
    plt.close(fig)
    return p


def _right_margin_blocks(fig, key_header, key_lines, fam_handles, fam_title):
    """Numbered real-name key (top) and family legend (bottom), both in the reserved
    right margin (axes end at 0.52) so neither can collide with data or each other."""
    x0 = 0.535
    txt = key_header + "\n" + "\n".join(key_lines)
    fig.text(x0, 0.90, txt, ha="left", va="top", fontsize=BASE, family="monospace",
             linespacing=1.35,
             bbox=dict(boxstyle="round,pad=0.5", fc="#F7F7F7", ec="#CCCCCC"))
    fig.legend(handles=fam_handles, loc="lower left", bbox_to_anchor=(x0, 0.12),
               fontsize=BASE, title=fam_title, title_fontsize=BASE, framealpha=0.95)


# ----------------------------------------------------------------- fig5a: anolyte
def fig_anolyte(rows):
    ano = [r for r in rows if r["pool"] == "anolyte" and not r["ctrl"]]
    ctrls = [r for r in rows if r["pool"] == "anolyte" and r["ctrl"]]
    front = sorted(pareto(ano, want_low_E=True), key=lambda r: r["E"])
    front_ids = {r["id"] for r in front}

    fig, ax = plt.subplots(figsize=(16.0, 8.4))
    fig.subplots_adjust(left=0.085, right=0.52, top=0.90, bottom=0.13)

    for r in ano:
        if r["id"] in front_ids:
            continue
        ax.scatter(r["E"], r["cap"], s=70, c=COL[r["fk"]], alpha=0.42,
                   edgecolors="none", zorder=2)
    for r in ctrls:  # 1e controls: hollow grey, excluded from the front
        ax.scatter(r["E"], r["cap"], s=80, facecolors="none", edgecolors="#777777",
                   linewidths=1.6, zorder=3)

    xs = [r["E"] for r in front]
    ys = [r["cap"] for r in front]
    ax.plot(xs, ys, "-", color=ANO_FRONT, lw=2.6, alpha=0.9, zorder=4)
    ax.errorbar(xs, ys, xerr=SIGMA_E, fmt="none", ecolor=ANO_FRONT, elinewidth=1.6,
                capsize=5, alpha=0.6, zorder=5)
    key_lines = []
    for i, r in enumerate(front, 1):
        ax.scatter(r["E"], r["cap"], s=440, c=COL[r["fk"]], edgecolors="black",
                   linewidths=1.6, marker="o", zorder=6)
        ax.annotate(str(i), (r["E"], r["cap"]), ha="center", va="center",
                    fontsize=BASE, fontweight="bold", color="white", zorder=7)
        key_lines.append(f"{i}. {r['lab']}  ({_mm(r['E'])} V, {r['cap']:.0f})")

    ax.set_xlim(-1.78, 0.06)
    ax.margins(y=0.12)
    ax.axvline(0.0, color="#555555", lw=1.1, ls="--", alpha=0.7, zorder=1)
    ax.text(0.993, 0.5, "catholyte →", transform=ax.transAxes, rotation=90,
            ha="right", va="center", fontsize=BASE, color="#999999")
    ax.set_xlabel(XLABEL)
    ax.set_ylabel(YLABEL)
    ax.set_title("Anolyte pool — 85 n-type candidates")

    fam_handles = [
        Line2D([0], [0], marker="o", ls="", mfc=COL["quinone"], mec="none", ms=13,
               label="quinone (AQ/NQ/BQ)"),
        Line2D([0], [0], marker="o", ls="", mfc=COL["imide"], mec="none", ms=13,
               label="imide (PMDI/NDI)"),
        Line2D([0], [0], marker="o", ls="", mfc=COL["viologen"], mec="none", ms=13,
               label="viologen (pyridine-multi-e)"),
        Line2D([0], [0], marker="o", ls="", mfc="none", mec="#777777", ms=13,
               label="1e control (excluded)"),
        Line2D([0], [0], marker="o", ls="-", color=ANO_FRONT, mfc="w", mec="black",
               ms=13, label="Pareto front"),
    ]
    _right_margin_blocks(fig, "Pareto front   " + WHISKER_LABEL,
                         key_lines, fam_handles, "family")
    return _save(fig, "fig5a_anolyte_pareto"), front, key_lines


# ----------------------------------------------------------------- fig5b: catholyte
def fig_catholyte(rows):
    cat = [r for r in rows if r["pool"] == "catholyte" and not r["ctrl"]]
    front = sorted(pareto(cat, want_low_E=False), key=lambda r: -r["E"])
    front_ids = {r["id"] for r in front}

    fig, ax = plt.subplots(figsize=(16.0, 8.4))
    fig.subplots_adjust(left=0.085, right=0.52, top=0.90, bottom=0.13)

    ax.set_xlim(-0.02, 2.0)   # full catholyte window (0 -> +2.0 V vs Fc) -> sparsity visible
    ax.set_ylim(74, 96)

    cat_max_E = max(r["E"] for r in cat)
    ax.axvspan(cat_max_E, 2.0, color="#DDDDDD", alpha=0.35, zorder=0)
    ax.text((cat_max_E + 2.0) / 2, 94.6, "sparse /\nlargely unexplored", ha="center",
            va="top", fontsize=BASE, color="#888888", style="italic")

    for r in cat:
        if r["id"] in front_ids:
            continue
        ax.scatter(r["E"], r["cap"], s=150, c=COL[r["fk"]], alpha=0.9,
                   edgecolors="none", zorder=2)

    xs = [r["E"] for r in front]
    ys = [r["cap"] for r in front]
    ax.plot(xs, ys, "-", color=CAT_FRONT, lw=2.6, alpha=0.9, zorder=4)
    ax.errorbar(xs, ys, xerr=SIGMA_E, fmt="none", ecolor=CAT_FRONT, elinewidth=1.6,
                capsize=5, alpha=0.6, zorder=5)
    key_lines = []
    for i, r in enumerate(front, 1):
        ax.scatter(r["E"], r["cap"], s=470, c=COL[r["fk"]], edgecolors="black",
                   linewidths=1.6, marker="s", zorder=6)
        ax.annotate(str(i), (r["E"], r["cap"]), ha="center", va="center",
                    fontsize=BASE, fontweight="bold", color="white", zorder=7)
        key_lines.append(f"{i}. {r['lab']}  ({_mm(r['E'])} V, {r['cap']:.0f})")

    ax.axvline(0.0, color="#555555", lw=1.1, ls="--", alpha=0.7, zorder=1)
    ax.set_xlabel(XLABEL)
    ax.set_ylabel(YLABEL)
    ax.set_title("Catholyte pool — 6 p-type candidates")
    ax.text(0.03, 0.045, "6 candidates vs 85 anolytes\n→ catholyte space is the gap",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=BASE,
            color=CAT_FRONT, fontweight="bold")

    fam_handles = [
        Line2D([0], [0], marker="o", ls="", mfc=COL["phenazine"], mec="none", ms=13,
               label="DMPZ (dihydrophenazine)"),
        Line2D([0], [0], marker="o", ls="", mfc=COL["phenothiazine"], mec="none", ms=13,
               label="PTZ (phenothiazine)"),
        Line2D([0], [0], marker="s", ls="-", color=CAT_FRONT, mfc="w", mec="black",
               ms=13, label="Pareto front"),
    ]
    _right_margin_blocks(fig, "Pareto front   " + WHISKER_LABEL,
                         key_lines, fam_handles, "family")
    return _save(fig, "fig5b_catholyte_pareto"), front, key_lines


def main():
    apply_style()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    rows, dropped = load_frozen()

    n_ano = sum(r["pool"] == "anolyte" for r in rows)
    n_cat = sum(r["pool"] == "catholyte" for r in rows)
    fams = {}
    for r in rows:
        fams[r["family"]] = fams.get(r["family"], 0) + 1
    print(f"frozen pool: {len(rows)} candidates  ({n_ano} anolyte / {n_cat} catholyte)")
    print("families:", fams)
    if dropped:
        print("DROPPED (live-campaign extra, not in frozen 91):",
              [(r["id"], r["name"], r["E_mean"]) for r in dropped])
    assert len(rows) == 91, f"expected frozen 91, got {len(rows)}"
    assert (n_ano, n_cat) == (85, 6), f"expected 85/6, got {n_ano}/{n_cat}"

    (pa, fa, ka) = fig_anolyte(rows)
    (pc, fc, kc) = fig_catholyte(rows)
    print("\nANOLYTE FRONT:", len(fa))
    for line in ka:
        print("  ", line)
    print("\nCATHOLYTE FRONT:", len(fc))
    for line in kc:
        print("  ", line)
    print("\nwrote:", pa, pc)


if __name__ == "__main__":
    main()
