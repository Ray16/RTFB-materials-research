#!/usr/bin/env python
"""Workflow diagram: how a molecule flows through the redox-screening pipeline.

Single panel, hand-laid on a 0-1 canvas so nothing overlaps. Every stage box
carries its method/level-of-theory; the two inputs (candidate monomers and the
known-E validation set) share the *identical* pipeline, and the validation gate
is what licenses the final ranking. Read left column top->bottom.

  python scripts/plot_pipeline.py   ->  results/figures/pipeline.png
"""
from __future__ import annotations
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_style import apply_style, C  # noqa: E402

import matplotlib.pyplot as plt          # noqa: E402
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIGDIR = ROOT / "results" / "figures"

TITLE_FS = 26
STAGE_FS = 20      # stage title (bold)
DETAIL_FS = 18     # method/detail line (18 pt floor)
PHASE_FS = 18

# ------------------------------------------------------------------ geometry
CX = 0.50          # main column centre
CW = 0.60          # card width
CH = 0.098         # card height
XL = CX - CW / 2   # card left edge
XR = CX + CW / 2   # card right edge


def card(ax, cy, title, detail, fc, ec, *, h=CH, w=CW, cx=CX,
         title_fs=STAGE_FS, detail_fs=DETAIL_FS):
    """Rounded stage box with a bold title line and a detail line, centred at (cx, cy)."""
    x = cx - w / 2
    box = FancyBboxPatch((x, cy - h / 2), w, h,
                         boxstyle="round,pad=0.006,rounding_size=0.014",
                         mutation_aspect=1.0, fc=fc, ec=ec, lw=2.2, zorder=3)
    ax.add_patch(box)
    if detail:
        ax.text(cx, cy + h * 0.20, title, ha="center", va="center",
                fontsize=title_fs, fontweight="bold", color="#111111", zorder=4)
        ax.text(cx, cy - h * 0.24, detail, ha="center", va="center",
                fontsize=detail_fs, color="#222222", zorder=4)
    else:
        ax.text(cx, cy, title, ha="center", va="center",
                fontsize=title_fs, fontweight="bold", color="#111111", zorder=4)
    return box


def arrow(ax, y0, y1, x=CX, label=None, color="#3A3A3A"):
    a = FancyArrowPatch((x, y0), (x, y1), arrowstyle="-|>", mutation_scale=26,
                        lw=2.6, color=color, shrinkA=0, shrinkB=0, zorder=2)
    ax.add_patch(a)
    if label:
        ax.text(x + 0.015, (y0 + y1) / 2, label, ha="left", va="center",
                fontsize=DETAIL_FS, color="#333333", zorder=5)


def phase_bracket(ax, y_top, y_bot, text, color):
    """A thin vertical phase label on the far left."""
    xb = XL - 0.075
    ax.plot([xb, xb], [y_bot, y_top], color=color, lw=3.2,
            solid_capstyle="round", zorder=2)
    for yy in (y_top, y_bot):
        ax.plot([xb, xb + 0.02], [yy, yy], color=color, lw=3.2,
                solid_capstyle="round", zorder=2)
    ax.text(xb - 0.018, (y_top + y_bot) / 2, text, ha="center", va="center",
            rotation=90, fontsize=PHASE_FS, fontweight="bold", color=color, zorder=5)


def main():
    apply_style()
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(13.5, 18.5))
    ax = fig.add_axes([0.02, 0.01, 0.96, 0.98])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    # ---- title
    ax.text(CX, 0.982, "Redox-potential screening pipeline",
            ha="center", va="top", fontsize=TITLE_FS, fontweight="bold")
    ax.text(CX, 0.952,
            "gas-phase ML pre-optimization  →  DFT geometry & energy in implicit solvent  →  "
            "E$^\\circ$ vs Fc/Fc$^+$",
            ha="center", va="top", fontsize=DETAIL_FS, color="#555555")

    # ---- two inputs sharing one pipeline
    iy = 0.895
    iw, ih = 0.375, 0.072
    card(ax, iy, "Candidate monomers", "benzylic-site redox groups",
         "#EDEDED", "#7A7A7A", w=iw, h=ih, cx=0.285,
         title_fs=19, detail_fs=17)
    card(ax, iy, "Validation set", "cores with measured E$^\\circ$",
         "#E6E0EF", "#8C6BA8", w=iw, h=ih, cx=0.715,
         title_fs=19, detail_fs=17)
    ax.text(CX, iy, "+", ha="center", va="center", fontsize=30,
            fontweight="bold", color="#555555", zorder=5)

    # every state carries explicit charge & spin
    merge_y = 0.815
    for sx in (0.285, 0.715):
        arrow(ax, iy - ih / 2, merge_y + 0.004, x=sx)
    ax.text(CX, 0.838, "each redox state  =  (charge $q$, spin multiplicity)  —  "
                       "never assume neutral singlet",
            ha="center", va="center", fontsize=DETAIL_FS, color="#333333",
            style="italic", zorder=6,
            bbox=dict(boxstyle="round,pad=0.35", fc="#FFF9E6", ec="#E0C36A", lw=1.4))

    # ---- main flow cards (centre column)
    ys = [0.740, 0.610, 0.480, 0.350, 0.220]
    card(ax, ys[0], "3D structure generation",
         "RDKit ETKDG conformer ensemble · stereo verified · all cores",
         "#DCE6F1", "#4B7BB0")
    card(ax, ys[1], "UMA ML pre-optimization  (gas phase)",
         "uma-s-1p2p1 · per (charge, spin) · GPU-batched",
         "#FBE3DE", C["uma"])
    card(ax, ys[2], "DFT geometry optimization  in SMD",
         "r2SCAN-D4 / def2-SVP(D) · acetonitrile ($\\epsilon$=37.5)",
         "#CBD9EC", C["dft"])
    card(ax, ys[3], "DFT single-point energy  in SMD",
         "$\\omega$B97M-V / def2-TZVP(D) + VV10 · diffuse basis for anions",
         "#CBD9EC", C["dft"])
    card(ax, ys[4], "Redox potential",
         "E$^\\circ$ = $-\\Delta G/nF$ $-$ E(Fc) · vs Fc/Fc$^+$",
         "#DCEBDD", C["good"])

    labels = ["conformer + charge/spin",
              "gas-phase geometry",
              "solvated geometry",
              "G in solvent",
              None]
    arrow(ax, merge_y, ys[0] + CH / 2, label=None)
    for i in range(len(ys) - 1):
        arrow(ax, ys[i] - CH / 2, ys[i + 1] + CH / 2, label=labels[i])

    # ---- validation gate + outputs
    gate_y = 0.095
    gw, gh = 0.72, 0.090
    card(ax, gate_y, "Validation gate:  computed  vs  measured E$^\\circ$",
         "MAE 0.18 V · within $\\pm$0.15 V band  →  trust the ranking",
         "#FDF0D5", "#C89A2B", w=gw, h=gh, title_fs=18.5, detail_fs=16)
    arrow(ax, ys[4] - CH / 2, gate_y + gh / 2)

    apply_footer(ax)
    fig.savefig(FIGDIR / "pipeline.png")
    print(f"  -> {FIGDIR/'pipeline.png'}")


def apply_footer(ax):
    ax.text(CX, 0.028,
            "No fitting: computed potentials are never rescaled to experiment; "
            "the gate exposes error, it does not hide it.",
            ha="center", va="center", fontsize=14.5, style="italic", color="#666666")


if __name__ == "__main__":
    main()
