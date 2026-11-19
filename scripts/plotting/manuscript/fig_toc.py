#!/usr/bin/env python
"""Table-of-contents graphic (JCIM: 3.25 in x 1.75 in at final size).

Drawn at 2x the final size (6.5 in x 3.5 in) with 18-pt text, so text prints at 9 pt.
Left: the pipeline in four steps. Right: potential vs capacity for the 91-candidate pool
(anolytes grey) with the grafted 5,10-dihydrophenazine (DMPZ) catholytes highlighted, read
from the scored table, and one proposed cell drawn as an arrow. PNG only.

  python scripts/plotting/manuscript/fig_toc.py
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_style import apply_style  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
CSV = ROOT / "data/manuscript/full_pool_scored_legacyG_MeCN_TBAPF6.csv"
OUT = ROOT / "manuscript/figures/toc_graphic.png"
FS = 18
DMPZ_C, ANO_C, PTZ_C = "#E69F00", "#9AA3AD", "#762A83"
CELL = ("exp_dmpz_mb_ethyl", "exp_aq_2oxy_6methoxy")   # 1.91 V pairing used in the text


def main():
    apply_style()
    rows = list(csv.DictReader(open(CSV)))
    by = {r["id"]: r for r in rows}
    fig = plt.figure(figsize=(6.5, 3.5))
    # ---- left: pipeline steps
    axl = fig.add_axes([0.0, 0.0, 0.40, 1.0]); axl.axis("off")
    axl.set_xlim(0, 1); axl.set_ylim(0, 1)
    steps = ["grafted motif", "DFT in solvent", "quality checks", "potential vs Fc"]
    ys = [0.86, 0.62, 0.38, 0.14]
    for i, (s, y) in enumerate(zip(steps, ys)):
        axl.add_patch(FancyBboxPatch((0.06, y - 0.085), 0.88, 0.17, boxstyle="round,pad=0.01",
                                     fc="#EEF1F4", ec="#555555", lw=1.2))
        axl.text(0.5, y, s, ha="center", va="center", fontsize=FS)
        if i < len(steps) - 1:
            axl.add_patch(FancyArrowPatch((0.5, y - 0.09), (0.5, ys[i + 1] + 0.09),
                                          arrowstyle="-|>", mutation_scale=14, lw=1.4, color="#555555"))
    # ---- right: design map
    ax = fig.add_axes([0.53, 0.20, 0.45, 0.74])
    ano = [r for r in rows if r["pool"] == "anolyte" and r["n"] == "2"]
    dm = [r for r in rows if "dmpz" in r["id"]]
    ax.scatter([float(r["E_mean"]) for r in ano], [float(r["cap_active_mAh_g"]) for r in ano],
               s=14, c=ANO_C, edgecolors="none", label="anolytes")
    ax.scatter([float(r["E_mean"]) for r in dm], [float(r["cap_active_mAh_g"]) for r in dm],
               s=46, c=DMPZ_C, edgecolors="black", linewidths=0.6, label="DMPZ")
    c, a = (by[i] for i in CELL)
    ec, ea = float(c["E_mean"]), float(a["E_mean"])
    y_arrow = 222
    ax.add_patch(FancyArrowPatch((ea, y_arrow), (ec, y_arrow), arrowstyle="<|-|>",
                                 mutation_scale=12, lw=1.6, color="#B2182B"))
    ax.text((ea + ec) / 2, y_arrow + 6, f"cell ≈ {ec - ea:.1f} V", ha="center", va="bottom",
            fontsize=FS, color="#B2182B")
    ax.text(0.42, 120, "DMPZ", ha="center", va="top", fontsize=FS, color="#8a5a00")
    ax.set_xlim(-1.8, 1.0); ax.set_ylim(80, 250)
    ax.set_xticks([-1.5, -0.5, 0.5])
    ax.set_yticks([])
    ax.set_xlabel("E vs Fc (V)", fontsize=FS)
    ax.set_ylabel("capacity", fontsize=FS)
    ax.tick_params(labelsize=FS)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.savefig(OUT, dpi=600)
    plt.close(fig)
    print("wrote", OUT, f"(cell {ec - ea:.2f} V: {c['name']} / {a['name']})")


if __name__ == "__main__":
    main()
