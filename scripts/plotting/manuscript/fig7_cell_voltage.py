#!/usr/bin/env python
"""Manuscript Figure 7 - proposed full-cell voltage ladder (pairing diagram).

A single-panel energy ladder on a potential axis (V vs Fc/Fc+). The two separate
screening pools are shown as bands: the catholyte pool (p-type DMPZ /
phenothiazine, orange, near the top) and the anolyte pool (n-type quinone /
imide / pyridine, blue, near the bottom). A handful of REAL, verifiable
candidate pairings are drawn as vertical arrows annotated with the resulting
nominal cell voltage, E_cell = E(catholyte) - E(anolyte).

HARD CAVEAT: every arrow is a COMPUTED PROPOSAL from two validated half-cell
potentials. It is NOT a demonstrated device - nothing here is measured capacity,
cycling, Coulombic efficiency, membrane crossover, or synthetic feasibility.

Style contract (manuscript/figures/STYLE_CONTRACT.md): PNG only, 18 pt floor on
EVERY text element, zero overlaps, DMPZ orange / anolyte blue, true minus signs,
Fc/Fc+ reference, abbreviations unified (DMPZ, PTZ, AQ, PMDI).

CPU / matplotlib only. No GPU, no DFT, no QM. Reads the frozen scored pool CSV
and does nothing but arithmetic on its E_mean column.
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "scripts" / "plotting"))
from plot_style import apply_style, OKABE_ITO  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402

CSV = REPO / "data/manuscript/full_pool_scored_legacyG_MeCN_TBAPF6.csv"
OUT = REPO / "manuscript/figures/fig7_cell_voltage"

# Family colours (STYLE_CONTRACT section 2).
DMPZ_ORANGE = OKABE_ITO["orange"]      # #E69F00  phenazine / DMPZ (p-type)
DMPZ_TEXT = "#8A5A00"                   # darker orange for legible 18 pt text
PTZ_PURPLE = "#762A83"                 # phenothiazine (p-type)
ANOL_BLUE = OKABE_ITO["blue"]          # #0072B2  n-type anolyte pool
ANOL_TEXT = "#0A4A73"
WARN = "#B2182B"
FS = 18                                 # hard floor - every text element

# Compact, real-chemistry labels for the catholyte levels (keyed by candidate id).
CAT_LABEL = {
    "exp_ptz_37dimethoxy": "3,7-dimethoxy PTZ",
    "hyp_dmpz_mb_trifluoroethyl": "N'-CF$_3$CH$_2$ DMPZ",
    "hyp_dmpz_mb_methoxyethyl": "N'-(2-methoxyethyl) DMPZ",
    "exp_dmpz_mb_methyl": "N'-methyl DMPZ",
    "exp_dmpz_mb_ethyl": "N'-ethyl DMPZ",
    "hyp_dmpz_27dimethyl_ethyl": "2,7-dimethyl DMPZ",
    "hyp_dmpz_27dimethoxy_ethyl": "2,7-dimethoxy DMPZ",
}

# Representative, real pairings (ids verified present in the CSV). Cell voltages
# are NOT hard-coded - they are computed below as E_cat - E_an so the annotation
# can never drift from the data.  (cath_id, anol_id, anol_tag, highlight, star)
PAIRINGS = [
    ("exp_dmpz_mb_ethyl", "exp_aq_2oxy_6methoxy", "alkoxy-AQ", True, False),
    ("exp_dmpz_mb_ethyl", "search_pmdi_87840c94153a", "ether-PMDI", False, False),
    ("hyp_dmpz_mb_trifluoroethyl", "search_pmdi_87840c94153a", "ether-PMDI", False, True),
]


def _mv(v: float) -> str:
    """Signed potential, 2 dp, with a true minus sign."""
    return f"{v:+.2f} V".replace("-", "−")


def main() -> None:
    apply_style()
    df = pd.read_csv(CSV).set_index("id")

    cath = df[df["pool"].str.contains("cath", case=False, na=False)].sort_values(
        "E_mean", ascending=False)
    anol = df[df["pool"].str.contains("anol", case=False, na=False)].sort_values("E_mean")

    fig, ax = plt.subplots(figsize=(17.0, 12.0))
    xmax = 16.5
    ax.set_xlim(0, xmax)
    ax.set_ylim(-1.98, 1.10)

    # --- pool bands -------------------------------------------------------
    c_lo, c_hi = cath["E_mean"].min(), cath["E_mean"].max()
    a_lo, a_hi = anol["E_mean"].min(), anol["E_mean"].max()
    ax.axhspan(c_lo - 0.03, c_hi + 0.03, color=DMPZ_ORANGE, alpha=0.13, zorder=0)
    ax.axhspan(a_lo - 0.03, a_hi + 0.03, color=ANOL_BLUE, alpha=0.11, zorder=0)

    # Fc/Fc+ reference
    ax.axhline(0.0, color="#555555", ls="--", lw=1.4, zorder=1)
    ax.text(0.25, -0.075, "Fc/Fc$^+$  (0 V)", ha="left", va="center",
            fontsize=FS, color="#555555")

    # pool headers
    ax.text(0.15, 1.015, "catholyte pool  (p-type: DMPZ / phenothiazine)",
            ha="left", va="center", fontsize=FS, color=DMPZ_TEXT, fontweight="bold")
    ax.text(0.15, a_lo - 0.14, "anolyte pool  (n-type: quinone / imide / pyridine, 85)",
            ha="left", va="center", fontsize=FS, color=ANOL_TEXT, fontweight="bold")

    # --- anolyte rug + range ----------------------------------------------
    rug_x0, rug_x1 = 0.4, 2.3
    for e in anol["E_mean"]:
        ax.plot([rug_x0, rug_x1], [e, e], color=ANOL_BLUE, alpha=0.32, lw=1.5, zorder=2)
    ax.text(rug_x1 + 0.25, (a_lo + a_hi) / 2,
            f"{_mv(a_lo)} to {_mv(a_hi)}", ha="left", va="center",
            fontsize=FS, color=ANOL_TEXT)

    # --- catholyte levels with decluttered leader-line labels --------------
    lvl_x0, lvl_x1 = 0.4, 1.7
    lbl_x = 1.9
    slots = np.linspace(0.88, 0.04, len(cath))   # >=0.14 V apart -> no 18 pt overlap
    for (cid, row), ys in zip(cath.iterrows(), slots):
        e = row["E_mean"]
        is_ptz = cid == "exp_ptz_37dimethoxy"
        line_c = PTZ_PURPLE if is_ptz else DMPZ_ORANGE
        txt_c = PTZ_PURPLE if is_ptz else DMPZ_TEXT
        ax.plot([lvl_x0, lvl_x1], [e, e], color=line_c, lw=3.0, zorder=3,
                solid_capstyle="round")
        ax.plot([lvl_x1, lbl_x - 0.07], [e, ys], color=line_c, lw=1.0,
                alpha=0.65, zorder=3)
        ax.text(lbl_x, ys, f"{CAT_LABEL.get(cid, cid)}  ({_mv(e)})",
                ha="left", va="center", fontsize=FS, color=txt_c)

    # --- pairing arrows ----------------------------------------------------
    arrow_x = [9.2, 11.9, 14.6]
    computed = []
    for (cid, aid, an_tag, hi, star), x in zip(PAIRINGS, arrow_x):
        e_cat = float(df.loc[cid, "E_mean"])
        e_an = float(df.loc[aid, "E_mean"])
        v = e_cat - e_an
        computed.append((cid, aid, e_cat, e_an, v))

        arr_c = PTZ_PURPLE if star else "#555555"
        ax.plot([x - 0.3, x + 0.3], [e_cat, e_cat], color=DMPZ_ORANGE, lw=3.4, zorder=4)
        ax.plot([x - 0.3, x + 0.3], [e_an, e_an], color=ANOL_BLUE, lw=3.4, zorder=4)
        ax.annotate("", xy=(x, e_cat), xytext=(x, e_an),
                    arrowprops=dict(arrowstyle="<|-|>", color=arr_c, lw=2.8,
                                    shrinkA=0, shrinkB=0), zorder=5)
        # voltage box on the arrow (white fill covers the line)
        mid = (e_cat + e_an) / 2
        lbl = f"{v:.2f} V" + ("  ★" if star else "")
        fc = "#FFF3CF" if hi else "white"
        ec = PTZ_PURPLE if star else (DMPZ_TEXT if hi else "#555555")
        ax.text(x, mid, lbl, ha="center", va="center", fontsize=FS, fontweight="bold",
                color="#222222",
                bbox=dict(boxstyle="round,pad=0.3", fc=fc, ec=ec, lw=1.8), zorder=6)
        # partner tags: DMPZ above the top cap, anolyte below the bottom cap
        ax.text(x, e_cat + 0.055, CAT_LABEL[cid], ha="center", va="bottom",
                fontsize=FS, color=DMPZ_TEXT, zorder=6)
        ax.text(x, e_an - 0.055, an_tag, ha="center", va="top",
                fontsize=FS, color=ANOL_TEXT, zorder=6)

    # formula header over the arrow region
    ax.text((arrow_x[0] + arrow_x[-1]) / 2, 1.015,
            "proposed pairings:  $E_{cell}=E_{cat}-E_{an}$",
            ha="center", va="center", fontsize=FS, fontweight="bold", color="#333333")

    # on-figure caveat, placed in the empty central band
    ax.text(5.6, -0.35,
            "COMPUTED PROPOSALS\n(not demonstrated capacity, cycling,\n"
            "crossover, or synthetic feasibility)",
            ha="center", va="center", fontsize=FS, style="italic", color=WARN,
            bbox=dict(boxstyle="round,pad=0.4", fc="#FBEEEE", ec=WARN, lw=1.4),
            zorder=6)

    ax.set_ylabel("Redox potential, $E$ vs Fc/Fc$^+$ (V)", fontsize=FS)
    ax.set_yticks(np.arange(-1.5, 1.01, 0.5))
    ax.set_yticklabels([_mv(t).replace(" V", "") for t in np.arange(-1.5, 1.01, 0.5)],
                       fontsize=FS)
    ax.set_xticks([])
    ax.spines["bottom"].set_visible(False)
    ax.set_axisbelow(True)
    ax.grid(True, axis="y", color="#C2C2C2", lw=0.7, alpha=0.35)

    fig.tight_layout()
    fig.savefig(f"{OUT}.png", dpi=300)
    plt.close(fig)

    print("Pairings (E_cell = E_cat - E_an):")
    for cid, aid, ec, ea, v in computed:
        print(f"  {cid} ({ec:+.4f}) x {aid} ({ea:+.4f}) = {v:.4f} V")
    print(f"wrote {OUT}.png")


if __name__ == "__main__":
    main()
