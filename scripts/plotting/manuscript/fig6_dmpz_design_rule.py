#!/usr/bin/env python
"""Manuscript Figure 6 - grafted 5,10-dihydrophenazine (DMPZ) catholyte design rule.

A matched single-substituent expansion of one scaffold (DMPZ, 4-methylbenzyl /
Merrifield-proxy cap, two-electron {neutral -> +1 -> +2} ladder) shows two
independent, oppositely-signed handles on the mean redox potential vs Fc/Fc+:

  * an electron-WITHDRAWING N'-substituent RAISES E (N'-2,2,2-trifluoroethyl
    sits +0.20 V above the N'-ethyl reference);
  * an electron-DONATING ring substituent LOWERS E (2,7-dimethyl -0.13 V,
    2,7-dimethoxy -0.22 V vs the same N'-ethyl parent).

The two handles shift the potential in opposite directions. Whether they add when
combined in one molecule is untested (no member carries both).

fig6  : ordered potential ladder of the DMPZ series, Delta-vs-N'-ethyl-reference
        annotations, electron-withdrawing (N') vs electron-donating (ring)
        effects grouped with arrows, and the distinct 3,7-dimethoxyphenothiazine
        (PTZ) as a reach extension.
fig6b : the catholyte sub-space (E vs Fc/Fc+ against usable capacity) with the
        same DMPZ family highlighted, showing the gap filled by a structured
        family rather than isolated outliers.

Style contract (manuscript/figures/STYLE_CONTRACT.md): PNG only, 300 dpi, 18 pt
floor on EVERY text element, zero overlaps, family colours fixed across figures
(DMPZ orange #E69F00, PTZ purple #762A83). CPU / matplotlib only. No number is
altered - potentials are read live from the scored CSV.

All potentials are at the pool's ranking thermal tier (GFN2-xTB RRHO), like every
value in the scored pool; none of the DMPZ members has been confirmed with the
native-Hessian tier, so no member is marked differently from the others.

  python scripts/plotting/manuscript/fig6_dmpz_design_rule.py
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_style import apply_style, C, grid_x, grid_xy  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import FancyArrowPatch  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
CSV = ROOT / "results" / "validation" / "full_pool_scored_legacyG_MeCN_TBAPF6.csv"
OUT = ROOT / "manuscript" / "figures"

FS = 18  # pt - hard floor for every text element (STYLE_CONTRACT sec. 7)
REF_ID = "exp_dmpz_mb_ethyl"  # N'-ethyl DMPZ: the matched reference (E = +0.43 V)

# Family colours, fixed across the whole figure set (STYLE_CONTRACT sec. 2).
DMPZ_C = "#E69F00"   # phenazine / DMPZ (p-type) - orange
PTZ_C = "#762A83"    # phenothiazine / PTZ (p-type) - purple
RAISE_C = C["warn"]  # red   - direction annotation: substituent raises E
LOWER_C = C["good"]  # green - direction annotation: substituent lowers E

# Short, real-chemistry labels (no invented shorthand) + mechanistic role.
# role: "ref" reference | "raise" N'-EWG raises E | "lower" ring donor lowers E |
#       "ptz" distinct phenothiazine scaffold (not part of the DMPZ delta rule)
META = {
    "exp_dmpz_mb_ethyl":         dict(label="N'-ethyl (reference)",      fam="dmpz", role="ref",   provisional=False),
    "exp_dmpz_mb_methyl":        dict(label="N'-methyl",                 fam="dmpz", role="raise", provisional=False),
    "hyp_dmpz_mb_methoxyethyl":  dict(label="N'-(2-methoxyethyl)",       fam="dmpz", role="raise", provisional=False),
    "hyp_dmpz_mb_trifluoroethyl":dict(label="N'-(2,2,2-trifluoroethyl)", fam="dmpz", role="raise", provisional=False),
    "hyp_dmpz_27dimethyl_ethyl": dict(label="2,7-dimethyl, N'-ethyl",    fam="dmpz", role="lower", provisional=False),
    "hyp_dmpz_27dimethoxy_ethyl":dict(label="2,7-dimethoxy, N'-ethyl",   fam="dmpz", role="lower", provisional=False),
    "exp_ptz_37dimethoxy":       dict(label="3,7-dimethoxyphenothiazine",fam="ptz",  role="ptz",   provisional=False),
}

FAM_C = {"dmpz": DMPZ_C, "ptz": PTZ_C}


def _fmt(v):
    """Signed potential with a true unicode minus, 2 decimals."""
    s = f"{abs(v):.2f}"
    return ("+" if v >= 0 else "−") + s


def load_rows():
    rows = {r["id"]: r for r in csv.DictReader(open(CSV))}
    out = {}
    for cid, m in META.items():
        r = rows[cid]
        out[cid] = dict(E=float(r["E_mean"]), cap=float(r["cap_active_mAh_g"]),
                        couples=r["E_couples"], **m)
    return out, rows


def _marker_kw(d):
    """Filled for confirmed anchors, open (white face) for provisional members."""
    color = FAM_C[d["fam"]]
    if d["provisional"]:
        return dict(facecolor="white", edgecolor=color, linewidth=3.0)
    return dict(facecolor=color, edgecolor=color, linewidth=3.0)


# ---------------------------------------------------------------------------
def make_fig6(data):
    """Ordered potential ladder with matched Delta-vs-reference annotations."""
    ref_E = data[REF_ID]["E"]
    # order bottom -> top by ascending E (ladder reads low -> high catholyte)
    items = sorted(data.items(), key=lambda kv: kv[1]["E"])
    ids = [k for k, _ in items]
    dy = 1.25
    ypos = {cid: i * dy for i, cid in enumerate(ids)}
    ytop = (len(ids) - 1) * dy

    fig, ax = plt.subplots(figsize=(16.0, 12.0))

    # reference potential guide
    # drawn only across the data rows so it never crosses the arrow captions
    ax.plot([ref_E, ref_E], [-0.55, ytop + 0.55], color=C["measured"], lw=1.8, ls="--",
            zorder=1, alpha=0.8)

    for cid in ids:
        d = data[cid]
        y = ypos[cid]
        fam_c = FAM_C[d["fam"]]
        # connector from the reference to the member, coloured by shift direction
        if d["role"] in ("raise", "lower"):
            dir_c = RAISE_C if d["role"] == "raise" else LOWER_C
            ax.plot([ref_E, d["E"]], [y, y], color=dir_c, lw=4.0, alpha=0.55,
                    solid_capstyle="round", zorder=2)
        ax.scatter([d["E"]], [y], s=440, marker="o", zorder=4, **_marker_kw(d))

        # far-side annotation block: absolute potential + matched delta.
        if d["role"] == "raise":
            dE = d["E"] - ref_E
            txt = f"{_fmt(d['E'])} V\n({_fmt(dE)} vs ref)"
            ax.annotate(txt, (d["E"], y), xytext=(16, 0), textcoords="offset points",
                        ha="left", va="center", fontsize=FS, color=RAISE_C)
        elif d["role"] == "lower":
            dE = d["E"] - ref_E
            txt = f"{_fmt(d['E'])} V\n({_fmt(dE)} vs ref)"
            ax.annotate(txt, (d["E"], y), xytext=(-16, 0), textcoords="offset points",
                        ha="right", va="center", fontsize=FS, color=LOWER_C)
        elif d["role"] == "ref":
            ax.annotate(f"{_fmt(d['E'])} V\n(reference)", (d["E"], y),
                        xytext=(20, 0), textcoords="offset points",
                        ha="left", va="center", fontsize=FS, color=C["measured"])
        else:  # ptz
            ax.annotate(f"{_fmt(d['E'])} V\n(PTZ)", (d["E"], y),
                        xytext=(16, 0), textcoords="offset points",
                        ha="left", va="center", fontsize=FS, color=PTZ_C)

    # substituent names as y-tick labels (18 pt); mark provisional with *
    ax.set_yticks([ypos[c] for c in ids])
    ax.set_yticklabels([data[c]["label"] + (" *" if data[c]["provisional"] else "")
                        for c in ids], fontsize=FS)

    # directional design-rule arrows in the top / bottom margins (tied to the x axis)
    ax.add_patch(FancyArrowPatch((ref_E, ytop + 0.78), (0.86, ytop + 0.78),
                 arrowstyle="-|>", mutation_scale=28, lw=3.2, color=RAISE_C, zorder=3))
    ax.text((ref_E + 0.86) / 2, ytop + 1.02,
            "electron-withdrawing N'-substituent $\\rightarrow$ raises E",
            ha="center", va="bottom", fontsize=FS, color=RAISE_C)
    ax.add_patch(FancyArrowPatch((ref_E, -1.05), (0.12, -1.05),
                 arrowstyle="-|>", mutation_scale=28, lw=3.2, color=LOWER_C, zorder=3))
    ax.text((ref_E + 0.12) / 2, -1.30,
            "electron-donating ring substituent $\\rightarrow$ lowers E",
            ha="center", va="top", fontsize=FS, color=LOWER_C)

    ax.set_ylim(-1.75, ytop + 1.55)
    ax.set_xlim(-0.12, 1.16)
    ax.set_xlabel("Redox potential, E vs Fc/Fc⁺ (V)", fontsize=FS)
    ax.set_title("DMPZ catholyte design rule: two opposite-sign handles",
                 fontsize=FS + 2, pad=18)
    ax.tick_params(axis="x", labelsize=FS)
    ax.tick_params(axis="y", length=0)
    grid_x(ax, alpha=0.35)

    legend = [
        Line2D([], [], marker="o", ls="none", mfc=DMPZ_C, mec=DMPZ_C, ms=15,
               label="DMPZ (dihydrophenazine)"),
        Line2D([], [], marker="o", ls="none", mfc=PTZ_C, mec=PTZ_C, ms=15,
               label="phenothiazine (PTZ)"),
    ]
    ax.legend(handles=legend, loc="lower right", fontsize=FS, framealpha=0.95,
              borderpad=0.6)

    fig.tight_layout()
    fig.savefig(OUT / "fig6_dmpz_design_rule.png")
    plt.close(fig)


# ---------------------------------------------------------------------------
def make_fig6b(data):
    """The catholyte sub-space: E vs Fc/Fc+ against usable capacity."""
    fig, ax = plt.subplots(figsize=(15.0, 11.0))

    for cid, d in data.items():
        ax.scatter([d["E"]], [d["cap"]], s=460, marker="o", zorder=4, **_marker_kw(d))

    # label placement (dx_pts, dy_pts, ha) chosen so no two labels / markers collide
    lab = {
        "exp_dmpz_mb_ethyl":          (-16, 2, "right"),
        "exp_dmpz_mb_methyl":         (0, 30, "center"),
        "hyp_dmpz_mb_methoxyethyl":   (18, 2, "left"),
        "hyp_dmpz_mb_trifluoroethyl": (0, -34, "center"),
        "hyp_dmpz_27dimethyl_ethyl":  (0, 30, "center"),
        "hyp_dmpz_27dimethoxy_ethyl": (0, 30, "center"),
        "exp_ptz_37dimethoxy":        (-16, 10, "right"),
    }
    for cid, d in data.items():
        dx, dy, ha = lab[cid]
        star = " *" if d["provisional"] else ""
        ax.annotate(d["label"] + star, (d["E"], d["cap"]),
                    xytext=(dx, dy), textcoords="offset points", ha=ha,
                    va="bottom" if dy >= 0 else "top", fontsize=FS,
                    color=FAM_C[d["fam"]])

    # span annotation across the grafted DMPZ series
    dmpz_E = [d["E"] for d in data.values() if d["fam"] == "dmpz"]
    lo, hi = min(dmpz_E), max(dmpz_E)
    caps = [d["cap"] for d in data.values()]
    y_arrow = max(caps) + 9.0
    ax.annotate("", xy=(hi, y_arrow), xytext=(lo, y_arrow),
                arrowprops=dict(arrowstyle="<->", color="#555555", lw=2.0))
    ax.text((lo + hi) / 2, y_arrow + 1.0,
            f"grafted DMPZ series: {_fmt(lo)} to {_fmt(hi)} V",
            ha="center", va="bottom", fontsize=FS, color="#333333")

    ax.set_xlabel("Redox potential, E vs Fc/Fc⁺ (V)", fontsize=FS)
    ax.set_ylabel("Usable capacity (mAh g⁻¹, active mass)", fontsize=FS)
    ax.set_title("Catholyte sub-space populated by a structured DMPZ family",
                 fontsize=FS + 2, pad=14)
    ax.tick_params(axis="both", labelsize=FS)
    ax.set_xlim(0.05, 0.98)
    ax.set_ylim(min(caps) - 9.0, max(caps) + 15.0)
    grid_xy(ax, alpha=0.3)

    legend = [
        Line2D([], [], marker="o", ls="none", mfc=DMPZ_C, mec=DMPZ_C, ms=15,
               label="DMPZ (dihydrophenazine)"),
        Line2D([], [], marker="o", ls="none", mfc=PTZ_C, mec=PTZ_C, ms=15,
               label="phenothiazine (PTZ)"),
    ]
    ax.legend(handles=legend, loc="upper right", fontsize=FS, framealpha=0.95,
              borderpad=0.6)

    fig.tight_layout()
    fig.savefig(OUT / "fig6b_dmpz_catholyte_space.png")
    plt.close(fig)


def main():
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    data, _ = load_rows()
    make_fig6(data)
    make_fig6b(data)
    print("wrote (PNG only):")
    for n in ("fig6_dmpz_design_rule.png", "fig6b_dmpz_catholyte_space.png"):
        print("  ", OUT / n)


if __name__ == "__main__":
    main()
