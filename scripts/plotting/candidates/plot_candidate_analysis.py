#!/usr/bin/env python
"""Candidate redox landscape (results/figures/candidates/candidate_E_redox_landscape.png):
computed E deg vs Fc/Fc+ per candidate (sorted), with sigma_E error bars and specific capacity
annotated — the "where do they sit" plot. (The former lambda-vs-SA panel is archived:
archive/scripts/plotting/candidates/plot_candidate_analysis.py; SA is not a ranking objective.)

Physics, not fitting — cosmetics only; no number is altered.

  results/scorecard.csv            per-candidate E, lambda(+-sigma), SA, capacity, family
  results/pareto_shortlist.csv     pool + pipeline's own pareto_optimal flags (cross-check)

  python scripts/plotting/candidates/plot_candidate_analysis.py
"""
from __future__ import annotations
import csv
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_style import apply_style, C, grid_xy, grid_x  # noqa: E402

import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
FIGDIR = ROOT / "results" / "figures" / "candidates"
FIGDIR.mkdir(parents=True, exist_ok=True)

# Family colours consistent with candidates_2x2 (viologen pink, imide green, quinone blue).
# Family owns the FILL colour (kept stable across the whole figure set). The second
# categorical dimension - which batch a candidate came from - is deliberately NOT given its
# own fill hue: that would overload colour with two meanings. It gets composite encoding
# instead (marker shape + an accent ring), per the categorical-palette rules.
# Validated with the palette checker: 4 fills + accent, all inside the lightness band, above
# the chroma floor, normal-vision dE >= 16.4, CVD worst-adjacent dE 7.6 (legal in the 6-8
# floor band precisely BECAUSE shape + direct labels provide the secondary encoding).
FAM = {
    "pyridine-multi-e": "#CC79A7",
    "imide (n-type)":   "#009E73",
    "quinone (n-type)": "#0072B2",
    "phenol (p-type)":  "#762A83",
}
BATCH_ACCENT = "#D55E00"          # ring colour marking the new merrifield_multi cohort
BATCH_MARKER = {"starting": "o", "merrifield_multi": "D"}
BATCH_LABEL  = {"starting": "starting candidates", "merrifield_multi": "new (Merrifield multi-e)"}
FAM_LABEL = {
    "pyridine-multi-e": "viologen",
    "imide (n-type)":   "imide",
    "quinone (n-type)": "quinone",
    "phenol (p-type)":  "phenol",
}
# Short display names for annotation.
SHORT = {
    "viologen": "benzyl-methyl viologen",
    "ethylviologen": "benzyl-ethyl viologen",
    "pmdi": "PMDI",
    "ndi_ammonium": "ammonium-NDI",
    "mophquinone": "MeO-phenyl quinone",
    "dmophquinone": "(MeO)$_2$-phenyl quinone",
    # batch 2
    "bisviologen": "bis-benzyl viologen",
    "aq_benzyloxy": "BnO-anthraquinone",
    "nq_benzyloxy": "BnO-naphthoquinone",
    "aq_benzylamino": "BnNH-anthraquinone",
    "dtbc_phenol": "BnO-DTB-phenol",
}


def _f(row, k):
    v = (row.get(k) or "").strip()
    try:
        return float(v)
    except ValueError:
        return None


def load_candidates():
    p = ROOT / "results" / "scorecard.csv"
    rows = [r for r in csv.DictReader(p.open()) if r.get("status") == "candidate"]
    out = []
    for r in rows:
        out.append(dict(
            id=r["id"], family=r["family"], role=r.get("role", ""),
            E=_f(r, "E_anolyte_V"), sE=_f(r, "sigma_E_V"),
            lam=_f(r, "lambda_i_eV"), slam=_f(r, "sigma_lambda_eV"),
            lam_qc=(r.get("lambda_qc") or ""), batch=(r.get("batch") or ""),
            SA=_f(r, "SA_score"), cap=_f(r, "specific_capacity_mAh_g"),
            n=int(float(r.get("n_accessible") or 0)),
        ))
    return out


def plot_redox_landscape(cands):
    apply_style()
    fig, ax = plt.subplots(figsize=(11.5, 7.5))

    # anolyte landscape: a catholyte-only candidate has no anolyte potential -> not plotted here
    cs = sorted((c for c in cands if c["E"] is not None), key=lambda c: c["E"])  # most reducing at bottom
    y = np.arange(len(cs))
    E = np.array([c["E"] for c in cs])
    sE = np.array([c["sE"] or 0.0 for c in cs])
    cols = [FAM[c["family"]] for c in cs]

    ax.barh(y, E, height=0.58, color=cols, alpha=0.85, edgecolor="#222222", lw=1.2, zorder=3)
    ax.errorbar(E, y, xerr=sE, fmt="none", ecolor="#333333", elinewidth=1.8,
                capsize=5, capthick=1.8, zorder=4)

    xr = 0.30   # right margin (positive side) reserved for the capacity annotations
    for yi, c in zip(y, cs):
        # capacity + electron count in the clear right margin (bars are negative)
        ax.text(0.06, yi, f"{c['cap']:.0f} mAh/g · {c['n']}e",
                ha="left", va="center", fontsize=13, color="#444444")

    ax.axvline(0, color="#888888", lw=1.2, ls=":")
    ax.set_yticks(y)
    ax.set_yticklabels([SHORT.get(c["id"], c["id"]) for c in cs], fontsize=15)
    ax.set_xlabel("computed E$^\\circ$  (V vs Fc/Fc$^+$, DFT+SMD in MeCN)"
                  "\nmore negative $\\rightarrow$ higher cell voltage")
    ax.set_title("Candidate redox landscape — anolyte potentials", pad=14)
    ax.set_xlim(min(E) - 0.15, xr)
    ax.set_ylim(-0.7, len(cs) - 0.3)
    grid_x(ax)

    from matplotlib.lines import Line2D
    fam_handles = [Line2D([0], [0], marker="s", ls="", ms=13, mfc=FAM[f], mec="#222222",
                          label=FAM_LABEL[f]) for f in FAM if f in [c["family"] for c in cs]]
    ax.legend(handles=fam_handles, title="redox family", loc="upper left", fontsize=15,
              title_fontsize=15, framealpha=0.95)

    out = FIGDIR / "candidate_E_redox_landscape.png"
    fig.savefig(out); plt.close(fig)
    print(f"  -> {out}")


def main():
    cands = load_candidates()
    plot_redox_landscape(cands)


if __name__ == "__main__":
    main()
