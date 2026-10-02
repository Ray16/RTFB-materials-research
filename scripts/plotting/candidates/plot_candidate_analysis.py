#!/usr/bin/env python
"""Candidate payoff figures — the plots that actually rank and present the six monomers.

Fills the gaps in results/figures/candidates/:
  1. pareto_front.png    sigma-aware Pareto front in lambda-vs-SA (both minimized), no
                         scalarized figure-of-merit (project preference); raw axes,
                         family-coloured, capacity as bubble size, non-dominated set marked.
  2. redox_landscape.png computed E deg vs Fc/Fc+ per candidate (sorted), with sigma_E error
                         bars and specific capacity annotated — the "where do they sit" plot.

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


def pareto_min2d(xs, ys):
    """Boolean mask of the non-dominated set for (min x, min y)."""
    xs = np.asarray(xs, float); ys = np.asarray(ys, float)
    keep = np.ones(len(xs), bool)
    for i in range(len(xs)):
        for j in range(len(xs)):
            if j == i:
                continue
            # j dominates i if j is <= on both and < on at least one
            if xs[j] <= xs[i] and ys[j] <= ys[i] and (xs[j] < xs[i] or ys[j] < ys[i]):
                keep[i] = False
                break
    return keep


# --------------------------------------------------------------- Pareto figure
def plot_pareto(cands):
    apply_style()
    fig, ax = plt.subplots(figsize=(11, 8.5))

    # A candidate whose every lambda was QC-flagged now has lambda_i_eV empty (see
    # redox.screening.scorecard): it has no trustworthy kinetics axis, so it cannot be placed on a
    # lambda-vs-SA plane at all. Plotting it at some fallback value would invent data.
    # Drop it from the plane and say so on the figure.
    plotted = [c for c in cands if c["lam"] is not None and c["SA"] is not None]
    dropped = [c for c in cands if c not in plotted]
    if not plotted:
        print("  !! no candidate has a QC-clean lambda — Pareto plane not drawn")
        plt.close(fig)
        return

    lam = np.array([c["lam"] for c in plotted], float)
    sa = np.array([c["SA"] for c in plotted], float)
    slam = np.array([c["slam"] or 0.0 for c in plotted], float)
    cap = np.array([c["cap"] for c in plotted], float)
    fams = [c["family"] for c in plotted]
    batches = [c["batch"] or "starting" for c in plotted]

    smin, smax = cap.min(), cap.max()
    sizes = 300 + 850 * (cap - smin) / (smax - smin + 1e-9)

    front = pareto_min2d(lam, sa)

    ax.set_xlim(lam.min() - 0.06, lam.max() + 0.10)
    ax.set_ylim(sa.min() - 0.13, sa.max() + 0.16)
    ax.axhspan(sa.min() - 0.13, np.median(sa), xmin=0, xmax=1, color="#F3F6F2", zorder=0)

    fx, fy = lam[front], sa[front]
    order = np.argsort(fx)
    ax.plot(fx[order], fy[order], color="#555555", lw=2.0, ls="--", zorder=2,
            label="Pareto front (non-dominated)")

    # NOTE on lambda error bars: the ABSOLUTE lambda uncertainty (~0.21 eV vs D3TaLES) is
    # LARGER than the whole candidate spread (~0.09 eV), so a per-point bar would span the
    # entire axis and mislead. It is stated in text instead; the axis resolves the RELATIVE
    # ordering at fixed level of theory (where systematic method error cancels).
    texts = []
    for i, (c, x, y, s_, e) in enumerate(zip(plotted, lam, sa, sizes, slam)):
        col = FAM.get(c["family"], "#777777")
        b = c["batch"] or "starting"
        on = front[i]
        new_cohort = (b == "merrifield_multi")
        # family -> fill; batch -> marker shape + accent ring; Pareto front -> dark ring.
        # A new-cohort point on the front keeps the dark front ring and is identified by
        # its diamond, so the two encodings never collide.
        ax.scatter(x, y, s=s_, color=col, marker=BATCH_MARKER.get(b, "o"),
                   alpha=0.85 if on else 0.45,
                   edgecolors=("#222222" if on else (BATCH_ACCENT if new_cohort else "none")),
                   linewidths=(2.4 if on else (2.6 if new_cohort else 0)), zorder=4)
        # label beside the bubble, offset by its radius (adjustText ignores marker size and
        # parked labels on top of the big bubbles); flip to the left near the right edge
        r_pt = (s_ ** 0.5) / 2 + 7
        right = x < lam.min() + 0.75 * (lam.max() - lam.min())
        ax.annotate(SHORT.get(c["id"], c["id"]), (x, y), xytext=((r_pt if right else -r_pt), 0),
                    textcoords="offset points", ha=("left" if right else "right"),
                    va="center", fontsize=15, color="#222222", zorder=6)

    ax.set_xlabel("inner-sphere reorganization energy  $\\lambda_i$  (eV)   —  lower is faster")
    ax.set_ylabel("synthetic accessibility (SA)   —  lower is easier")
    ax.set_title("Candidate trade-off: reorganization vs. synthesizability", pad=14)
    grid_xy(ax)

    from matplotlib.lines import Line2D
    fam_handles = [Line2D([0], [0], marker="o", ls="", ms=13, mfc=FAM[f], mec="none",
                          label=FAM_LABEL.get(f, f)) for f in FAM if f in fams]
    batch_handles = [Line2D([0], [0], marker=BATCH_MARKER[b], ls="", ms=13, mfc="#BBBBBB",
                            mec=(BATCH_ACCENT if b == "merrifield_multi" else "none"),
                            mew=2.6, label=BATCH_LABEL[b])
                     for b in ("starting", "merrifield_multi") if b in batches]
    front_handle = Line2D([0], [0], marker="o", ls="", ms=15, mfc="#BBBBBB",
                          mec="#222222", mew=2.2, label="on Pareto front")
    line_handle = Line2D([0], [0], color="#555555", lw=2.0, ls="--", label="front")
    # Three legends + a caption; keep them in separate corners so nothing overlaps
    # (house rule: no overlapping elements).
    leg1 = ax.legend(handles=fam_handles, title="redox family", loc="lower left",
                     frameon=True, fontsize=14, title_fontsize=14, framealpha=0.95)
    ax.add_artist(leg1)
    if batch_handles:
        leg2 = ax.legend(handles=batch_handles, title="cohort", loc="upper left",
                         frameon=True, fontsize=14, title_fontsize=14, framealpha=0.95)
        ax.add_artist(leg2)
    ax.legend(handles=[front_handle, line_handle], loc="upper right", fontsize=13,
              frameon=True, framealpha=0.95)

    # Caption lives BELOW the axes (figure coords), never on top of the data or a legend.
    note = ("bubble size $\\propto$ specific capacity   ·   "
            "absolute $\\lambda_i$ uncertainty $\\pm$"
            f"{slam[0]:.2f} eV $>$ spread $\\rightarrow$ ranking is relative only")
    if dropped:
        note += ("\nexcluded, no QC-clean $\\lambda_i$: "
                 + ", ".join(SHORT.get(c["id"], c["id"]) for c in dropped))
    fig.text(0.5, 0.0, note, ha="center", va="top", fontsize=14,
             color="#666666", style="italic")

    out = FIGDIR / "candidate_pareto_lambda_SA.png"
    fig.savefig(out); plt.close(fig)
    print(f"  -> {out}   ({len(plotted)} plotted, {len(dropped)} excluded for lambda QC)")


# --------------------------------------------------------------- redox landscape
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
    plot_pareto(cands)
    plot_redox_landscape(cands)


if __name__ == "__main__":
    main()
