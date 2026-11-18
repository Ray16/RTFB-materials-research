#!/usr/bin/env python
"""Figure 3 (manuscript, JCIM): the reliability-gating funnel.

NEW generator (no prior script). A four-tier funnel showing how the identity
registry narrows to the scored candidate pool through explicit QC gates:

    registered graphs -> state-complete -> QC-clean -> scored candidates

Honesty about the numbers (CPU-only; no QM is run here):
  * REAL, verified:
      - registered graphs = 912 (identity registry; docs/PUBLICATION_READINESS.md)
      - state-complete = 450 of 456 declared states (pool-wide
        redox.quality.readiness run, 2026-11-17; 450 have converged geometry +
        active SP, 6 missing). This check is path-independent (see below).
      - scored candidates = 91 (manuscript pool), with the maintained
        results/validation/full_pool_scored_legacyG_MeCN_TBAPF6.csv currently
        carrying 92 rows (one extra phenazine catholyte from the live expansion).
  * PLACEHOLDER (one tier, flagged on the figure and in the caption): QC-clean.
    Its final gate needs validated native-thermal free energies, whose
    certificates are path-locked to the original /nfs/.../0_redox checkout; from
    this archive-symlinked tree readiness reports 351 spurious
    missing_validated_native_thermal, so the naive "103 ready" is a relocation
    artifact and is NOT shown. A faithful count needs
    redox.quality.relocated_native (or a checkout at the original path).

The funnel spans THREE grains, annotated on the left: GRAPH grain (identity
registry counts molecular graphs, incl. per-state and held/uncomputed graphs),
STATE grain (the two middle tiers, manifest redox states), and CANDIDATE grain
(scored motifs). Bar widths are a fixed decreasing schedule for legibility and do
NOT encode counts proportionally.

  python scripts/plotting/manuscript/fig3_qc_funnel.py
    -> manuscript/figures/fig3_qc_funnel.png (+ .caption.md)
"""
from __future__ import annotations
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_style import apply_style, C  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Polygon  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
FIGDIR = ROOT / "manuscript" / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)
CSV = ROOT / "results" / "validation" / "full_pool_scored_legacyG_MeCN_TBAPF6.csv"


def scored_count() -> int:
    try:
        return int(len(pd.read_csv(CSV)))
    except Exception:
        return -1


def main():
    apply_style()
    n_csv = scored_count()

    # State-complete is REAL and path-independent: from a pool-wide
    # `redox.quality.readiness` run (manifest grain, 2026-11-17), 456 redox
    # states are declared and 450 have a converged geometry + active single-point
    # result (only 6 are missing_geometry_or_result). This check reads result.json
    # + opt.xyz + composition/connectivity/spin; it does NOT touch the native-
    # thermal certificate, so it is unaffected by the archive relocation.
    #
    # QC-clean stays a flagged PLACEHOLDER: the final QC gate includes validated
    # native-thermal free energies, whose certificates are path-locked to the
    # original /nfs/.../0_redox checkout. Read from this archive-symlinked tree
    # the readiness run reports 351 spurious `missing_validated_native_thermal`,
    # so a faithful QC-clean count needs redox.quality.relocated_native (or a
    # checkout at the original path); the 103 "ready" here is a relocation
    # artifact and is deliberately NOT shown.
    N_DECLARED, N_STATE_COMPLETE = 456, 450

    # (label, count-text, gate description [wrapped], width, color, placeholder?)
    tiers = [
        ("Registered graphs",
         "912", "identity registry\n(all graphs, incl.\nheld / uncomputed)",
         1.00, C["measured"], False),
        ("State-complete",
         str(N_STATE_COMPLETE),
         f"converged geometry + SP\n({N_STATE_COMPLETE} of {N_DECLARED}\ndeclared states)",
         0.80, C["measured"], False),
        ("QC-clean",
         "placeholder", "charge/spin, ⟨S²⟩, geometry,\ncoverage, cross-points pass",
         0.62, "#9AA3AD", True),
        ("Scored candidates",
         "91", "ranked on σ-aware\nPareto fronts",
         0.46, C["good"], False),
    ]

    fig, ax = plt.subplots(figsize=(16.5, 10.5))
    ax.set_xlim(-1.12, 1.45)
    n = len(tiers)
    top, h, gap = 1.0, 0.9 / n * 0.82, 0.9 / n * 0.18
    y = top
    centers = []
    for i, (label, count, desc, w, color, placeholder) in enumerate(tiers):
        wnext = tiers[i + 1][3] if i + 1 < n else w * 0.86
        yb, yt = y - h, y
        xl, xr = -w / 2, w / 2
        xln, xrn = -wnext / 2, wnext / 2
        # trapezoid funnel segment
        poly = Polygon([(xl, yt), (xr, yt), (xrn, yb), (xln, yb)],
                       closed=True, facecolor=color,
                       edgecolor=(C["warn"] if placeholder else "#333333"),
                       lw=(2.4 if placeholder else 1.6),
                       hatch=("///" if placeholder else None), alpha=0.92, zorder=2)
        ax.add_patch(poly)
        ycen = 0.5 * (yb + yt)
        centers.append((ycen, yb))
        # in-bar count (big) + label
        ctxt = count if not placeholder else "placeholder"
        ax.text(0, ycen + 0.028, ctxt, ha="center", va="center",
                fontsize=(30 if not placeholder else 18),
                fontweight="bold",
                color=("white" if not placeholder else "#7a1414"), zorder=3)
        ax.text(0, ycen - 0.040, label, ha="center", va="center", fontsize=18,
                fontweight="bold",
                color=("white" if not placeholder else "#333333"), zorder=3)
        # gate description to the right, via a leader line
        ax.annotate(desc, xy=(xr if w >= wnext else xrn, ycen),
                    xytext=(0.62, ycen), textcoords="data", ha="left",
                    va="center", fontsize=18, color="#333333",
                    arrowprops=dict(arrowstyle="-", lw=1.1, color="#AAAAAA"))
        y = yb - gap

    # "held, not dropped" side note
    ax.text(-1.08, centers[2][0], "failures are\nHELD,\nnot dropped",
            ha="left", va="center", fontsize=18, color=C["warn"],
            fontweight="bold")

    # grain annotation (rotated, far left) — the funnel spans three grains
    ax.annotate("graph grain", xy=(-0.72, centers[0][0]), ha="center",
                va="center", fontsize=18, color="#666666", rotation=90)
    ax.annotate("state grain", xy=(-0.72, 0.5 * (centers[1][0] + centers[2][0])),
                ha="center", va="center", fontsize=18, color="#666666", rotation=90)
    ax.annotate("candidate grain", xy=(-0.72, centers[3][0]), ha="center",
                va="center", fontsize=18, color="#666666", rotation=90)

    ax.set_title("Reliability-gating funnel: registry → scored pool",
                 fontsize=22, pad=14)
    sub = ("Real: 912 registered (identity registry), 450/456 state-complete "
           "(readiness run),\nand 91 scored (maintained CSV currently "
           f"{n_csv} rows). The hatched QC-clean tier is a\nplaceholder: its "
           "native-thermal gate is path-locked to the original checkout and "
           "needs\nredox.quality.relocated_native for a faithful count.")
    fig.text(0.5, 0.012, sub, ha="center", va="bottom", fontsize=18,
             color="#555555")
    ax.set_ylim(0.02, 1.02)
    ax.axis("off")
    fig.subplots_adjust(left=0.02, right=0.98, top=0.92, bottom=0.16)
    fig.savefig(FIGDIR / "fig3_qc_funnel.png", dpi=300)
    plt.close(fig)
    print("scored CSV rows:", n_csv)
    print("  ->", FIGDIR / "fig3_qc_funnel.png")


if __name__ == "__main__":
    main()
