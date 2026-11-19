#!/usr/bin/env python
"""Figure 3 (manuscript, JCIM): the reliability-gating funnel.

NEW generator (no prior script). A four-tier funnel showing how the identity
registry narrows to the scored candidate pool through explicit QC gates:

    registered graphs -> state-complete -> QC-clean -> scored candidates

All counts are computed from committed readiness records (CPU only; no QM):
  * manuscript/evidence/2026-11-18_cpu_checks/readiness_registry.csv
      one row per declared redox state of every registered molecule
      (redox.quality.readiness over the whole registry).
      registered molecules = distinct ids; state-complete = states without
      missing_geometry_or_result.
  * manuscript/evidence/2026-11-18_cpu_checks/readiness_pool91.csv
      the same check restricted to the 91 scored candidates.
      QC-clean = states with no issue other than
      missing_validated_native_thermal. That certificate is excluded because the
      pool is ranked on the xTB thermal tier, not the native tier (the
      manuscript Methods and SI S6 state this exclusion).
  * scored candidates = distinct ids in readiness_pool91.csv (91).

The funnel spans THREE grains, annotated on the left: MOLECULE grain
(registered molecules), STATE grain (the two middle tiers, manifest redox states), and CANDIDATE grain
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
EVID = ROOT / "manuscript" / "evidence" / "2026-11-18_cpu_checks"
REGISTRY = EVID / "readiness_registry.csv"
POOL = EVID / "readiness_pool91.csv"
EXCLUDED_GATE = "missing_validated_native_thermal"   # native tier; pool uses xTB tier


def funnel_counts() -> dict:
    reg = pd.read_csv(REGISTRY).fillna("")
    pool = pd.read_csv(POOL).fillna("")
    qc_clean = pool["issues"].isin(["", EXCLUDED_GATE])
    return {
        "molecules": reg["id"].nunique(),
        "declared": len(reg),
        "complete": int((~reg["issues"].str.contains("missing_geometry_or_result")).sum()),
        "pool_states": len(pool),
        "qc_clean": int(qc_clean.sum()),
        "flagged": int((~qc_clean).sum()),
        "scored": pool["id"].nunique(),
    }


def main():
    apply_style()
    k = funnel_counts()
    N_DECLARED, N_STATE_COMPLETE = k["declared"], k["complete"]
    N_POOL_STATES, N_QC_CLEAN = k["pool_states"], k["qc_clean"]

    # (label, count-text, gate description [wrapped], width, color, placeholder?)
    tiers = [
        ("Registered molecules",
         str(k["molecules"]), f"candidate registry\n({N_DECLARED} declared\nredox states)",
         1.00, C["measured"], False),
        ("State-complete",
         str(N_STATE_COMPLETE),
         f"converged geometry + SP\n({N_STATE_COMPLETE} of {N_DECLARED}\ndeclared states)",
         0.84, C["measured"], False),
        ("QC-clean",
         str(N_QC_CLEAN), f"all ranking gates pass\n({N_QC_CLEAN} of {N_POOL_STATES} states\nof the {k['scored']} candidates)",
         0.70, C["measured"], False),
        ("Scored candidates",
         str(k["scored"]), f"ranked on Pareto fronts\n({k['flagged']} carry a flagged state;\nneither on a front)",
         0.62, C["good"], False),
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
    ax.annotate("molecule grain", xy=(-0.72, centers[0][0]), ha="center",
                va="center", fontsize=18, color="#666666", rotation=90)
    ax.annotate("state grain", xy=(-0.72, 0.5 * (centers[1][0] + centers[2][0])),
                ha="center", va="center", fontsize=18, color="#666666", rotation=90)
    ax.annotate("candidate grain", xy=(-0.72, centers[3][0]), ha="center",
                va="center", fontsize=18, color="#666666", rotation=90)

    ax.set_title("Reliability-gating funnel: registry → scored pool",
                 fontsize=22, pad=14)
    sub = ("State counts from redox.quality.readiness (2026-11-18). Bar widths are "
           "for legibility\nand do not scale with the counts; the tiers count "
           "molecules, states, and candidates.")
    fig.text(0.5, 0.012, sub, ha="center", va="bottom", fontsize=18,
             color="#555555")
    ax.set_ylim(0.02, 1.02)
    ax.axis("off")
    fig.subplots_adjust(left=0.02, right=0.98, top=0.92, bottom=0.16)
    fig.savefig(FIGDIR / "fig3_qc_funnel.png", dpi=300)
    plt.close(fig)
    print("funnel counts:", k)
    print("  ->", FIGDIR / "fig3_qc_funnel.png")


if __name__ == "__main__":
    main()
