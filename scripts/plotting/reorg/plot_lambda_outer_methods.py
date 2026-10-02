#!/usr/bin/env python
"""Compare outer-sphere lambda methods (Born single-sphere vs molecular-cavity PCM) across a
shape-spanning set, in MeCN — so the agreement/disagreement is visible at a glance.

Left: grouped bars (Born vs PCM) per molecule, ordered by size (SASA radius), with the PCM/Born
ratio annotated. Right: Born-vs-PCM parity. PNG only.

  python scripts/plotting/reorg/plot_lambda_outer_methods.py
"""
from __future__ import annotations
import csv
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
import sys; sys.path.insert(0, str(ROOT / "scripts" / "plotting"))
import plot_style as ps  # noqa: E402

SRC = ROOT / "results/reorg_anchors/lambda_out_mecn_shapes_results.csv"
OUT = ROOT / "results/figures/reorg/lambda_outer_methods_mecn.png"
C_BORN, C_PCM = "#2B6CB0", "#DD6B20"


def main():
    ps.apply_style()
    rows = list(csv.DictReader(open(SRC)))
    rows.sort(key=lambda r: float(r["sasa_radius_A"]))
    ids = [r["id"] for r in rows]
    born = np.array([float(r["born_lambda_o_eV"]) for r in rows])
    pcm = np.array([float(r["pcm_lambda_o_eV"]) for r in rows])
    ratio = pcm / born

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(22, 7.8), gridspec_kw={"width_ratios": [1.9, 1], "wspace": 0.22})

    # --- grouped bars ---
    x = np.arange(len(ids)); w = 0.38
    ax.bar(x - w/2, born, w, color=C_BORN, label="Born single-sphere")
    ax.bar(x + w/2, pcm, w, color=C_PCM, label="molecular-cavity PCM")
    for i in range(len(ids)):
        top = max(born[i], pcm[i])
        ax.text(i, top + 0.02, f"×{ratio[i]:.2f}", ha="center", va="bottom", fontsize=15, color="#333")
    ax.set_xticks(x)
    ax.set_xticklabels([f"{ids[i]}\n({rows[i]['shape'].replace('-', chr(10))})"
                        for i in range(len(ids))], fontsize=15, rotation=0)
    ax.set_ylabel("outer-sphere $\\lambda_o$  (eV)")
    ax.set_title("Born vs molecular-cavity PCM (1-body)")
    ax.set_ylim(0, max(born.max(), pcm.max()) * 1.18)
    ax.grid(True, axis="y", color=ps.C["grid"], lw=0.6, alpha=0.35); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.legend(frameon=False, fontsize=16, loc="upper right")

    # --- parity ---
    lo = min(born.min(), pcm.min()) * 0.9
    hi = max(born.max(), pcm.max()) * 1.05
    ax2.plot([lo, hi], [lo, hi], color="#333", lw=1.4, zorder=1, label="y = x")
    ax2.scatter(born, pcm, s=90, c=C_PCM, edgecolors="white", linewidths=1.0, zorder=3)
    # hand-placed offsets (points): phenothiazine and anthraquinone sit ~15 meV apart
    offs = {"phenothiazine": (-10, 12, "right"), "anthraquinone": (10, -16, "left"),
            "benzoquinone": (-10, 0, "right")}
    for i in range(len(ids)):
        dx, dy, ha = offs.get(ids[i], (10, -4, "left"))
        ax2.annotate(ids[i], (born[i], pcm[i]), textcoords="offset points", xytext=(dx, dy),
                     fontsize=15, color="#1A202C", ha=ha, va="center")
    mad = np.abs(pcm - born).mean()
    ax2.set_xlim(lo, hi); ax2.set_ylim(lo, hi); ax2.set_aspect("equal")
    ax2.set_xlabel("Born $\\lambda_o$ (eV)"); ax2.set_ylabel("PCM $\\lambda_o$ (eV)")
    ax2.set_title(f"parity  (MAD = {mad*1000:.0f} meV)")
    ax2.grid(True, color=ps.C["grid"], lw=0.6, alpha=0.3); ax2.set_axisbelow(True)
    for s in ("top", "right"):
        ax2.spines[s].set_visible(False)
    ax2.legend(frameon=False, fontsize=16, loc="upper left")

    fig.suptitle("Outer-sphere reorganization energy in MeCN: continuum method comparison",
                 fontsize=20, y=1.0)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=300, bbox_inches="tight")   # PNG only
    print(f"wrote {OUT}  (n={len(ids)}, PCM-vs-Born MAD={mad*1000:.0f} meV, ratio {ratio.min():.2f}-{ratio.max():.2f})")


if __name__ == "__main__":
    main()
