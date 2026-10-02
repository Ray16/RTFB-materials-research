#!/usr/bin/env python
"""Born vs molecular-cavity PCM outer-sphere lambda for OUR library molecules (MeCN).
Shows the PCM/Born ratio per couple, colored by category (grafted / bare _sa-parent / small),
revealing that Born UNDERESTIMATES lambda_out for grafted candidates (inert handle inflates the
SASA sphere while charge stays on the core). PNG only.

  python scripts/plotting/reorg/plot_lambda_outer_our_molecules.py
"""
from __future__ import annotations
import csv, glob
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[3]
import sys; sys.path.insert(0, str(ROOT / "scripts" / "plotting"))
import plot_style as ps  # noqa: E402

OUT = ROOT / "results/figures/reorg/lambda_outer_our_molecules_mecn.png"
CAT_COLOR = {"grafted (benzyl handle)": "#C0392B", "bare (_sa/parent)": "#2B6CB0",
             "small/validation": "#2E7D32"}


def category(i):
    if "_sa_" in i or i.endswith("_sa") or "parent" in i:
        return "bare (_sa/parent)"
    if i.startswith(("ferrocene", "methyl_viologen", "methylpyridinium")):
        return "small/validation"
    return "grafted (benzyl handle)"


STATES = {"neu", "ox", "ox1", "ox2", "red", "red1", "red2", "rad"}


def _label(i):
    """'ndi_ammonium_sa_red1_ox' -> 'ndi_ammonium_sa  (red1/ox)'."""
    parts = i.split("_")
    if len(parts) >= 3 and parts[-1] in STATES and parts[-2] in STATES:
        return f"{'_'.join(parts[:-2])}  ({parts[-2]}/{parts[-1]})"
    return i


def main():
    ps.apply_style()
    rows = []
    for f in sorted(glob.glob(str(ROOT / "results/reorg_anchors/our_chunks/chunk*_results.csv"))):
        rows += list(csv.DictReader(open(f)))
    rows = [r for r in rows if r.get("pcm_over_born") not in ("", "nan", None)]
    # de-dup by id (keep first)
    seen = {}
    for r in rows:
        seen.setdefault(r["id"], r)
    rows = list(seen.values())
    for r in rows:
        r["cat"] = category(r["id"])
        r["ratio"] = float(r["pcm_over_born"])
    rows.sort(key=lambda r: r["ratio"])

    labels = [_label(r["id"]) for r in rows]
    ratios = [r["ratio"] for r in rows]
    colors = [CAT_COLOR[r["cat"]] for r in rows]

    fig, ax = plt.subplots(figsize=(13, max(8, 0.52 * len(rows))))
    y = np.arange(len(rows))
    ax.barh(y, ratios, color=colors, height=0.7)
    ax.axvline(1.0, color="#333", lw=1.4, ls="--", label="Born = PCM")
    ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=15)
    ax.set_ylim(-0.7, len(rows) - 0.3)
    ax.set_xlabel("PCM / Born  (outer-sphere $\\lambda_o$ ratio, MeCN)")
    ax.set_xlim(0.95, max(ratios) * 1.07)
    ax.set_title("Born underestimates $\\lambda_o$ for grafted candidates\n"
                 "(molecular-cavity PCM / Born single-sphere, MeCN)", fontsize=20)
    for i, r in enumerate(rows):
        ax.text(r["ratio"] + 0.005, i, f"{r['ratio']:.2f}", va="center", fontsize=14, color="#333")
    ax.grid(True, axis="x", color=ps.C["grid"], lw=0.6, alpha=0.35); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    # legend by category
    from matplotlib.patches import Patch
    handles = [Patch(color=c, label=k) for k, c in CAT_COLOR.items()]
    handles.append(plt.Line2D([0], [0], color="#333", ls="--", label="Born = PCM"))
    ax.legend(handles=handles, frameon=False, fontsize=16, loc="lower right")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=300, bbox_inches="tight")   # PNG only
    import statistics as st
    print(f"wrote {OUT}  (n={len(rows)})")
    for c in CAT_COLOR:
        rr = [r["ratio"] for r in rows if r["cat"] == c]
        if rr:
            print(f"  {c:26s} n={len(rr):2d} mean PCM/Born={st.mean(rr):.2f}  range {min(rr):.2f}-{max(rr):.2f}")


if __name__ == "__main__":
    main()
