#!/usr/bin/env python
"""Where do the systematically identified D3TaLES candidates sit relative to the registered
candidates? Reads results/discovery/ (scripts/mining/identify_candidates.py). PNG only.

  discovery_capacity_vs_sa.png   nominal capacity vs grafted SA — both computed by the SAME
                                 code for pool and registered candidates (comparable axes)

  python scripts/plotting/candidates/plot_discovery.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_style import apply_style, grid_xy  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np               # noqa: E402
import pandas as pd              # noqa: E402
from rdkit import Chem, DataStructs, RDLogger  # noqa: E402
from rdkit.Chem import rdFingerprintGenerator  # noqa: E402
from verified_overlay import load as load_overlay,finite

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).resolve().parents[3]
D = ROOT / "results" / "discovery"
OUT = ROOT / "results" / "figures" / "candidates"
FAM = {"p-quinone": "#0072B2", "o-quinone": "#56B4E9", "bis-imide": "#009E73",
       "viologen": "#CC79A7", "phenol": "#999999"}
FAM_LABEL = {"p-quinone": "p-quinone", "o-quinone": "o-quinone", "bis-imide": "aromatic bis-imide",
             "viologen": "viologen", "": "1e⁻ (phenol)"}


def _load():
    pool = pd.read_csv(D / "d3tales_pool.csv")
    pool = pool[(pool.status == "ok") & ~pool.parent_in_library.astype(bool)].copy()
    front = pd.read_csv(D / "d3tales_prefilter_front.csv")
    cur = pd.read_csv(D / "current_candidates.csv")
    cur = cur[cur.in_scorecard].reset_index(drop=True)
    cur["num"] = np.arange(1, len(cur) + 1)
    return pool, front, cur


def _key(fig, cur, y=0.0):
    """Numbered key with the registered candidates' real names (config `name`), two columns,
    placed under the axes (18 pt like every other text)."""
    lines = [f"{n}  {nm}" for n, nm in zip(cur.num, cur.name)]
    half = (len(lines) + 1) // 2
    for col, chunk in enumerate((lines[:half], lines[half:])):
        fig.text(0.02 + 0.5 * col, y, "\n".join(chunk), va="top", ha="left", fontsize=18)


def _draw_current(ax, cur, x, y):
    """Star at the TRUE position; the number label is moved off overlaps by adjustText and
    joined to its star by a leader line (positions are never altered, only labels)."""
    from adjustText import adjust_text
    texts = []
    for _, r in cur.iterrows():
        c = FAM.get(r.families, FAM["phenol"])
        ax.scatter(r[x], r[y], s=420, marker="*", facecolor=c, edgecolor="black", lw=1.2,
                   zorder=6)
        texts.append(ax.text(r[x], r[y], str(r.num), fontsize=18, fontweight="bold",
                             color=c, zorder=7))
    adjust_text(texts, ax=ax, expand=(1.6, 1.8), force_text=(0.6, 0.8),
                arrowprops=dict(arrowstyle="-", color="#555555", lw=1.2))


def _legend(ax, pool, extra=()):
    def n_mol(f):
        return int(pool[pool.families == f].source_id.nunique())
    hs = [plt.Line2D([], [], ls="", marker="o", ms=10, color=FAM[f], alpha=0.6,
                     label=f"D3TaLES {FAM_LABEL[f]}: {n_mol(f)} molecules")
          for f in ("p-quinone", "o-quinone") if (pool.families == f).any()]
    hs += list(extra)
    hs.append(plt.Line2D([], [], ls="", marker="*", ms=20, mfc="#888888", mec="black",
                         label="registered candidate (number = key)"))
    ax.legend(handles=hs, loc="upper right", fontsize=18)


def fig_capacity_sa(pool, front, cur, feed=None):
    fig = plt.figure(figsize=(17, 15))
    ax = fig.add_axes([0.09, 0.36, 0.88, 0.60])
    for f, g in pool.groupby("families"):
        ax.scatter(g.capacity_nominal_mAh_g, g.SA_grafted, s=46, color=FAM.get(f, "#999"),
                   alpha=0.45, edgecolors="none", zorder=2)
    fr = front[front.source_id.isin(pool.source_id)]
    ax.scatter(fr.capacity_nominal_mAh_g, fr.SA_grafted, s=150, facecolor="none",
               edgecolor="#D55E00", lw=2.0, zorder=4)
    _draw_current(ax, cur, "capacity_nominal_mAh_g", "SA_grafted")
    payload,expansion,admitted=load_overlay(feed)
    expansion=[r for r in expansion if r.get('descriptors_status')=='verified exact-grafted graph/arithmetic'
               and finite(r.get('capacity_nominal_mAh_g')) and finite(r.get('SA_grafted'))]
    if expansion:
        ax.scatter([r['capacity_nominal_mAh_g'] for r in expansion],
                   [r['SA_grafted'] for r in expansion],s=180,marker='D',facecolor='#D55E00',
                   edgecolor='black',linewidth=.9,zorder=5,alpha=.8)
    ax.set_xlabel("nominal specific capacity (mAh g$^{-1}$, grafted monomer)")
    ax.set_ylabel("synthetic accessibility of grafted species (Ertl; lower = easier)")
    grid_xy(ax)
    _legend(ax, pool, [plt.Line2D([], [], ls="", marker="o", ms=12, mfc="none", mec="#D55E00",
                                  mew=2, label="pre-filter front (per family)"),
                      plt.Line2D([],[],ls='',marker='D',ms=11,mfc='#D55E00',mec='black',
                                 label=f'new expansion: {len(expansion)} exact grafts')])
    _key(fig, cur, y=0.235)
    fig.text(.02,.275,'Nominal capacity and SA are graph-based descriptors; QM validation and synthetic feasibility remain separate.',fontsize=16)
    p = OUT / "discovery_capacity_vs_sa.png"; fig.savefig(p);fig.savefig(p.with_suffix('.pdf'))
    box=ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted()).expanded(1.03,1.06)
    fig.savefig(OUT/'discovery_capacity_vs_sa_slide.png',bbox_inches=box,dpi=220)
    plt.close(fig); print(" ->", p)


def main(feed=None):
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    pool, front, cur = _load()
    fig_capacity_sa(pool, front, cur,feed)


if __name__ == "__main__":
    main()
