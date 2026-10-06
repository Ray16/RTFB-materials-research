#!/usr/bin/env python
"""Where do the systematically identified D3TaLES candidates sit relative to the registered
candidates? Reads results/discovery/ (scripts/mining/identify_candidates.py). PNG only.

  discovery_capacity_vs_sa.png   nominal capacity vs grafted SA — both computed by the SAME
                                 code for pool and registered candidates (comparable axes)
  discovery_chemspace_pca.png    chemical-space map: PCA of Morgan fingerprints (r=2, 2048 bit)
                                 of the GRAFTED species, pool + registered candidates
  discovery_d3tales_lambda.png   D3TaLES-reported electron lambda of the PARENT (their level:
                                 LC-wHPBE/def2-SVP gas) vs nominal capacity — pre-filter view;
                                 registered candidates appear only where their precursor is in
                                 D3TaLES (our production lambda is a different level, not mixed in)

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


def fig_capacity_sa(pool, front, cur):
    fig = plt.figure(figsize=(17, 15))
    ax = fig.add_axes([0.09, 0.36, 0.88, 0.60])
    for f, g in pool.groupby("families"):
        ax.scatter(g.capacity_nominal_mAh_g, g.SA_grafted, s=46, color=FAM.get(f, "#999"),
                   alpha=0.45, edgecolors="none", zorder=2)
    fr = front[front.source_id.isin(pool.source_id)]
    ax.scatter(fr.capacity_nominal_mAh_g, fr.SA_grafted, s=150, facecolor="none",
               edgecolor="#D55E00", lw=2.0, zorder=4)
    _draw_current(ax, cur, "capacity_nominal_mAh_g", "SA_grafted")
    ax.set_xlabel("nominal specific capacity (mAh g$^{-1}$, grafted monomer)")
    ax.set_ylabel("synthetic accessibility of grafted species (Ertl; lower = easier)")
    grid_xy(ax)
    _legend(ax, pool, [plt.Line2D([], [], ls="", marker="o", ms=12, mfc="none", mec="#D55E00",
                                  mew=2, label="pre-filter front (per family)")])
    _key(fig, cur, y=0.27)
    p = OUT / "discovery_capacity_vs_sa.png"; fig.savefig(p); plt.close(fig); print(" ->", p)


def _fps(smiles):
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    X = np.zeros((len(smiles), 2048), dtype=np.float32)
    for i, s in enumerate(smiles):
        fp = gen.GetFingerprint(Chem.MolFromSmiles(s))
        arr = np.zeros((2048,), dtype=np.float32)
        DataStructs.ConvertToNumpyArray(fp, arr); X[i] = arr
    return X


def fig_chemspace(pool, front, cur):
    smi = list(pool.grafted_smiles) + list(cur.grafted_smiles)
    X = _fps(smi)
    Xc = X - X.mean(0)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    P = Xc @ Vt[:2].T
    ev = (S ** 2) / (S ** 2).sum()
    n = len(pool)
    pool = pool.assign(pc1=P[:n, 0], pc2=P[:n, 1]); cur = cur.assign(pc1=P[n:, 0], pc2=P[n:, 1])
    fig = plt.figure(figsize=(17, 15))
    ax = fig.add_axes([0.09, 0.36, 0.88, 0.60])
    for f, g in pool.groupby("families"):
        ax.scatter(g.pc1, g.pc2, s=46, color=FAM.get(f, "#999"), alpha=0.45, edgecolors="none")
    fr = pool[pool.source_id.isin(front.source_id)]
    ax.scatter(fr.pc1, fr.pc2, s=150, facecolor="none", edgecolor="#D55E00", lw=2.0, zorder=4)
    _draw_current(ax, cur, "pc1", "pc2")
    ax.set_xlabel(f"PC1 ({ev[0]*100:.0f}% of fingerprint variance)")
    ax.set_ylabel(f"PC2 ({ev[1]*100:.0f}%)")
    grid_xy(ax)
    _legend(ax, pool, [plt.Line2D([], [], ls="", marker="o", ms=12, mfc="none", mec="#D55E00",
                                  mew=2, label="pre-filter front (per family)")])
    _key(fig, cur, y=0.27)
    p = OUT / "discovery_chemspace_pca.png"; fig.savefig(p); plt.close(fig); print(" ->", p)


def fig_d3_lambda(pool, front, cur):
    g = pool[pool.d3_lambda_valid.astype(bool)]          # physical D3TaLES lambda only
    fig = plt.figure(figsize=(17, 15))
    ax = fig.add_axes([0.09, 0.36, 0.88, 0.60])
    for f, gg in g.groupby("families"):
        ax.scatter(gg.capacity_nominal_mAh_g, gg.d3_electron_reorg_eV, s=46,
                   color=FAM.get(f, "#999"), alpha=0.45, edgecolors="none")
    fr = g[g.source_id.isin(front.source_id)]
    ax.scatter(fr.capacity_nominal_mAh_g, fr.d3_electron_reorg_eV, s=150, facecolor="none",
               edgecolor="#D55E00", lw=2.0, zorder=4)
    # registered candidates whose un-grafted precursor is itself in D3TaLES
    rd = pd.read_csv(D / "rediscovery_check.csv")
    d3 = pd.read_csv(ROOT / "data/raw/validation/D3TaLES/d3tales_public.csv", low_memory=False,
                     usecols=["smiles", "electron_reorganization_energy"]).dropna()
    key = lambda s: (Chem.MolToInchiKey(Chem.MolFromSmiles(s))[:14]
                     if Chem.MolFromSmiles(str(s)) else None)
    d3k = {key(s): v for s, v in zip(d3.smiles, d3.electron_reorganization_energy)}
    hits = cur.merge(rd[["id", "precursor"]], on="id")
    hits["d3_lambda"] = [d3k.get(key(p)) for p in hits.precursor]
    hits = hits.dropna(subset=["d3_lambda"])
    _draw_current(ax, hits, "capacity_nominal_mAh_g", "d3_lambda")
    ax.set_xlabel("nominal specific capacity (mAh g$^{-1}$, grafted monomer)")
    ax.set_ylabel("D3TaLES electron λ of the parent (eV, their level)")
    ax.set_ylim(-0.05, 1.55)
    grid_xy(ax)
    _legend(ax, g, [plt.Line2D([], [], ls="", marker="o", ms=12, mfc="none", mec="#D55E00",
                               mew=2, label="pre-filter front (per family)")])
    hits = hits[hits.d3_lambda.between(0.0, 1.5)]
    fig.text(0.02, 0.30, f"only physical D3TaLES λ (0–1.5 eV): {g.source_id.nunique()} molecules. "
             "Registered candidates shown where their precursor is in D3TaLES with a physical "
             f"λ: {', '.join(map(str, hits.num)) or 'none'}", fontsize=18, va="top", wrap=True)
    _key(fig, cur, y=0.25)
    p = OUT / "discovery_d3tales_lambda.png"; fig.savefig(p); plt.close(fig); print(" ->", p)


def main():
    apply_style()
    OUT.mkdir(parents=True, exist_ok=True)
    pool, front, cur = _load()
    fig_capacity_sa(pool, front, cur)
    fig_chemspace(pool, front, cur)
    fig_d3_lambda(pool, front, cur)


if __name__ == "__main__":
    main()
