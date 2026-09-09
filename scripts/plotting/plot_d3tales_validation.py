#!/usr/bin/env python
"""Validation figure: we INDEPENDENTLY reproduce D3TaLES's SA scores and reorganization
energies for a subset of D3TaLES structures (with viologen / quinone / imide highlighted).

  a) SA score  : our RDKit sascorer(SMILES) vs D3TaLES `sa_score`.
  b) Reorg     : our Nelsen 4-point sum (relaxation_gs->ion + relaxation_ion->gs) vs D3TaLES
                 reported hole/electron_reorganization_energy.

Both are parity plots against y = x, annotated with n / R^2 / MAE. Target redox-flow families
(viologen, quinone, imide) are colored; everything else is a light-grey background cloud.

  PYTHONPATH=src python scripts/plotting/plot_d3tales_validation.py
"""
from __future__ import annotations
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
CSV = ROOT / "data" / "raw" / "validation" / "D3TaLES" / "d3tales_public.csv"
OUT = ROOT / "results" / "figures" / "validation" / "d3tales_validation.png"
OUT.parent.mkdir(parents=True, exist_ok=True)

FAMS = {
    "viologen": ["[n+]1ccc(cc1)-c1cc[n+]cc1", "c1cc[n]cc1-c1cc[n]cc1"],
    "quinone":  ["O=C1C=CC(=O)C=C1", "O=C1c2ccccc2C(=O)c2ccccc21", "O=C1C=CC(=O)c2ccccc21"],
    "imide":    ["O=C[NX3]C=O"],
}
COLOR = {"viologen": "#CC79A7", "quinone": "#0072B2", "imide": "#009E73"}


def _style(plt):
    plt.rcParams.update({
        "font.size": 18, "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "axes.linewidth": 1.4, "axes.edgecolor": "#2b2b2b",
        "xtick.major.width": 1.3, "ytick.major.width": 1.3,
        "xtick.major.size": 6, "ytick.major.size": 6,
    })


def classify(df):
    from rdkit import Chem
    from rdkit import RDLogger; RDLogger.DisableLog("rdApp.*")
    pat = {k: [Chem.MolFromSmarts(p) for p in v if Chem.MolFromSmarts(p)] for k, v in FAMS.items()}
    fams, mols = [], []
    for smi in df["smiles"]:
        m = Chem.MolFromSmiles(smi) if isinstance(smi, str) else None
        mols.append(m)
        lab = "other"
        if m is not None:
            for k, ps in pat.items():
                if any(m.HasSubstructMatch(p) for p in ps):
                    lab = k; break
        fams.append(lab)
    return np.array(fams), mols


def _stats(x, y):
    x, y = np.asarray(x), np.asarray(y)
    r2 = np.corrcoef(x, y)[0, 1] ** 2
    return len(x), r2, np.abs(x - y).mean()


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from rdkit.Chem import RDConfig
    sys.path.append(f"{RDConfig.RDContribDir}/SA_Score"); import sascorer
    _style(plt)

    df = pd.read_csv(CSV, low_memory=False).reset_index(drop=True)
    fam, mols = classify(df)
    df["fam"] = fam
    print("family members:", {k: int((fam == k).sum()) for k in FAMS})

    # background sample (other) + all family members, for a legible cloud
    rng = np.random.default_rng(42)
    other_idx = rng.choice(np.where(fam == "other")[0], size=900, replace=False)
    fam_idx = np.where(fam != "other")[0]
    idx = np.concatenate([other_idx, fam_idx])

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(15.5, 7.6))
    fig.subplots_adjust(left=0.075, right=0.985, top=0.9, bottom=0.14, wspace=0.24)

    # ---- (a) SA score parity ----
    sa_ref, sa_rec, sa_fam = [], [], []
    for i in idx:
        if mols[i] is None or pd.isna(df.at[i, "sa_score"]):
            continue
        try:
            sa_rec.append(sascorer.calculateScore(mols[i]))
            sa_ref.append(float(df.at[i, "sa_score"])); sa_fam.append(fam[i])
        except Exception:
            pass
    sa_ref, sa_rec, sa_fam = np.array(sa_ref), np.array(sa_rec), np.array(sa_fam)
    lim = (min(sa_ref.min(), sa_rec.min()) - 0.2, max(sa_ref.max(), sa_rec.max()) + 0.2)
    axa.plot(lim, lim, "-", color="#888", lw=1.5, zorder=1)
    axa.scatter(sa_ref[sa_fam == "other"], sa_rec[sa_fam == "other"], s=22, color="#c8ccd2",
                edgecolor="none", zorder=2, label="other D3TaLES")
    for f in ["quinone", "imide", "viologen"]:
        mask = sa_fam == f
        axa.scatter(sa_ref[mask], sa_rec[mask], s=70, color=COLOR[f], edgecolor="white",
                    linewidth=0.8, zorder=3, label=f)
    n, r2, mae = _stats(sa_ref, sa_rec)
    axa.text(0.05, 0.95, f"n = {n}\nR² = {r2:.4f}\nMAE = {mae:.3f}", transform=axa.transAxes,
             va="top", ha="left", fontsize=17,
             bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#cccccc"))
    axa.set_xlim(*lim); axa.set_ylim(*lim)
    axa.set_xlabel("D3TaLES SA score", fontsize=19)
    axa.set_ylabel("Our SA score (RDKit)", fontsize=19)
    axa.set_title("Synthetic accessibility", fontsize=20, fontweight="bold", pad=12)
    axa.text(0.03, 1.06, "a", transform=axa.transAxes, fontsize=24, fontweight="bold")
    axa.spines[["top", "right"]].set_visible(False); axa.tick_params(labelsize=17)

    # ---- (b) reorganization-energy parity (hole + electron) ----
    comps = [("relaxation_groundState_cation1", "relaxation_cation1_groundState",
              "hole_reorganization_energy"),
             ("relaxation_groundState_anion1", "relaxation_anion1_groundState",
              "electron_reorganization_energy")]
    rr_ref, rr_rec, rr_fam = [], [], []
    sub = df.iloc[idx]
    for c1, c2, rep in comps:
        s = sub[[c1, c2, rep, "fam"]].copy()
        s[[c1, c2, rep]] = s[[c1, c2, rep]].apply(pd.to_numeric, errors="coerce")
        s = s.dropna(subset=[c1, c2, rep])
        s = s[(s[rep] >= 0) & (s[rep] < 2.5)]        # drop the rare anion-column outliers
        rr_ref += list(s[rep]); rr_rec += list(s[c1] + s[c2]); rr_fam += list(s["fam"])
    rr_ref, rr_rec, rr_fam = np.array(rr_ref), np.array(rr_rec), np.array(rr_fam)
    lim2 = (-0.03, max(rr_ref.max(), rr_rec.max()) * 1.05)
    axb.plot(lim2, lim2, "-", color="#888", lw=1.5, zorder=1)
    axb.scatter(rr_ref[rr_fam == "other"], rr_rec[rr_fam == "other"], s=22, color="#c8ccd2",
                edgecolor="none", zorder=2, label="other D3TaLES")
    for f in ["quinone", "imide", "viologen"]:
        mask = rr_fam == f
        axb.scatter(rr_ref[mask], rr_rec[mask], s=70, color=COLOR[f], edgecolor="white",
                    linewidth=0.8, zorder=3, label=f)
    n, r2, mae = _stats(rr_ref, rr_rec)
    axb.text(0.05, 0.95, f"n = {n}\nR² = {r2:.5f}\nMAE = {mae*1000:.2f} meV", transform=axb.transAxes,
             va="top", ha="left", fontsize=17,
             bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#cccccc"))
    axb.set_xlim(*lim2); axb.set_ylim(*lim2)
    axb.set_xlabel("D3TaLES reorganization energy (eV)", fontsize=19)
    axb.set_ylabel("Our 4-point λ (eV)", fontsize=19)
    axb.set_title("Reorganization energy", fontsize=20, fontweight="bold", pad=12)
    axb.text(0.03, 1.06, "b", transform=axb.transAxes, fontsize=24, fontweight="bold")
    axb.spines[["top", "right"]].set_visible(False); axb.tick_params(labelsize=17)

    axb.legend(loc="lower right", fontsize=15, frameon=False, markerscale=1.1)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=300, bbox_inches="tight", facecolor="white")
    print("wrote", OUT)


if __name__ == "__main__":
    main()
