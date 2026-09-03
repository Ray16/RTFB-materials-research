#!/usr/bin/env python
"""Presentation figures for the ΔG_solv experimental-validation (Ray's Sept-2 action item).

Reads results/solvation_validation.csv (validate_solvation.py --report) and makes 300-dpi,
house-style slide figures:

  solvation_parity.png      2-panel parity: neutrals (zoom) | ions colored by CHARGE CHARACTER.
  solvation_mae_summary.png MAE by category, with SMD accuracy targets.

The ion panel groups by charge character, not sign, because that is what governs accuracy:
delocalized ions (protonated amines, carboxylates, phenolates) hug experiment; small
"hard" cations with concentrated charge (protonated methanol, H2S, alcohols) fall off — the
continuum-solvation limit (cf. FINDINGS #0 and the viologen dication).

Physics, not fitting: nothing is rescaled; computed vs measured, as-is.
"""
from __future__ import annotations
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from plot_style import apply_style, C, grid_xy, grid_y  # noqa: E402

import matplotlib.pyplot as plt   # noqa: E402
import numpy as np                # noqa: E402

FIG = ROOT / "results" / "figures"
CSV = ROOT / "results" / "solvation_validation.csv"
MN = ROOT / "data/raw/validation/MNSol/extracted/MNSolDatabase_v2012/MNSol_alldata.txt"


def _load():
    formula = {r["FileHandle"]: r["Formula"]
               for r in csv.DictReader(MN.open(), delimiter="\t")} if MN.exists() else {}
    rows = list(csv.DictReader(CSV.open()))
    for r in rows:
        r["calc"] = float(r["calc_kcal"]); r["exp"] = float(r["exp_kcal"]); r["err"] = float(r["err"])
        cc = r["charge_class"]
        if cc == "neutral":
            r["char"] = "neutral"
        elif cc == "cation" and "N" not in formula.get(r["id"], "N"):
            r["char"] = "concentrated"     # hard O/S-protonated cation
        else:
            r["char"] = "delocalized"      # amine cation, carboxylate/phenolate anion
    return rows


def _mae(sub):
    return np.mean([abs(r["err"]) for r in sub]) if sub else float("nan")


def _parity_panel(ax, groups, lo, hi, band, title):
    ax.plot([lo, hi], [lo, hi], color="#333333", lw=1.6, zorder=1)
    ax.fill_between([lo, hi], [lo - band, hi - band], [lo + band, hi + band],
                    color=C["grid"], alpha=0.30, zorder=0, label=f"±{band:g} kcal/mol")
    for label, color, marker, rows in groups:
        if not rows:
            continue
        ax.scatter([r["exp"] for r in rows], [r["calc"] for r in rows], s=95, c=color,
                   marker=marker, edgecolor="white", linewidth=0.8, zorder=3,
                   label=f"{label}  (n={len(rows)}, MAE {_mae(rows):.1f})")
    ax.set_xlim(lo, hi); ax.set_ylim(lo, hi); ax.set_aspect("equal")
    ax.set_xlabel("experimental ΔG$_{solv}$  (kcal/mol)")
    ax.set_ylabel("computed ΔG$_{solv}$  (kcal/mol)")
    ax.set_title(title); grid_xy(ax)
    ax.legend(loc="upper left", fontsize=12.5, framealpha=0.95)


def fig_parity(rows):
    apply_style()
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(16, 7.8))

    neu = [r for r in rows if r["charge_class"] == "neutral"]
    _parity_panel(
        axL,
        [("water · FreeSolv", C["measured"], "o", [r for r in neu if r["solvent"] == "water"]),
         ("acetonitrile · MNSol", C["good"], "^", [r for r in neu if r["solvent"] == "acetonitrile"])],
        lo=-25, hi=5, band=1.0, title="Neutral solutes")

    ions = [r for r in rows if r["charge_class"] in ("cation", "anion")]
    deloc = [r for r in ions if r["char"] == "delocalized"]
    conc = [r for r in ions if r["char"] == "concentrated"]
    lo = 5 * round((min(r["exp"] for r in ions) - 8) / 5) if ions else -105
    _parity_panel(
        axR,
        [("delocalized ions\n(amine⁺, RCOO⁻, ArO⁻)", C["dft"], "o", deloc),
         ("concentrated charge\n(hard O/S–H⁺)", C["warn"], "D", conc)],
        lo=lo, hi=0, band=5.0, title="Acetonitrile ions")

    fig.suptitle("Computed vs. experimental solvation free energy — SMD / ωB97M-V",
                 fontsize=22, fontweight="bold", y=1.02)
    fig.subplots_adjust(wspace=0.28, top=0.88)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "solvation_parity.png"); plt.close(fig)
    print(f"  -> {FIG/'solvation_parity.png'}")


def fig_mae_summary(rows):
    apply_style()
    cats = [
        ("water\nneutral", lambda r: r["solvent"] == "water" and r["charge_class"] == "neutral", C["good"]),
        ("MeCN\nneutral", lambda r: r["solvent"] == "acetonitrile" and r["charge_class"] == "neutral", C["good"]),
        ("MeCN ions\ndelocalized", lambda r: r["char"] == "delocalized", C["dft"]),
        ("MeCN cations\nconcentrated", lambda r: r["char"] == "concentrated", C["warn"]),
    ]
    labels, maes, ns, colors = [], [], [], []
    for lab, sel, col in cats:
        sub = [r for r in rows if sel(r)]
        if sub:
            labels.append(lab); maes.append(_mae(sub)); ns.append(len(sub)); colors.append(col)
    fig, ax = plt.subplots(figsize=(10.5, 7))
    xs = np.arange(len(labels))
    ax.bar(xs, maes, width=0.62, color=colors, edgecolor="#333333", linewidth=1.1)
    for x, m, n in zip(xs, maes, ns):
        ax.text(x, m + 0.3, f"{m:.1f}\n(n={n})", ha="center", va="bottom",
                fontsize=14, fontweight="bold")
    ax.axhline(1.0, color=C["good"], ls="--", lw=1.8, alpha=0.8)
    ax.text(len(labels) - 0.5, 1.2, "SMD neutral target ≈1 kcal/mol",
            ha="right", va="bottom", fontsize=12.5, color=C["good"])
    ax.set_xticks(xs); ax.set_xticklabels(labels)
    ax.set_ylabel("mean absolute error  (kcal/mol)")
    ax.set_title("ΔG$_{solv}$ vs. experiment — reliable except for concentrated charge")
    ax.set_ylim(0, max(maes) * 1.22 + 1); grid_y(ax)
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / "solvation_mae_summary.png"); plt.close(fig)
    print(f"  -> {FIG/'solvation_mae_summary.png'}")


def main():
    rows = _load()
    print(f"loaded {len(rows)} validation rows")
    fig_parity(rows)
    fig_mae_summary(rows)


if __name__ == "__main__":
    main()
