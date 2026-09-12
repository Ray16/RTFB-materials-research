#!/usr/bin/env python
"""Before- vs after-grafting structures of the SIX starting candidates, in one figure:
each row = one candidate, left = bare redox core (before grafting), right = Merrifield-grafted
monomer (after grafting onto the 4-methylbenzyl benzylic site).

  PYTHONPATH=src python scripts/plotting/plot_before_after_grafting.py
    -> results/figures/candidates/before_after_grafting.png
"""
from __future__ import annotations
import csv
import io
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import RDLogger; RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[2]
MANIFEST = ROOT / "library" / "manifest.csv"
OUT = ROOT / "results" / "figures" / "candidates" / "before_after_grafting.png"

# (before-core id, after-grafted id, pretty name, family color)
FAM = {"pyridine-multi-e": "#CC79A7", "imide (n-type)": "#009E73", "quinone (n-type)": "#0072B2"}
PAIRS = [
    ("methyl_viologen",  "viologen",       "Methyl viologen",                "pyridine-multi-e"),
    ("ethylviologen_sa", "ethylviologen",  "Ethyl viologen",                 "pyridine-multi-e"),
    ("pmdi_sa",          "pmdi",           "Pyromellitic diimide (PMDI)",    "imide (n-type)"),
    ("ndi_ammonium_sa",  "ndi_ammonium",   "Ammonium naphthalene diimide",   "imide (n-type)"),
    ("mophquinone_sa",   "mophquinone",    "Methoxyphenyl-benzoquinone",     "quinone (n-type)"),
    ("dmophquinone_sa",  "dmophquinone",   "Dimethoxyphenyl-methoxyquinone", "quinone (n-type)"),
]
TILE_PX = 560


def _smiles():
    out = {}
    for r in csv.DictReader(MANIFEST.open()):
        out.setdefault(r["id"], r["smiles"])
    return out


def _draw(smiles):
    mol = Chem.MolFromSmiles(smiles)
    Chem.rdDepictor.SetPreferCoordGen(True)
    Chem.rdDepictor.Compute2DCoords(mol)
    d = rdMolDraw2D.MolDraw2DCairo(TILE_PX, TILE_PX)
    d.drawOptions().bondLineWidth = 2
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    return mpimg.imread(io.BytesIO(d.GetDrawingText()), format="png")


def main():
    smi = _smiles()
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"]})
    nrows, ncols = len(PAIRS), 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.0 * ncols, 4.6 * nrows), squeeze=False)
    fig.suptitle("Before vs. after grafting  (six starting candidates)",
                 fontsize=26, fontweight="bold", y=0.997)
    # column headers
    axes[0][0].set_title("Before grafting\n(bare redox core)", fontsize=20, fontweight="bold", pad=14)
    axes[0][1].set_title("After grafting\n(Merrifield-grafted monomer)", fontsize=20, fontweight="bold", pad=14)

    for r, (bid, aid, name, fam) in enumerate(PAIRS):
        col = FAM.get(fam, "#888888")
        for c, sid in ((0, bid), (1, aid)):
            ax = axes[r][c]
            ax.imshow(_draw(smi[sid])); ax.set_xticks([]); ax.set_yticks([])
            for s in ("top", "bottom", "left", "right"):
                ax.spines[s].set_visible(True); ax.spines[s].set_color(col); ax.spines[s].set_linewidth(3.0)
        # candidate name as a left-edge row label
        axes[r][0].set_ylabel(name, fontsize=17, fontweight="bold", color=col, rotation=90,
                              labelpad=14, va="center")
    fig.tight_layout(rect=[0, 0, 1, 0.975])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=200, facecolor="white")
    plt.close(fig)
    print(f"  -> {OUT}")


if __name__ == "__main__":
    main()
