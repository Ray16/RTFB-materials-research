#!/usr/bin/env python
"""Gallery of the SIX starting candidates (from starting_candidates/Candidates.xlsx), as the
Merrifield-grafted monomers we built. Landscape 3x2 for slides.

  PYTHONPATH=src python scripts/plot_starting_candidates.py  ->  results/figures/molecules_starting_candidates.png
"""
from __future__ import annotations
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import RDLogger; RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "library" / "manifest.csv"
OUT = ROOT / "results" / "figures" / "candidates" / "molecules_starting_candidates.png"
OUT.parent.mkdir(parents=True, exist_ok=True)

# the six xlsx candidates, in order, as the ids we built (methylviologen -> viologen)
ORDER = ["viologen", "ethylviologen", "pmdi", "ndi_ammonium", "mophquinone", "dmophquinone"]
PRETTY = {
    "viologen": "Methyl viologen", "ethylviologen": "Ethyl viologen",
    "pmdi": "Pyromellitic diimide (PMDI)", "ndi_ammonium": "Ammonium naphthalene diimide",
    "mophquinone": "Methoxyphenyl-benzoquinone", "dmophquinone": "Dimethoxyphenyl-methoxyquinone",
}
FAM_COLOR = {"pyridine-multi-e": "#CC79A7", "imide (n-type)": "#009E73",
             "quinone (n-type)": "#0072B2"}
TILE_PX = 560


def _load():
    mols = {}
    for row in csv.DictReader(MANIFEST.open()):
        rid = row["id"]
        m = mols.setdefault(rid, dict(id=rid, family=row["family"], smiles=row["smiles"], charges=set()))
        m["charges"].add(int(row["charge"]))
    return mols


def _draw(smiles):
    mol = Chem.MolFromSmiles(smiles)
    Chem.rdDepictor.SetPreferCoordGen(True)
    Chem.rdDepictor.Compute2DCoords(mol)
    d = rdMolDraw2D.MolDraw2DCairo(TILE_PX, TILE_PX)
    d.drawOptions().bondLineWidth = 2
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    import io
    return mpimg.imread(io.BytesIO(d.GetDrawingText()), format="png")


def main():
    mols = _load()
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"]})
    ncols, nrows = 3, 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.4 * ncols, 6.0 * nrows), squeeze=False)
    fig.suptitle("Starting candidates  (Merrifield-grafted monomers)",
                 fontsize=26, fontweight="bold", y=0.995)
    for idx, rid in enumerate(ORDER):
        r, c = divmod(idx, ncols)
        ax = axes[r][c]
        m = mols[rid]
        ax.imshow(_draw(m["smiles"])); ax.axis("off")
        col = FAM_COLOR.get(m["family"], "#888888")
        for s in ("top", "bottom", "left", "right"):
            ax.spines[s].set_visible(True); ax.spines[s].set_color(col); ax.spines[s].set_linewidth(3.0)
        ax.set_frame_on(True); ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(PRETTY[rid], fontsize=18, fontweight="bold", pad=8, color="#111111")
        charges = ", ".join(f"{q:+d}" if q else "0" for q in sorted(m["charges"], reverse=True))
        ax.text(0.5, -0.05, f"{m['family']}   ·   charge states:  {charges}",
                transform=ax.transAxes, ha="center", va="top", fontsize=14, color=col, style="italic")
    fig.tight_layout(rect=[0, 0, 1, 0.975])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=200)
    print("wrote", OUT)


if __name__ == "__main__":
    main()
