#!/usr/bin/env python
"""Before- vs after-grafting structures of the SIX starting candidates, in one figure:
each row = one candidate, left = the MINIMAL CAPPED ANALOGUE (the redox core with the tether
position capped by a methyl), right = the Merrifield benzyl-tether model (grafted onto the
4-methylbenzyl benzylic site).

NOT a literal before/after of the experimental molecule: the left column is a methyl-capped
analogue, not the actual pre-grafting precursor. For PMDI the sheet compound is an
N-2-pentyl/PEG-substituted diimide while the computed reference is N,N'-dimethyl PMDI, so the
comparison isolates a SUBSTITUENT effect (methyl -> benzyl) rather than the grafting reaction
itself. Compute the true precursor if a literal grafting comparison is needed.

  PYTHONPATH=src python scripts/plotting/candidates/plot_before_after_grafting.py
    -> results/figures/candidates/before_after_grafting.png
"""
from __future__ import annotations
import csv
import io
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import RDLogger; RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
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


def _charges():
    """id -> set of charge states, from the manifest."""
    out = {}
    for r in csv.DictReader(MANIFEST.open()):
        out.setdefault(r["id"], set()).add(int(r["charge"]))
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
    charges = _charges()
    plt.rcParams.update({"font.family": "sans-serif",
                         "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"]})
    nrows, ncols = len(PAIRS), 2
    fig, axes = plt.subplots(nrows, ncols, figsize=(6.0 * ncols + 3.2, 4.6 * nrows), squeeze=False)
    fig.suptitle("Before vs. after grafting  (six starting candidates)",
                 fontsize=26, fontweight="bold", y=0.997)
    # column headers
    axes[0][0].set_title("Minimal capped analogue\n(methyl at the tether site)",
                         fontsize=20, fontweight="bold", pad=14)
    axes[0][1].set_title("After grafting\n(Merrifield-grafted monomer)", fontsize=20, fontweight="bold", pad=14)

    for r, (bid, aid, name, fam) in enumerate(PAIRS):
        col = FAM.get(fam, "#888888")
        for c, sid in ((0, bid), (1, aid)):
            ax = axes[r][c]
            ax.imshow(_draw(smi[sid])); ax.set_xticks([]); ax.set_yticks([])
            for s in ("top", "bottom", "left", "right"):
                ax.spines[s].set_visible(True); ax.spines[s].set_color(col); ax.spines[s].set_linewidth(3.0)
        # candidate name as a HORIZONTAL row label left of the tiles (rotated multi-line
        # labels grow sideways into the structure tile); wrapped to fit the left margin
        label = textwrap.fill(name.replace("-", "- "), width=14, break_long_words=False).replace("- ", "-")
        axes[r][0].set_ylabel(label, fontsize=18, fontweight="bold", color=col, rotation=0,
                              labelpad=18, va="center", ha="right",
                              multialignment="right", linespacing=1.15)
        # family + redox charge states under the grafted (after) tile — carried over from
        # the (now-retired) standalone gallery so this figure supersedes it.
        qs = ", ".join(f"{q:+d}" if q else "0" for q in sorted(charges[aid], reverse=True))
        axes[r][1].set_xlabel(f"{fam}   ·   charge states:  {qs}",
                              fontsize=15, color=col, style="italic", labelpad=8)
    fig.tight_layout(rect=[0.0, 0, 1, 0.975])
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=200, facecolor="white")
    plt.close(fig)
    print(f"  -> {OUT}")


if __name__ == "__main__":
    main()
