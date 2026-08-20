#!/usr/bin/env python
"""2D structure gallery of every molecule in the calculation.

Reads the SINGLE source of truth (`library/manifest.csv`) so the gallery can
never drift from what the pipeline actually computes. Each unique molecule is
drawn once (RDKit 2D depiction) and labelled with its name, redox family, and
the set of charge states it is evaluated at. Two figures, matching the house
style (300 dpi, family colours, no overlap):

  molecules_candidates.png   the 6 decorated benzylic-site monomers (the screen)
  molecules_validation.png   the 6 known-E cores used for the validation gate

Depictions are schematic 2D graphs — e.g. ferrocene renders as Fe + two Cp
anions because RDKit cannot embed a metallocene sandwich in 2D. That is a
drawing limitation only; the 3D pipeline special-cases ferrocene's geometry.

  python scripts/plot_molecules.py

Physics, not fitting: nothing here touches a number; these are cosmetics.
"""
from __future__ import annotations
import csv
import io
import sys
from collections import OrderedDict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_style import apply_style, FAMILY_COLOR  # noqa: E402

import matplotlib.pyplot as plt          # noqa: E402
import matplotlib.image as mpimg         # noqa: E402
import numpy as np                       # noqa: E402

from rdkit import Chem                   # noqa: E402
from rdkit.Chem import Draw              # noqa: E402
from rdkit.Chem.Draw import rdMolDraw2D  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
FIGDIR = ROOT / "results" / "figures"
MANIFEST = ROOT / "library" / "manifest.csv"

TILE_PX = 560          # per-molecule depiction size (square), high-res for print
NCOLS = 2              # gallery columns

# Compact, readable display names (keyed by id) — the manifest names are verbose.
_PRETTY = {
    "pyridinium":           "N-benzylpyridinium",
    "cyanopyridinium":      "N-benzyl-4-cyanopyridinium",
    "viologen":             "Benzyl–methyl viologen",
    "phenothiazine":        "N-benzylphenothiazine",
    "tempo":                "4-(benzyloxy)-TEMPO",
    "anthraquinone":        "2-(benzyloxymethyl)anthraquinone",
    "ferrocene":            "Ferrocene  (internal reference)",
    "methyl_viologen":      "Methyl viologen",
    "tempo_parent":         "TEMPO",
    "phenothiazine_parent": "Phenothiazine (10H)",
    "anthraquinone_parent": "9,10-Anthraquinone",
    "methylpyridinium":     "N-methylpyridinium",
}


def _fmt_charge(q: int) -> str:
    q = int(q)
    if q == 0:
        return "0"
    return f"{q:+d}"


def _load_molecules():
    """Group manifest rows by id -> {name, family, smiles, charges[sorted]}."""
    mols = OrderedDict()
    with MANIFEST.open() as f:
        for row in csv.DictReader(f):
            rid = row["id"]
            m = mols.setdefault(rid, {
                "id": rid,
                "family": row["family"],
                "smiles": row["smiles"],
                "charges": set(),
            })
            try:
                m["charges"].add(int(row["charge"]))
            except (TypeError, ValueError):
                pass
    for m in mols.values():
        m["charges"] = sorted(m["charges"], reverse=True)  # most oxidized first
    return mols


def _draw_mol_png(smiles: str) -> np.ndarray:
    """Render a SMILES to an RGBA image array via RDKit's Cairo drawer."""
    mol = Chem.MolFromSmiles(smiles)
    d = rdMolDraw2D.MolDraw2DCairo(TILE_PX, TILE_PX)
    opts = d.drawOptions()
    opts.bondLineWidth = 2
    opts.padding = 0.10
    opts.clearBackground = True
    if mol is None:
        # Fallback: blank tile (should not happen for our curated set).
        d.FinishDrawing()
    else:
        Chem.rdDepictor.Compute2DCoords(mol)
        rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
        d.FinishDrawing()
    png = d.GetDrawingText()
    return mpimg.imread(io.BytesIO(png), format="png")


def _gallery(mols, title, out_name):
    apply_style()
    n = len(mols)
    ncols = NCOLS
    nrows = int(np.ceil(n / ncols))

    fig, axes = plt.subplots(
        nrows, ncols, figsize=(6.6 * ncols, 6.2 * nrows), squeeze=False)
    fig.suptitle(title, fontsize=24, fontweight="bold", y=0.995)

    for idx, m in enumerate(mols):
        r, c = divmod(idx, ncols)
        ax = axes[r][c]
        img = _draw_mol_png(m["smiles"])
        ax.imshow(img)
        ax.axis("off")

        colour = FAMILY_COLOR.get(m["family"], "#888888")
        # coloured family strip framing each tile
        for spine in ("top", "bottom", "left", "right"):
            ax.spines[spine].set_visible(True)
            ax.spines[spine].set_color(colour)
            ax.spines[spine].set_linewidth(3.0)
        ax.set_frame_on(True)
        ax.set_xticks([]); ax.set_yticks([])

        name = _PRETTY.get(m["id"], m["id"])
        charges = ", ".join(_fmt_charge(q) for q in m["charges"])
        ax.set_title(name, fontsize=17, fontweight="bold", pad=8, color="#111111")
        ax.text(0.5, -0.045,
                f"{m['family']}   ·   charge states:  {charges}",
                transform=ax.transAxes, ha="center", va="top",
                fontsize=13.5, color=colour, style="italic")

    # blank any unused cells
    for idx in range(n, nrows * ncols):
        r, c = divmod(idx, ncols)
        axes[r][c].axis("off")

    fig.tight_layout(rect=[0, 0.0, 1, 0.975])
    FIGDIR.mkdir(parents=True, exist_ok=True)
    out = FIGDIR / out_name
    fig.savefig(out)
    plt.close(fig)
    print(f"  -> {out}   ({n} molecules)")


def main():
    mols = _load_molecules()
    candidates = [m for m in mols.values() if m["family"] != "validation"]
    validation = [m for m in mols.values() if m["family"] == "validation"]
    print(f"molecules: {len(mols)}  (candidates: {len(candidates)}  "
          f"validation: {len(validation)})")

    _gallery(candidates,
             "Candidate redox monomers  (benzylic-site decorations)",
             "molecules_candidates.png")
    _gallery(validation,
             "Validation cores  (known E$^\\circ$ — the accuracy gate)",
             "molecules_validation.png")


if __name__ == "__main__":
    main()
