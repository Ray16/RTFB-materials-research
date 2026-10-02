#!/usr/bin/env python
"""2D structures of the perfluoroacene electron/hole validation set (companion to the parity
figure), names only. PNG only.

  python scripts/plotting/reorg/plot_perfluoroacene_structures.py
"""
from __future__ import annotations
import csv, io
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
from rdkit import Chem
from rdkit.Chem import rdCoordGen
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
import sys; sys.path.insert(0, str(ROOT / "scripts" / "plotting"))
import plot_style as ps  # noqa: E402

OUT = ROOT / "results" / "figures" / "reorg" / "perfluoroacene_structures.png"
SRC = ROOT / "results" / "reorg_anchors" / "electron_acceptors.csv"
NAME = {"perfluorotetracene": "perfluorotetracene (PFT)", "perfluoropentacene": "perfluoropentacene (PFP)"}


def _draw(mol, px=560):
    rdCoordGen.AddCoords(mol)
    d = rdMolDraw2D.MolDraw2DCairo(px, px)
    o = d.drawOptions(); o.bondLineWidth = 2; o.padding = 0.10
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    return Image.open(io.BytesIO(d.GetDrawingText()))


def main():
    ps.apply_style()
    rows = list(csv.DictReader(open(SRC)))
    n = len(rows)
    fig, axes = plt.subplots(1, n, figsize=(4.2 * n, 4.2))
    if n == 1:
        axes = [axes]
    for ax, r in zip(axes, rows):
        img = _draw(Chem.MolFromSmiles(r["smiles"]))
        ax.imshow(img); ax.axis("off")
        ax.set_title(NAME.get(r["id"], r["id"]), fontsize=15, pad=4)
    fig.suptitle("Perfluoroacene validation set — 2D structures", fontsize=17, y=1.03)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=300, bbox_inches="tight")   # PNG only
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
