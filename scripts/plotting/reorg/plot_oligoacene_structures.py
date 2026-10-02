#!/usr/bin/env python
"""Companion figure to the parity plot: 2D structures of the oligoacene validation set,
labelled with molecule name and our computed inner-sphere hole lambda.

naphthalene, anthracene, tetracene, pentacene (our UB3LYP/6-311G**) + rubrene (UB3LYP/6-31G(d,p)).
Uses CoordGen for clean 2D layouts and composites with matplotlib so the lambda glyph renders.
PNG only (no PDF unless asked).

  python scripts/plotting/reorg/plot_oligoacene_structures.py
"""
from __future__ import annotations
import csv
import io
import json
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
import sys
sys.path.insert(0, str(ROOT / "scripts" / "plotting"))
import plot_style as ps  # noqa: E402

OUT = ROOT / "results" / "figures" / "reorg" / "oligoacene_structures.png"
ORDER = ["naphthalene", "anthracene", "tetracene", "pentacene", "rubrene"]
SMI_SRC = ROOT / "results/reorg_anchors/oligoacene_inner.csv"
ACENE_CALC = ROOT / "results/reorg_anchors/inner_acene/b3lyp_6311gdp/calc"
RUB_CALC = ROOT / "results/reorg_anchors/inner_acene/calc"
LEVEL = {"rubrene": "UB3LYP/6-31G(d,p)"}
DEFAULT_LEVEL = "UB3LYP/6-311G**"


def _our_hole(gid):
    calc = RUB_CALC if gid == "rubrene" else ACENE_CALC
    j = calc / f"{gid}.json"
    if j.exists():
        v = json.loads(j.read_text()).get("our_hole")
        return float(v) if v is not None else None
    return None


def _draw(mol, px=460):
    rdCoordGen.AddCoords(mol)                    # cleaner 2D layout than Compute2DCoords
    d = rdMolDraw2D.MolDraw2DCairo(px, px)
    o = d.drawOptions(); o.bondLineWidth = 2; o.padding = 0.10
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    return Image.open(io.BytesIO(d.GetDrawingText()))


def main():
    ps.apply_style()
    smiles = {r["id"]: r["smiles"] for r in csv.DictReader(open(SMI_SRC))}
    n = len(ORDER)
    fig, axes = plt.subplots(1, n, figsize=(3.0 * n, 3.4))
    for ax, gid in zip(axes, ORDER):
        img = _draw(Chem.MolFromSmiles(smiles[gid]))
        ax.imshow(img); ax.axis("off")
        ax.set_title(gid, fontsize=15, pad=4)
    fig.suptitle("Oligoacene validation set — 2D structures", fontsize=17, y=1.02)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=300, bbox_inches="tight")   # PNG only
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
