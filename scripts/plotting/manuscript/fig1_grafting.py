#!/usr/bin/env python
"""Figure 1 (manuscript, JCIM): the grafted-motif abstraction.

Schematic, CPU-only (RDKit + matplotlib), PNG only (300 dpi). Draws, for a single
real example (5,10-dihydrophenazine, DMPZ, the showcase p-type catholyte motif):

  (a) redox motif (both N-H) + 4-methylbenzyl chloride cap (Merrifield proxy)
      --graft--> the exact capped-graft model actually computed. Only ONE nitrogen
      (N5) carries the Merrifield backbone graft; the other (N10) carries a small
      ethyl TUNING cap (not a second graft).
  (b) the per-motif redox-state ladder {neutral, ox1 (+1), ox2 (+2)} as formal
      charge/spin states (NOT computed spin density).

The cap `Cc1ccc(C[*:1])cc1` and the motif fragment `[*:1]N1c2ccccc2N(CC)c2ccccc21`
are taken verbatim from config/custom_candidates.py (id exp_dmpz_mb_ethyl), so the
depicted capped graft is byte-identical to the registered model.

Style: plot_style.apply_style() + an 18 pt floor on every text element; zero
overlaps. Caveat on-figure (AGENTS.md): capped graft != soluble parent != oligomer.
"""
from __future__ import annotations
import io
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from plot_style import apply_style  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path(__file__).resolve().parents[3]
FIGDIR = ROOT / "manuscript" / "figures"
FIGDIR.mkdir(parents=True, exist_ok=True)

# --- real chemistry, verbatim from config/custom_candidates.py -------------
CAP_SCAFFOLD = "Cc1ccc(C[*:1])cc1"          # 4-methylbenzyl (Merrifield proxy)
MOTIF_FRAG = "[*:1]N1c2ccccc2N(CC)c2ccccc21"  # N'-ethyl-5,10-dihydrophenazine
CAP_CHLORIDE = "Cc1ccc(CCl)cc1"             # 4-methylbenzyl chloride (reagent)
MOTIF_PARENT = "N1c2ccccc2Nc2ccccc21"       # 5,10-dihydrophenazine (both N-H), flat 6-6-6

GRAFT_HL = (0.95, 0.65, 0.30)   # orange: the N5 Merrifield backbone graft
REDOX_HL = (0.70, 0.55, 0.85)   # purple: the two redox-active N centers


def _join(cap: str, frag: str) -> str:
    """Form the exact capped graft by fusing the two [*:1] dummy atoms."""
    combo = Chem.MolFromSmiles(cap + "." + frag)
    rw = Chem.RWMol(combo)
    dummies = [a.GetIdx() for a in rw.GetAtoms() if a.GetAtomicNum() == 0]
    assert len(dummies) == 2, "need exactly two attachment handles"
    nbrs = [n.GetIdx() for d in dummies for n in rw.GetAtomWithIdx(d).GetNeighbors()]
    rw.AddBond(nbrs[0], nbrs[1], Chem.BondType.SINGLE)
    for d in sorted(dummies, reverse=True):
        rw.RemoveAtom(d)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)


CAPPED_GRAFT = _join(CAP_SCAFFOLD, MOTIF_FRAG)


def _graft_attachment(mol):
    """Return (N5 graft-nitrogen idx, benzylic CH2 idx) for the 4-methylbenzyl
    backbone graft; N5 is the ring N whose exocyclic CH2 bears an aromatic ring
    (the tolyl cap), distinguishing it from the N10 ethyl tuning cap."""
    for a in mol.GetAtoms():
        if a.GetAtomicNum() != 7 or not a.IsInRing():
            continue
        for nb in a.GetNeighbors():
            if nb.GetAtomicNum() == 6 and not nb.IsInRing():
                if any(nn.GetIdx() != a.GetIdx() and nn.GetIsAromatic()
                       for nn in nb.GetNeighbors()):
                    return a.GetIdx(), nb.GetIdx()
    raise ValueError("graft nitrogen not found")


def _ring_nitrogens(mol):
    return [a.GetIdx() for a in mol.GetAtoms()
            if a.GetAtomicNum() == 7 and a.IsInRing()]


def _render(mol, w=720, h=560, hl_atoms=None, hl_bonds=None, hl_color=GRAFT_HL):
    """RDKit -> RGBA numpy array (high-res raster for crisp embedding)."""
    d = rdMolDraw2D.MolDraw2DCairo(w, h)
    opt = d.drawOptions()
    opt.bondLineWidth = 2
    opt.minFontSize = 22
    opt.maxFontSize = 30
    opt.padding = 0.08
    kw = {}
    if hl_atoms:
        kw["highlightAtoms"] = list(hl_atoms)
        kw["highlightAtomColors"] = {i: hl_color for i in hl_atoms}
    if hl_bonds:
        kw["highlightBonds"] = list(hl_bonds)
        kw["highlightBondColors"] = {i: hl_color for i in hl_bonds}
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol, **kw)
    d.FinishDrawing()
    from PIL import Image
    img = Image.open(io.BytesIO(d.GetDrawingText()))
    return np.asarray(img.convert("RGBA"))


def _place(fig, left, bottom, width, height, img, title=None, sub=None):
    ax = fig.add_axes([left, bottom, width, height])
    ax.imshow(img)
    ax.axis("off")
    if title:
        ax.set_title(title, fontsize=18, fontweight="bold", pad=3)
    if sub:
        ax.text(0.5, -0.08, sub, transform=ax.transAxes, ha="center",
                va="top", fontsize=18)
    return ax


def _arrow(fig, x0, x1, y, label_top=None, label_bot=None):
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_zorder(-1)
    ax.annotate("", xy=(x1, y), xytext=(x0, y),
                arrowprops=dict(arrowstyle="-|>", lw=2.4, color="#3A3A3A"))
    xm = 0.5 * (x0 + x1)
    if label_top:
        ax.text(xm, y + 0.022, label_top, ha="center", va="bottom", fontsize=18)
    if label_bot:
        ax.text(xm, y - 0.024, label_bot, ha="center", va="top", fontsize=18,
                color="#555555")
    return ax


def main():
    apply_style()
    graft = Chem.MolFromSmiles(CAPPED_GRAFT)
    motif = Chem.MolFromSmiles(MOTIF_PARENT)
    cap = Chem.MolFromSmiles(CAP_CHLORIDE)

    n5, ch2 = _graft_attachment(graft)
    graft_bond = graft.GetBondBetweenAtoms(n5, ch2).GetIdx()
    ring_n = _ring_nitrogens(graft)

    fig = plt.figure(figsize=(16.0, 14.4))

    # ---- Panel a: grafting abstraction (top) ----
    fig.text(0.015, 0.980, "a  Exact grafted-motif model", fontsize=22,
             fontweight="bold", va="top")
    rowy, rowh = 0.735, 0.205
    _place(fig, 0.015, rowy, 0.235, rowh, _render(motif),
           title="Redox motif",
           sub="5,10-dihydrophenazine\n(DMPZ; both N–H)")
    _place(fig, 0.305, rowy, 0.235, rowh, _render(cap),
           title="Merrifield-proxy cap", sub="4-methylbenzyl chloride")
    _place(fig, 0.645, rowy, 0.335, rowh,
           _render(graft, w=880, hl_atoms=[n5, ch2], hl_bonds=[graft_bond],
                   hl_color=GRAFT_HL),
           title="Capped graft (DMPZ, computed)",
           sub="N-(4-methylbenzyl)-N′-ethyl-DMPZ\n(N5 graft highlighted)")
    _arrow(fig, 0.255, 0.300, rowy + 0.5 * rowh, label_top="+")
    _arrow(fig, 0.545, 0.640, rowy + 0.5 * rowh,
           label_top="graft", label_bot="−HCl")

    fig.text(0.50, rowy - 0.085,
             "4-methylbenzyl cap = Merrifield (poly(vinylbenzyl)) backbone proxy;\n"
             "ring methyl (para) stands in for the polymer chain.",
             ha="center", va="top", fontsize=18, color="#333333")

    # ---- Panel b: redox-state ladder (middle) ----
    fig.text(0.015, 0.540, "b  Per-motif redox states (p-type, vs Fc/Fc⁺)",
             fontsize=22, fontweight="bold", va="top")
    by, bh = 0.300, 0.185
    xs = [0.015, 0.360, 0.705]
    badges = ["0", "+1 •", "2+"]
    names = ["Neutral", "Radical cation", "Dication"]
    for x, name, badge in zip(xs, names, badges):
        axm = _place(fig, x, by, 0.265, bh,
                     _render(graft, hl_atoms=ring_n, hl_color=REDOX_HL),
                     title=name)
        axm.text(0.90, 0.90, badge, transform=axm.transAxes, ha="center",
                 va="center", fontsize=20, fontweight="bold", color="#4a2a6a",
                 bbox=dict(boxstyle="circle,pad=0.3", fc="#EDE4F5",
                           ec="#762A83", lw=1.8))
    _arrow(fig, 0.283, 0.358, by + 0.5 * bh, label_top="−e⁻")
    _arrow(fig, 0.628, 0.703, by + 0.5 * bh, label_top="−e⁻")

    # ---- mechanistic note (between panel b and banner) ----
    fig.text(0.50, 0.270,
             "Both N are substituted (N,N′-disubstituted): this blocks N–H "
             "deprotonation on oxidation,\ngiving a reversible 2e⁻ couple. "
             "N5 bears the backbone graft; N10 bears a small tuning cap.",
             ha="center", va="top", fontsize=18, color="#333333")

    # ---- caveat banner (bottom) ----
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off"); ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.add_patch(plt.Rectangle((0.03, 0.020), 0.94, 0.135, fill=True,
                 facecolor="#FFF9E6", edgecolor="#E0C36A", lw=1.8, zorder=2))
    ax.text(0.5, 0.118,
            "capped graft  ≠  soluble parent  ≠  oligomer   —   "
            "three distinct models",
            ha="center", va="center", fontsize=20, fontweight="bold",
            color="#5a4a10", zorder=3)
    ax.text(0.5, 0.058,
            "The capped graft fixes local electronic structure and attachment chemistry; it does not\n"
            "represent chain packing, solubility, or cooperative effects. Badges are formal overall\n"
            "charge/spin labels, not QM spin density.",
            ha="center", va="center", fontsize=18, color="#555555", zorder=3)

    fig.savefig(FIGDIR / "fig1_grafting.png", dpi=300)
    plt.close(fig)
    print("capped graft SMILES:", CAPPED_GRAFT)
    print("N5 (graft) idx:", n5, "| CH2 idx:", ch2, "| ring N:", ring_n)
    print("  ->", FIGDIR / "fig1_grafting.png")


if __name__ == "__main__":
    main()
