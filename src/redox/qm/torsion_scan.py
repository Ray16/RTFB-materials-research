"""Optional torsion sampling for the reorg pipeline (opt-in via run_batch(torsion_scan=True)).

Floppy molecules whose redox-relevant soft mode is a rotatable bond joining two ring systems
(biaryl quinones, inter-pyridinium viologens) can have conformer-dependent reorganization
energies: a single-shot optimization may land the neutral and the ion in different torsional
minima, inflating lambda_i. We protect against this WITHOUT a bespoke scan engine by seeding the
existing multi-conformer machinery (run_batch: SMD-optimize N seeds, keep lowest E_smd) with
geometries rotated around that torsion. Each state then finds its global-min conformer using the
production DFT+SMD+diffuse settings, so reorg.py's 4-point uses global-min geometries.

Detection is geometry-based (bond perception from the seed's 3D coords), so it needs no SMILES
atom-order matching and is charge/spin agnostic. Rigid molecules (no ring-ring rotatable bond)
yield no torsion -> the caller falls back to normal seeding (a no-op).
"""
from __future__ import annotations
from pathlib import Path
import numpy as np


def _load_perceived(xyz_path):
    """RDKit mol with a conformer and geometry-perceived connectivity, or None on any failure."""
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
    except Exception:
        return None
    try:
        mol = Chem.MolFromXYZFile(str(xyz_path))
        if mol is None:
            return None
        mol = Chem.RWMol(mol)
        rdDetermineBonds.DetermineConnectivity(mol)   # bonds from 3D distances (order-free)
        Chem.GetSSSR(mol)
        return mol
    except Exception:
        return None


def primary_conjugated_torsion(mol):
    """(a, i, j, b) dihedral atom indices for the primary rotatable bond joining two ring
    systems (biaryl / inter-ring), or None if the molecule has no such bond (i.e. rigid)."""
    ri = mol.GetRingInfo()
    rings = [set(r) for r in ri.AtomRings()]
    if len(rings) < 2:
        return None
    best, best_score = None, -1
    for bond in mol.GetBonds():
        if bond.IsInRing():
            continue
        ai, bi = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        ra = next((r for r in rings if ai in r), None)
        rb = next((r for r in rings if bi in r), None)
        if ra is None or rb is None or ra is rb:
            continue
        na = next((n.GetIdx() for n in mol.GetAtomWithIdx(ai).GetNeighbors()
                   if n.GetIdx() in ra), None)
        nb = next((n.GetIdx() for n in mol.GetAtomWithIdx(bi).GetNeighbors()
                   if n.GetIdx() in rb), None)
        if na is None or nb is None:
            continue
        score = len(ra) + len(rb)      # prefer the largest ring-ring linkage
        if score > best_score:
            best_score, best = score, (na, ai, bi, nb)
    return best


def _write_xyz(mol, conf, path, comment=""):
    syms = [a.GetSymbol() for a in mol.GetAtoms()]
    pos = conf.GetPositions()
    lines = [str(len(syms)), comment]
    for s, p in zip(syms, pos):
        lines.append(f"{s} {p[0]:.8f} {p[1]:.8f} {p[2]:.8f}")
    Path(path).write_text("\n".join(lines) + "\n")


def torsion_rotated_seeds(base_xyz, n_angles, outdir):
    """Write `n_angles` copies of `base_xyz` with the primary ring-ring torsion set to evenly
    spaced angles; return their paths. Returns [base_xyz] unchanged if no torsion is detected
    or RDKit bond perception is unavailable (safe no-op for rigid molecules)."""
    base_xyz = Path(base_xyz)
    mol = _load_perceived(base_xyz)
    if mol is None:
        return [base_xyz]
    dih = primary_conjugated_torsion(mol)
    if dih is None:
        return [base_xyz]
    from rdkit.Chem import rdMolTransforms
    conf = mol.GetConformer()
    outdir = Path(outdir); outdir.mkdir(parents=True, exist_ok=True)
    paths = []
    for k, ang in enumerate(np.linspace(0.0, 360.0, int(n_angles), endpoint=False)):
        rdMolTransforms.SetDihedralDeg(conf, *dih, float(ang))  # absolute set, deterministic
        p = outdir / f"_tseed_{k:02d}.xyz"
        _write_xyz(mol, conf, p, comment=f"torsion seed {k} dih={dih} ang={ang:.0f}")
        paths.append(p)
    return paths
