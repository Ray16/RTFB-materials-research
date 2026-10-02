#!/usr/bin/env python
"""Can our machinery reproduce D3TaLES electron_reorg for OTHER single-electron molecules
when we use THEIR protocol (single conformer -> B3LYP/6-31G* gas-opt of neutral + unbound
anion -> Nelsen 4-point lambda_i)? Tests rigid vs floppy molecules to localize the failure.

Single-shot on purpose: one RDKit conformer, no scan -- faithfully mimics a high-throughput
pipeline. Rigid molecules have one minimum (should match D3TaLES); floppy ones have several
(may diverge, like the two biaryl quinones).

  CUDA_VISIBLE_DEVICES=<idx> python scripts/diagnostics/diag_smiles_reorg.py --key anthraquinone
"""
from __future__ import annotations
import argparse
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from pyscf import gto
from pyscf.geomopt.geometric_solver import optimize
from gpu4pyscf import dft
from ase import Atoms

HARTREE_EV = 27.211386245988
XC, BASIS = "b3lyp", "6-31g*"

# D3TaLES quinones: id, SMILES, electron_reorg, rigidity tag
PANEL = {
    "anthraquinone":  dict(smiles="O=C1c2ccccc2C(=O)c2ccccc21",              d3=0.3035, kind="rigid (fused)"),
    "duroquinone":    dict(smiles="CC1=C(C)C(=O)C(C)=C(C)C1=O",             d3=0.5725, kind="rigid (methyls)"),
    "chlorophenylBQ": dict(smiles="O=C1C=CC(=O)C(c2ccc(Cl)cc2)=C1",         d3=0.4841, kind="floppy (biaryl)"),
    "styrylBQ":       dict(smiles="O=C1C=CC(=O)C(C=Cc2ccccc2)=C1",          d3=1.0130, kind="floppy (styryl)"),
}


def _embed(smiles, charge, spin):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = AllChem.ETKDGv3(); p.randomSeed = 0xC0FFEE
    AllChem.EmbedMolecule(m, p)
    AllChem.MMFFOptimizeMolecule(m)
    conf = m.GetConformer()
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    pos = np.array([list(conf.GetAtomPosition(i)) for i in range(m.GetNumAtoms())])
    return Atoms(symbols=syms, positions=pos)


def _mol(atoms, charge, spin):
    astr = "\n".join(f"{s} {p[0]} {p[1]} {p[2]}"
                     for s, p in zip(atoms.get_chemical_symbols(), atoms.positions))
    return gto.M(atom=astr, basis=BASIS, charge=charge, spin=spin, verbose=0)


def _mf(mol):
    mf = (dft.RKS if mol.spin == 0 else dft.UKS)(mol)
    mf.xc = XC; mf.conv_tol = 1e-9; mf.max_cycle = 200
    return mf


def _atoms(mol):
    B = 0.52917721067
    return Atoms(symbols=[mol.atom_symbol(i) for i in range(mol.natm)],
                 positions=mol.atom_coords() * B)


def _gas_opt(start, charge, spin):
    mol_opt = optimize(_mf(_mol(start, charge, spin)), maxsteps=100)
    at = _atoms(mol_opt)
    return at, float(_mf(_mol(at, charge, spin)).kernel()) * HARTREE_EV


def _e(atoms, charge, spin):
    return float(_mf(_mol(atoms, charge, spin)).kernel()) * HARTREE_EV


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--key", required=True, choices=list(PANEL))
    a = ap.parse_args(); info = PANEL[a.key]
    print(f"=== {a.key}  [{info['kind']}]  {info['smiles']} ===", flush=True)
    start = _embed(info["smiles"], 0, 0)
    print("  gas-opt neutral (q=0)...", flush=True)
    neu_at, E_O_at_O = _gas_opt(start, 0, 0)
    print("  gas-opt anion (q=-1, unbound)...", flush=True)
    an_at, E_R_at_R = _gas_opt(neu_at, -1, 1)
    E_O_at_R = _e(an_at, 0, 0)      # neutral @ anion geom
    E_R_at_O = _e(neu_at, -1, 1)    # anion   @ neutral geom
    lam = (E_O_at_R - E_O_at_O) + (E_R_at_O - E_R_at_R)
    print(f"  our B3LYP/6-31G* single-shot lambda_i = {lam:.3f} eV", flush=True)
    print(f"  D3TaLES electron_reorg                = {info['d3']:.3f} eV  "
          f"(diff {lam-info['d3']:+.3f})", flush=True)


if __name__ == "__main__":
    main()
