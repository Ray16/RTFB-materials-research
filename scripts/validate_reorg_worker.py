#!/usr/bin/env python
"""Resumable worker: recompute D3TaLES reorganization energies under THEIR protocol
(single RDKit conformer -> B3LYP/6-31G* gas-opt of neutral + ion -> Nelsen 4-point lambda_i)
for a stride-sharded slice of results/d3tales_reorg_validation/molecules.csv.

One JSON per molecule in .../calc/<id>.json (skipped if already present), so parallel workers
and reruns never redo finished molecules. Errors are recorded, never fatal.

  CUDA_VISIBLE_DEVICES=<idx> python scripts/validate_reorg_worker.py --offset 0 --stride 6
"""
from __future__ import annotations
import argparse, json, traceback
from pathlib import Path
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")
from pyscf import gto
from pyscf.geomopt.geometric_solver import optimize
from gpu4pyscf import dft
from ase import Atoms

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results" / "d3tales_reorg_validation"
CALC = BASE / "calc"
HARTREE_EV = 27.211386245988
XC, BASIS = "b3lyp", "6-31g*"


def _embed(smiles):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = AllChem.ETKDGv3(); p.randomSeed = 0xC0FFEE
    if AllChem.EmbedMolecule(m, p) != 0:
        if AllChem.EmbedMolecule(m, AllChem.ETKDGv2()) != 0:
            raise RuntimeError("embed failed")
    AllChem.MMFFOptimizeMolecule(m)
    c = m.GetConformer()
    return Atoms(symbols=[a.GetSymbol() for a in m.GetAtoms()],
                 positions=np.array([list(c.GetAtomPosition(i)) for i in range(m.GetNumAtoms())]))


def _mol(atoms, charge, spin):
    astr = "\n".join(f"{s} {p[0]} {p[1]} {p[2]}"
                     for s, p in zip(atoms.get_chemical_symbols(), atoms.positions))
    return gto.M(atom=astr, basis=BASIS, charge=charge, spin=spin, verbose=0)


def _mf(mol):
    # density_fit (RI-J): ~5x faster on GPU with IDENTICAL reorg energy -- the DF error cancels
    # in same-molecule energy differences (verified: duroquinone lambda_e 0.512 == 0.512, 219s->44s).
    mf = (dft.RKS if mol.spin == 0 else dft.UKS)(mol).density_fit()
    mf.xc = XC; mf.conv_tol = 1e-9; mf.max_cycle = 300
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


def _four_point(neu_at, E_O_at_O, ion_charge, ion_spin):
    ion_at, E_R_at_R = _gas_opt(neu_at, ion_charge, ion_spin)
    E_O_at_R = _e(ion_at, 0, 0)
    E_R_at_O = _e(neu_at, ion_charge, ion_spin)
    return (E_O_at_R - E_O_at_O) + (E_R_at_O - E_R_at_R)


def process(row):
    out = CALC / f"{row['id']}.json"
    if out.exists():
        try:
            if json.loads(out.read_text()).get("status") in ("ok", "partial"):
                return "skip"
        except Exception:
            pass
    rec = dict(id=row["id"], smiles=row["smiles"], family=row["family"],
               d3_electron=row["d3_electron"], d3_hole=row["d3_hole"],
               n_atoms=int(row["n_atoms"]), n_rot=int(row["n_rot"]),
               our_electron=None, our_hole=None, status="ok", error="")
    try:
        neu = _embed(row["smiles"])
        neu_at, E_O_at_O = _gas_opt(neu, 0, 0)
        if pd.notna(row["d3_electron"]):
            try:
                rec["our_electron"] = round(_four_point(neu_at, E_O_at_O, -1, 1), 4)
            except Exception as e:
                rec["status"] = "partial"; rec["error"] += f"anion:{type(e).__name__} "
        if pd.notna(row["d3_hole"]):
            try:
                rec["our_hole"] = round(_four_point(neu_at, E_O_at_O, +1, 1), 4)
            except Exception as e:
                rec["status"] = "partial"; rec["error"] += f"cation:{type(e).__name__} "
    except Exception as e:
        rec["status"] = "fail"; rec["error"] = f"{type(e).__name__}: {str(e)[:120]}"
        rec["trace"] = traceback.format_exc()[-400:]
    CALC.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=2))
    return rec["status"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--offset", type=int, default=0)
    ap.add_argument("--stride", type=int, default=1)
    a = ap.parse_args()
    df = pd.read_csv(BASE / "molecules.csv")
    mine = df.iloc[a.offset::a.stride]
    print(f"[worker {a.offset}/{a.stride}] {len(mine)} molecules", flush=True)
    for i, (_, row) in enumerate(mine.iterrows()):
        st = process(row)
        print(f"[w{a.offset}] {i+1}/{len(mine)} {row['id']} {row['family']} -> {st}", flush=True)


if __name__ == "__main__":
    main()
