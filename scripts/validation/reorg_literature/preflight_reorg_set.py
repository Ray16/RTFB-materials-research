#!/usr/bin/env python
"""Pre-flight validation for a reorganization-energy set BEFORE any DFT is run.

Checks, per molecule:
  - CSV parses (pandas) with the expected columns
  - SMILES is valid, NEUTRAL (formal charge 0), and has ZERO radical electrons (closed-shell)
  - neutral electron count is EVEN (-> singlet); the ion for the couple is well defined
  - 3D embed (ETKDGv3 + MMFF) succeeds
  - reports formula + InChIKey so identity can be eyeballed / cross-checked
  - prints the exact (charge, spin) that will be used for neutral and ion

Exits nonzero if ANY molecule fails, so it can gate a run.

  python scripts/validation/reorg_literature/preflight_reorg_set.py --source <csv> --couple electron
"""
from __future__ import annotations
import argparse
import sys

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors as D
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")


def n_electrons(mol):
    return sum(a.GetAtomicNum() for a in mol.GetAtoms()) - Chem.GetFormalCharge(mol)


def embed_ok(smiles):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = AllChem.ETKDGv3(); p.randomSeed = 0xC0FFEE
    if AllChem.EmbedMolecule(m, p) != 0:
        p.useRandomCoords = True
        if AllChem.EmbedMolecule(m, p) != 0:
            return False
    try:
        AllChem.MMFFOptimizeMolecule(m)
    except Exception:
        pass
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True)
    ap.add_argument("--couple", choices=["hole", "electron"], required=True)
    a = ap.parse_args()

    try:
        df = pd.read_csv(a.source)
    except Exception as e:
        print(f"FAIL: CSV did not parse: {e}"); sys.exit(1)
    need = {"id", "smiles"}
    if not need.issubset(df.columns):
        print(f"FAIL: missing columns {need - set(df.columns)}"); sys.exit(1)
    print(f"CSV OK: {len(df)} rows, columns={list(df.columns)}\n")

    q_ion = 1 if a.couple == "hole" else -1
    hdr = f"{'id':20s} {'formula':10s} {'chg':>3s} {'rad':>3s} {'e-':>4s} {'embed':>5s} {'InChIKey':27s} status"
    print(hdr); print("-" * len(hdr))
    all_ok = True
    for _, r in df.iterrows():
        gid, smi = str(r["id"]), str(r["smiles"])
        m = Chem.MolFromSmiles(smi)
        if m is None:
            print(f"{gid:20s} {'-':10s} {'-':>3s} {'-':>3s} {'-':>4s} {'-':>5s} {'-':27s} FAIL:bad_smiles")
            all_ok = False; continue
        chg = Chem.GetFormalCharge(m)
        rad = D.CalcNumRadicalElectrons(m) if hasattr(D, "CalcNumRadicalElectrons") else \
              sum(a.GetNumRadicalElectrons() for a in m.GetAtoms())
        ne = n_electrons(m)
        emb = embed_ok(smi)
        probs = []
        if chg != 0:  probs.append("not_neutral")
        if rad != 0:  probs.append("has_radical")
        if ne % 2:    probs.append("odd_e-(neutral_not_singlet)")
        if not emb:   probs.append("embed_failed")
        status = "ok" if not probs else "FAIL:" + ",".join(probs)
        if probs: all_ok = False
        print(f"{gid:20s} {D.CalcMolFormula(m):10s} {chg:>3d} {rad:>3d} {ne:>4d} "
              f"{str(emb):>5s} {Chem.MolToInchiKey(m):27s} {status}")

    print("\nCouple setup that WILL be used:")
    print(f"  neutral : charge=0  spin=0 (singlet)")
    print(f"  ion     : charge={q_ion:+d} spin=1 (doublet)   [{a.couple}]")
    print(f"\n{'ALL CHECKS PASSED' if all_ok else 'SOME CHECKS FAILED -- do NOT run calculations'}")
    sys.exit(0 if all_ok else 2)


if __name__ == "__main__":
    main()
