#!/usr/bin/env python
"""Recompute the Nelsen 4-point inner-sphere lambda_i from the dihedral-scan GLOBAL minima
(vs the single-shot optimizations), to test whether the D3TaLES gap is a missed conformer.

Reads results/dihedral_scan/<gid>_{neu,anion}.json (each has the state's global-min geometry
and its B3LYP/6-31G* gas energy), computes the two CROSS single points at the same level, and
prints the conformer-converged lambda_i for each molecule that has both files.
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
from pyscf import gto
from gpu4pyscf import dft
from ase import Atoms

ROOT = Path(__file__).resolve().parents[2]
HARTREE_EV = 27.211386245988
XC, BASIS = "b3lyp", "6-31g*"
SCAN = ROOT / "results" / "dihedral_scan"

# (charge, spin) per state; D3TaLES electron_reorg + our single-shot gas B3LYP for reference
META = {
    "mophquinone_sa":  dict(name="2-(4-methoxyphenyl)-1,4-benzoquinone (80LERT)",
                            d3=0.8716, single_shot=0.502),
    "dmophquinone_sa": dict(name="2-(2,5-dimethoxyphenyl)-5-methoxy-1,4-benzoquinone (05TPTO)",
                            d3=0.4587, single_shot=0.726),
}
CS = {"neu": (0, 0), "anion": (-1, 1)}


def _atoms(xyz):
    return Atoms(symbols=[r[0] for r in xyz], positions=[r[1:] for r in xyz])


def _energy(atoms, charge, spin):
    astr = "\n".join(f"{s} {p[0]} {p[1]} {p[2]}"
                     for s, p in zip(atoms.get_chemical_symbols(), atoms.positions))
    mol = gto.M(atom=astr, basis=BASIS, charge=charge, spin=spin, verbose=0)
    mf = (dft.RKS if spin == 0 else dft.UKS)(mol)
    mf.xc = XC; mf.conv_tol = 1e-9; mf.max_cycle = 200
    return float(mf.kernel()) * HARTREE_EV


def main():
    for gid, m in META.items():
        fn = SCAN / f"{gid}_neu.json"; fa = SCAN / f"{gid}_anion.json"
        if not (fn.exists() and fa.exists()):
            print(f"[skip] {gid}: scans not both done yet"); continue
        neu = json.loads(fn.read_text())["global_min"]
        an = json.loads(fa.read_text())["global_min"]
        neu_at = _atoms(neu["xyz"]); an_at = _atoms(an["xyz"])

        E_O_at_O = neu["E_hartree"] * HARTREE_EV          # neutral @ neutral-min
        E_R_at_R = an["E_hartree"] * HARTREE_EV            # anion   @ anion-min
        E_O_at_R = _energy(an_at, *CS["neu"])              # neutral @ anion-min  (cross)
        E_R_at_O = _energy(neu_at, *CS["anion"])           # anion   @ neutral-min (cross)

        lam = (E_O_at_R - E_O_at_O) + (E_R_at_O - E_R_at_R)
        print(f"\n=== {m['name']} ===")
        print(f"  neutral global-min dihedral = {neu['dihedral']:6.1f} deg")
        print(f"  anion   global-min dihedral = {an['dihedral']:6.1f} deg")
        print(f"  conformer-converged lambda_i = {lam:.3f} eV")
        print(f"    vs our single-shot gas B3LYP = {m['single_shot']:.3f} eV")
        print(f"    vs D3TaLES electron_reorg    = {m['d3']:.3f} eV  (diff {lam-m['d3']:+.3f})")


if __name__ == "__main__":
    main()
