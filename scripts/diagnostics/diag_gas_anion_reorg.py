#!/usr/bin/env python
"""Diagnostic: can we reproduce D3TaLES's electron reorganization energy for the two quinones
by matching THEIR protocol -- gas-phase geometry optimization at B3LYP/6-31G* (small, no
diffuse) -- vs our default (SMD-opt geometry, wB97M-V/def2-TZVP(D))?

The quinone radical anions are UNBOUND in gas (adiabatic EA < 0 in D3TaLES). A small no-diffuse
basis artificially confines the escaping electron -> a 'bound-like' anion geometry -> large
relaxation. A diffuse basis lets it escape -> neutral-like geometry -> small relaxation. This
isolates geometry-source + basis as the cause of the ours(0.52)-vs-D3TaLES(0.87) gap.

  PYTHONPATH=src OMP_NUM_THREADS=6 python scripts/diagnostics/diag_gas_anion_reorg.py
"""
from __future__ import annotations
from pathlib import Path

import numpy as np
from pyscf import gto
from pyscf.geomopt.geometric_solver import optimize
from gpu4pyscf import dft           # GPU backend (~15x faster); clean SCF settings only
from ase.io import read

ROOT = Path(__file__).resolve().parents[2]
HARTREE_EV = 27.211386245988

CASES = {
    "mophquinone_sa": dict(name="2-(4-methoxyphenyl)benzoquinone", d3=0.8716),
    "dmophquinone_sa": dict(name="2-(2,5-dimethoxyphenyl)-5-methoxyquinone", d3=0.4587),
}
XC, BASIS = "b3lyp", "6-31g*"   # D3TaLES-like protocol


def _mol(atoms, charge, spin, basis=BASIS):
    astr = "\n".join(f"{s} {p[0]} {p[1]} {p[2]}"
                     for s, p in zip(atoms.get_chemical_symbols(), atoms.positions))
    return gto.M(atom=astr, basis=basis, charge=charge, spin=spin, verbose=0)


def _mf(mol):
    # Clean settings only: level_shift/damp crash gpu4pyscf. 6-31g* (no diffuse) already
    # confines the unbound-anion electron, so the SCF converges without them.
    mf = (dft.RKS if mol.spin == 0 else dft.UKS)(mol)
    mf.xc = XC; mf.conv_tol = 1e-9; mf.max_cycle = 200
    return mf


def _e(atoms, charge, spin):
    mf = _mf(_mol(atoms, charge, spin)); return float(mf.kernel())


def _atoms(mol):
    from ase import Atoms
    B = 0.52917721067
    return Atoms(symbols=[mol.atom_symbol(i) for i in range(mol.natm)],
                 positions=mol.atom_coords() * B)


def gas_opt(start_atoms, charge, spin):
    mol = _mol(start_atoms, charge, spin)
    mol_opt = optimize(_mf(mol), maxsteps=100)
    return _atoms(mol_opt)


def run(gid):
    # start from our SMD-optimized geometries (good starting points)
    neu0 = read(str(ROOT / "calcs" / "dft" / gid / "neu" / "opt.xyz"))
    red0 = read(str(ROOT / "calcs" / "dft" / gid / "red1" / "opt.xyz"))
    print(f"  gas-opt neutral (q=0)  at {XC}/{BASIS} ...", flush=True)
    gneu = gas_opt(neu0, 0, 0)
    print(f"  gas-opt anion   (q=-1) at {XC}/{BASIS} (unbound!) ...", flush=True)
    gred = gas_opt(red0, -1, 1)
    # 4-point on the gas-opt geometries (same level)
    E_neu_neu = _e(gneu, 0, 0); E_neu_red = _e(gred, 0, 0)
    E_red_red = _e(gred, -1, 1); E_red_neu = _e(gneu, -1, 1)
    lam = ((E_neu_red - E_neu_neu) + (E_red_neu - E_red_red)) * HARTREE_EV
    return lam


def main():
    print(f"protocol: {XC}/{BASIS}, GAS-phase optimization (matching D3TaLES)\n")
    for gid, info in CASES.items():
        print(f"=== {info['name']} ({gid}) ===")
        try:
            lam = run(gid)
        except Exception as e:
            print(f"  FAILED: {type(e).__name__}: {str(e)[:150]}"); continue
        print(f"  our gas-opt B3LYP/6-31G* lambda_i(e) = {lam:.3f} eV")
        print(f"  D3TaLES electron_reorg               = {info['d3']:.3f} eV   (diff {lam-info['d3']:+.3f})")
        print(f"  [ref] our SMD-geom wB97M-V value     = "
              + {"mophquinone_sa": "0.525", "dmophquinone_sa": "0.671"}[gid] + " eV\n")


if __name__ == "__main__":
    main()
