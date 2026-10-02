#!/usr/bin/env python
"""Relaxed dihedral scan of the biaryl (quinone-ring <-> phenyl) torsion for one
(molecule, charge-state), at B3LYP/6-31G* gas -- to find each state's GLOBAL minimum and
test whether the D3TaLES-vs-ours reorg discrepancy is a missed-conformer artifact.

For each angle in 0..180 deg we do a geomeTRIC dihedral-constrained optimization (all other
DOF relaxed), walking the scan from the previous geometry. We then run an UNCONSTRAINED opt
from the lowest scan point to get the true global minimum for that state.

Output -> results/dihedral_scan/<gid>_<state>.json  (angles, energies, min geom + energy).

Usage (one job per GPU):
  CUDA_VISIBLE_DEVICES=1 python scripts/scan_dihedral_reorg.py --gid mophquinone_sa --state neu
  CUDA_VISIBLE_DEVICES=2 python scripts/scan_dihedral_reorg.py --gid mophquinone_sa --state anion
"""
from __future__ import annotations
import argparse, json, tempfile, os
from pathlib import Path
import numpy as np
from pyscf import gto
from pyscf.geomopt.geometric_solver import optimize
from gpu4pyscf import dft
from ase.io import read
from ase import Atoms

ROOT = Path(__file__).resolve().parents[2]
HARTREE_EV = 27.211386245988
XC, BASIS = "b3lyp", "6-31g*"

# biaryl dihedral atoms (0-based, matches RDKit build order == opt.xyz order)
DIH = {
    "mophquinone_sa":  [4, 5, 6, 7],
    "dmophquinone_sa": [4, 6, 7, 8],
}
STATE = {  # (charge, spin, subdir with the starting geometry)
    "neu":   (0,  0, "neu"),
    "anion": (-1, 1, "red1"),
}


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


def _energy(atoms, charge, spin):
    return float(_mf(_mol(atoms, charge, spin)).kernel())


def _dihedral(c, idx):
    p0, p1, p2, p3 = (c[i] for i in idx)
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / np.linalg.norm(b1)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    x, y = np.dot(v, w), np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))


def _constrained_opt(start_atoms, charge, spin, dih_idx, angle_deg, maxsteps=60):
    # geomeTRIC constraint files use 1-BASED atom indices.
    a, b, c, d = (i + 1 for i in dih_idx)
    txt = f"$set\ndihedral {a} {b} {c} {d} {angle_deg:.4f}\n"
    with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
        f.write(txt); cpath = f.name
    try:
        mol0 = _mol(start_atoms, charge, spin)
        mol_opt = optimize(_mf(mol0), constraints=cpath, maxsteps=maxsteps)
    finally:
        os.unlink(cpath)
    at = _atoms(mol_opt)
    return at, _energy(at, charge, spin)


def _unconstrained_opt(start_atoms, charge, spin, maxsteps=100):
    mol_opt = optimize(_mf(_mol(start_atoms, charge, spin)), maxsteps=maxsteps)
    at = _atoms(mol_opt)
    return at, _energy(at, charge, spin)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gid", required=True, choices=list(DIH))
    ap.add_argument("--state", required=True, choices=list(STATE))
    ap.add_argument("--test", action="store_true", help="one angle only, verify constraint")
    args = ap.parse_args()

    charge, spin, sub = STATE[args.state]
    dih_idx = DIH[args.gid]
    start = read(str(ROOT / "calcs" / "dft" / args.gid / sub / "opt.xyz"))

    outdir = ROOT / "results" / "dihedral_scan"; outdir.mkdir(parents=True, exist_ok=True)

    if args.test:
        target = 90.0
        at, E = _constrained_opt(start, charge, spin, dih_idx, target, maxsteps=40)
        got = _dihedral(at.get_positions(), dih_idx)
        print(f"[test] {args.gid}/{args.state}: target={target:.1f} achieved={got:.1f} "
              f"(|err|={abs(abs(got)-target):.2f}) E={E:.6f}  "
              f"{'OK' if abs(abs(got)-target)<2 else 'CONSTRAINT FAILED'}", flush=True)
        return

    angles = list(np.arange(0.0, 180.1, 15.0))
    recs = []
    cur = start
    for ang in angles:
        at, E = _constrained_opt(cur, charge, spin, dih_idx, ang, maxsteps=60)
        got = _dihedral(at.get_positions(), dih_idx)
        recs.append({"angle": ang, "achieved": float(got), "E_hartree": E,
                     "xyz": [[s, *map(float, p)] for s, p in
                             zip(at.get_chemical_symbols(), at.get_positions())]})
        cur = at  # walk the scan
        print(f"  {args.gid}/{args.state}  ang={ang:5.1f}  got={got:7.2f}  E={E:.6f}", flush=True)

    imin = int(np.argmin([r["E_hartree"] for r in recs]))
    print(f"  lowest scan point: angle={recs[imin]['angle']:.1f} "
          f"E={recs[imin]['E_hartree']:.6f}; running unconstrained opt from it...", flush=True)
    min_start = Atoms(symbols=[r[0] for r in recs[imin]["xyz"]],
                      positions=[r[1:] for r in recs[imin]["xyz"]])
    gmin_at, gmin_E = _unconstrained_opt(min_start, charge, spin)
    gmin_dih = _dihedral(gmin_at.get_positions(), dih_idx)
    print(f"  GLOBAL MIN: dihedral={gmin_dih:.1f} E={gmin_E:.6f}", flush=True)

    out = {"gid": args.gid, "state": args.state, "charge": charge, "spin": spin,
           "xc": XC, "basis": BASIS, "dih_idx": dih_idx, "scan": recs,
           "global_min": {"dihedral": float(gmin_dih), "E_hartree": gmin_E,
                          "xyz": [[s, *map(float, p)] for s, p in
                                  zip(gmin_at.get_chemical_symbols(), gmin_at.get_positions())]}}
    fp = outdir / f"{args.gid}_{args.state}.json"
    fp.write_text(json.dumps(out, indent=2))
    print(f"  wrote {fp}", flush=True)


if __name__ == "__main__":
    main()
