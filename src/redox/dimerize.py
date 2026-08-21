"""Build radical-ion dimers to score dimerization stability: 2 R^q -> (R2)^(2q).

The dominant degradation of many redox organics (viologen radical cations above all) is
DIMERIZATION, not disproportionation. Two open-shell monomers (each doublet) pair their
SOMOs into a diamagnetic pi-dimer (pimer): closed-shell SINGLET, total charge 2q, held by
SOMO-SOMO pairing + dispersion against Coulomb repulsion and desolvation.

dG_dim = G(dimer) - 2 G(monomer R^q);  dG_dim > 0 => STABLE against dimerization.

Correctness notes:
  - product spin: 2 doublets -> singlet pimer (mult 1). We build/relax the singlet.
  - dispersion is essential for the pi-stack: opt uses D4, the wb97m-v SP uses VV10 (NLC).
  - solvation: the dimer is a concentrated 2q charge vs two separated q charges, so the
    continuum desolvation term is the least reliable piece (same q^2 issue as redox). The
    reaction is atom-conserving and both sides are like-charged, so cancellation is far
    better than an absolute potential, but a +2 pimer may still warrant explicit-solvation
    refinement (Tier 2). This is flagged, not hidden.
  - geometry: cofacial, antiparallel (180 deg about the stacking axis) offset stack is the
    standard viologen pimer motif. A single well-built start is a first estimate; multiple
    stack offsets should be sampled for a production number.

  PYTHONPATH=src python -m redox.dimerize --monomer-xyz calcs/dft/methyl_viologen/ox1/opt.xyz \
      --out calcs/uma/mv_dimer/seed.xyz
"""
from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.io import read, write


def plane_normal_and_centroid(atoms: Atoms):
    p = atoms.get_positions()
    sym = np.array(atoms.get_chemical_symbols())
    heavy = p[sym != "H"]
    c = heavy.mean(0)
    _, _, vt = np.linalg.svd(heavy - c)
    n = vt[2] / np.linalg.norm(vt[2])
    # in-plane long axis (largest principal component) for the longitudinal offset
    long_axis = vt[0] / np.linalg.norm(vt[0])
    return n, c, long_axis


def build_pi_dimer(monomer: Atoms, interplanar=3.4, offset=1.6, antiparallel=True):
    """Cofacial pi-stacked dimer: second monomer translated `interplanar` along the ring
    normal + `offset` along the in-plane long axis, optionally rotated 180 deg about the
    normal (antiparallel stack). Returns combined Atoms (monomer A first, then B)."""
    n, c, long_axis = plane_normal_and_centroid(monomer)
    A = monomer.copy()
    B = monomer.copy()
    posB = B.get_positions() - c                      # center at origin
    if antiparallel:
        # 180 deg rotation about the stacking normal n (Rodrigues, theta=pi): v -> 2(v.n)n - v
        posB = 2.0 * np.outer(posB @ n, n) - posB
    posB = posB + c + interplanar * n + offset * long_axis
    B.set_positions(posB)
    dimer = A + B
    return dimer


def min_inter_monomer_dist(dimer: Atoms, n_mono):
    p = dimer.get_positions()
    a, b = p[:n_mono], p[n_mono:]
    d = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1)
    return float(d.min())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--monomer-xyz", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--interplanar", type=float, default=3.4)
    ap.add_argument("--offset", type=float, default=1.6)
    ap.add_argument("--parallel", action="store_true", help="cofacial parallel (default antiparallel)")
    a = ap.parse_args()

    mono = read(a.monomer_xyz)
    dimer = build_pi_dimer(mono, a.interplanar, a.offset, antiparallel=not a.parallel)
    mind = min_inter_monomer_dist(dimer, len(mono))
    out = Path(a.out); out.parent.mkdir(parents=True, exist_ok=True)
    write(str(out), dimer)
    print(f"monomer atoms={len(mono)}  dimer atoms={len(dimer)}  "
          f"formula={dimer.get_chemical_formula()}")
    print(f"interplanar={a.interplanar} offset={a.offset} "
          f"{'parallel' if a.parallel else 'antiparallel'}  "
          f"min inter-monomer dist={mind:.2f} A")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
