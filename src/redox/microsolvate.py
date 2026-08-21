"""Explicit microsolvation (cluster-continuum) for redox couples in acetonitrile.

Places N explicit MeCN molecules in the first solvation shell around a solute's charged
region, oriented physically (MeCN's electron-rich N toward cationic sites; its methyl H+
toward anionic sites), then the caller UMA-relaxes and DFT+SMD-optimizes the cluster inside
the SMD continuum. The SAME explicit shell is applied to every redox state of a couple so the
bulk-solvent contribution cancels and we recover the *differential specific* solvation the pure
continuum misses (the charge-dependent q^2 error seen in the OROP benchmark).

Physics, not fitting: this changes the solvation MODEL, not any parameter tuned to data.
"""
from __future__ import annotations
import numpy as np
from ase import Atoms

# MeCN template (RDKit/MMFF), N-C-C(methyl) axis along +x; N is the last-listed atom region.
# atoms: C(methyl), C(nitrile), N, H, H, H
_MECN_SYM = ["C", "C", "N", "H", "H", "H"]
_MECN_XYZ = np.array([
    [-0.4863, -0.0043,  0.0040],   # methyl C
    [ 0.9755,  0.0086, -0.0081],   # nitrile C
    [ 2.1357,  0.0189, -0.0177],   # N (electron-rich end)
    [-0.8748, -0.7191, -0.7275],
    [-0.8644, -0.2883,  0.9906],
    [-0.8858,  0.9842, -0.2414],
])
_N_IDX, _METHYLC_IDX = 2, 0


def _mecn_orientation():
    """Unit vector from methyl-C toward N (the molecular long axis, N end forward)."""
    v = _MECN_XYZ[_N_IDX] - _MECN_XYZ[_METHYLC_IDX]
    return v / np.linalg.norm(v)


def _rot_matrix(a, b):
    """Rotation matrix taking unit vector a onto unit vector b."""
    a = a / np.linalg.norm(a); b = b / np.linalg.norm(b)
    v = np.cross(a, b); c = float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-8:
        return np.eye(3) if c > 0 else -np.eye(3)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * (1.0 / (1.0 + c))


def place_mecn(target, approach_dir, dist=3.2, point="N"):
    """One MeCN whose coordinating end sits `dist` A from `target` along `approach_dir`.

    approach_dir points FROM the solute site OUTWARD to where the MeCN center goes.
    point='N'  -> N end faces the target (use near CATIONS).
    point='CH3'-> methyl end faces the target (use near ANIONS).
    """
    approach_dir = np.asarray(approach_dir, float)
    approach_dir /= np.linalg.norm(approach_dir)
    axis = _mecn_orientation()
    # we want the coordinating end pointing back toward the target, i.e. molecular axis
    # anti-parallel (N end) or parallel (methyl end) to approach_dir.
    want = -approach_dir if point == "N" else approach_dir
    R = _rot_matrix(axis, want)
    xyz = (_MECN_XYZ - _MECN_XYZ.mean(0)) @ R.T
    coord_idx = _N_IDX if point == "N" else _METHYLC_IDX
    # translate so the coordinating atom is `dist` from target along approach_dir
    coord_pos_target = target + approach_dir * dist
    xyz += coord_pos_target - xyz[coord_idx]
    return Atoms(symbols=_MECN_SYM, positions=xyz)


def plane_normal(positions):
    heavy = np.asarray(positions)
    c = heavy.mean(0)
    _, _, vt = np.linalg.svd(heavy - c)
    n = vt[2]
    return n / np.linalg.norm(n)


def microsolvate_around_sites(solute: Atoms, sites, normal, per_site=2,
                              dist=3.2, point="N"):
    """Add `per_site` MeCN near each site coordinate, on alternating faces (+/- normal).

    sites: list of 3-vectors (e.g. the cationic N+ positions). normal: ring-plane normal.
    Returns a new Atoms = solute + explicit shell (solute atoms first, order preserved).
    """
    cluster = solute.copy()
    normal = np.asarray(normal, float); normal /= np.linalg.norm(normal)
    for site in sites:
        for k in range(per_site):
            sign = 1.0 if k % 2 == 0 else -1.0
            cluster += place_mecn(np.asarray(site, float), sign * normal,
                                  dist=dist, point=point)
    return cluster


def min_intermolecular_dist(cluster: Atoms, n_solute):
    """Smallest distance between a solute atom and any solvent atom (clash check)."""
    p = cluster.get_positions()
    su, sv = p[:n_solute], p[n_solute:]
    if len(sv) == 0:
        return np.inf
    d = np.linalg.norm(su[:, None, :] - sv[None, :, :], axis=-1)
    return float(d.min())
