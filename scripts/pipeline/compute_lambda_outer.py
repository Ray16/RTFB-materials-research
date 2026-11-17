#!/usr/bin/env python
"""Compute outer-sphere reorganization energy lambda_out for a redox couple, by BOTH continuum
methods (Born single-sphere + molecular-cavity nonequilibrium PCM), in any solvent (MeCN default).

Reusable: point it at a SMILES (or a CSV of them) and a solvent; it optimizes the reactant
geometry, computes the two redox-state densities at that fixed geometry, and reports lambda_out
from each method (1-body electrochemical AND 2-body self-exchange).

  python scripts/pipeline/compute_lambda_outer.py --smiles "O=C1C=CC(=O)C=C1" --couple reduction --solvent acetonitrile
  python scripts/pipeline/compute_lambda_outer.py --source results/reorg_anchors/lambda_out_set.csv --solvent acetonitrile

Notes / pitfalls respected:
  * lambda_in stays GAS-phase 4-point (computed elsewhere) -> orthogonal to this vertical lambda_out.
  * eps_op = n**2 used explicitly (MeCN 1.806); NOT defaulted to 1.
  * Born and PCM both give the 1-body value; --n-body 2 for self-exchange (x2, infinite-separation).
  * CPU PySCF is used so the PCM K/R operators and densities share AO ordering (correctness > speed).
"""
from __future__ import annotations
import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from redox.properties.lambda_outer import (born_lambda_outer, pcm_lambda_outer, sasa_radius_A,  # noqa: E402
                                solvent_constants, SOLVENTS)


def _atoms_from_smiles(smiles, seed=0xC0FFEE):
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = AllChem.ETKDGv3(); p.randomSeed = seed
    if AllChem.EmbedMolecule(m, p) != 0:
        p.useRandomCoords = True; AllChem.EmbedMolecule(m, p)
    AllChem.MMFFOptimizeMolecule(m)
    c = m.GetConformer()
    return [(a.GetSymbol(), tuple(c.GetAtomPosition(a.GetIdx()))) for a in m.GetAtoms()]


def _mol(atoms, charge, spin, basis):
    from pyscf import gto
    astr = "\n".join(f"{s} {x} {y} {z}" for s, (x, y, z) in atoms)
    return gto.M(atom=astr, basis=basis, charge=charge, spin=spin, verbose=0)


def _mf(mol, xc, backend="cpu"):
    if backend == "gpu":
        from gpu4pyscf import dft as gdft
        KS = gdft.RKS if mol.spin == 0 else gdft.UKS
        mf = KS(mol).density_fit()
    else:
        from pyscf import dft
        KS = dft.RKS if mol.spin == 0 else dft.UKS
        mf = KS(mol)
    mf.xc = xc; mf.conv_tol = 1e-9; mf.max_cycle = 200
    return mf


def _dm_numpy(mf):
    import numpy as np
    dm = mf.make_rdm1()
    return np.asarray(dm.get() if hasattr(dm, "get") else dm)   # cupy -> numpy for CPU PCM


def _opt_geom(atoms, charge, spin, xc, basis):
    from pyscf.geomopt.geometric_solver import optimize
    mopt = optimize(_mf(_mol(atoms, charge, spin, basis), xc), maxsteps=100)
    return [(mopt.atom_symbol(i), tuple(mopt.atom_coords()[i] * 0.52917721067))
            for i in range(mopt.natm)]


def _load_xyz(path):
    lines = Path(path).read_text().splitlines()
    n = int(lines[0].split()[0])
    atoms = []
    for ln in lines[2:2 + n]:
        p = ln.split()
        atoms.append((p[0], (float(p[1]), float(p[2]), float(p[3]))))
    return atoms


def _scf_dm(mol, xc, backend):
    """Converged density matrix; REJECTS an unconverged SCF (a non-stationary density gives a
    meaningless charge-transfer potential)."""
    mf = _mf(mol, xc, backend); mf.kernel()
    if not bool(mf.converged):
        raise RuntimeError(f"SCF not converged (q={mol.charge}, 2S={mol.spin}, {xc}/{mol.basis})")
    return _dm_numpy(mf)


def compute_couple(states, solv, smiles=None, xyz=None, xyz2=None, xc="b3lyp", basis="def2-svpd",
                   n_body=1, geom_method="mmff", backend="cpu"):
    """General outer-sphere lambda for an arbitrary couple.
    states = ((q1,m1),(q2,m2)); geometry from xyz (preferred) else MMFF/DFT from smiles.
    SCF on `backend` (gpu=gpu4pyscf); PCM contraction on validated CPU code (densities hand off)."""
    (q1, m1), (q2, m2) = states
    if xyz:
        geom = _load_xyz(xyz)
    else:
        atoms0 = _atoms_from_smiles(smiles)
        geom = atoms0 if geom_method == "mmff" else _opt_geom(atoms0, q1, m1, xc, basis)
    eps_s, eps_op = solv["eps_s"], solv["eps_op"]

    def pcm_at(g):
        mol1 = _mol(g, q1, m1, basis); dm1 = _scf_dm(mol1, xc, backend)
        mol2 = _mol(g, q2, m2, basis); dm2 = _scf_dm(mol2, xc, backend)
        return pcm_lambda_outer(mol1, dm1, dm2, eps_s, eps_op, n_body=n_body)

    # vertical lambda_o at BOTH equilibrium geometries when the second is given; the reported
    # value is their mean and the spread is kept as a diagnostic.
    pcm1 = pcm_at(geom)
    pcm2 = pcm_at(_load_xyz(xyz2)) if xyz2 else None
    pcm = pcm1 if pcm2 is None else 0.5 * (pcm1 + pcm2)
    a = sasa_radius_A(smiles) if smiles else _radius_from_geom(geom)
    dz = abs(q2 - q1)  # electrons transferred
    born = born_lambda_outer(a, eps_s, eps_op, z=dz, n_body=n_body)
    return dict(smiles=smiles or "", solvent=solv["name"], eps_s=eps_s, eps_op=eps_op,
                q1=q1, m1=m1, q2=q2, m2=m2, sasa_radius_A=round(a, 3),
                born_lambda_o_eV=round(born, 4), pcm_lambda_o_eV=round(pcm, 4),
                pcm_lambda_o_geom1_eV=round(pcm1, 4),
                pcm_lambda_o_geom2_eV=(round(pcm2, 4) if pcm2 is not None else ""),
                n_body=n_body, xc=xc, basis=basis, backend=backend,
                cavity="pyscf IEF-PCM, default radii/Lebedev", scf_converged=True,
                geom_source=("xyz:" + str(xyz) + (";" + str(xyz2) if xyz2 else "")) if xyz
                else f"smiles:{geom_method}")


def _radius_from_geom(geom):
    """crude effective radius from a geometry (A) when no SMILES for SASA: vdW-ish gyration."""
    import numpy as np
    xyz = np.array([p for _, p in geom]); c = xyz.mean(0)
    rg = np.sqrt(((xyz - c) ** 2).sum(1).mean())
    return rg + 1.5   # + ~vdW shell


def compute(smiles, couple, solv, xc="b3lyp", basis="def2-svpd", n_body=1,
            geom_method="mmff", backend="cpu"):
    """Neutral-singlet convenience wrapper (reduction/oxidation) for the generic shape set."""
    states = ((0, 0), (-1, 1)) if couple == "reduction" else ((0, 0), (1, 1))
    d = compute_couple(states, solv, smiles=smiles, xc=xc, basis=basis,
                       n_body=n_body, geom_method=geom_method, backend=backend)
    d["couple"] = couple
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smiles")
    ap.add_argument("--source", help="CSV with columns id,smiles,couple[,solvent]")
    ap.add_argument("--couple", default="reduction", choices=["reduction", "oxidation"])
    ap.add_argument("--solvent", default=None, help="name in SOLVENTS, or omit for config (MeCN)")
    ap.add_argument("--xc", default="b3lyp")
    ap.add_argument("--basis", default="def2-svpd", help="diffuse by default (anion densities)")
    ap.add_argument("--n-body", type=int, default=1)
    ap.add_argument("--geom", default="mmff", choices=["mmff", "dft"], help="geometry source")
    ap.add_argument("--backend", default="gpu", choices=["gpu", "cpu"], help="SCF backend")
    ap.add_argument("--out", default=None, help="write results CSV to this path")
    a = ap.parse_args()

    def solv_of(name):
        if not name:
            return solvent_constants()
        eps_s, n = SOLVENTS[name]
        return dict(eps_s=eps_s, eps_op=n * n, n=n, name=name)

    rows = []
    if a.smiles:
        rows.append({"id": "mol", "smiles": a.smiles, "couple": a.couple, "solvent": a.solvent})
    elif a.source:
        rows = list(csv.DictReader(open(a.source)))
    else:
        ap.error("give --smiles or --source")

    print(f"{'id':18s} {'couple':12s} {'solv':13s} {'a(A)':>6s} {'Born':>7s} {'PCM':>7s} {'PCM/Born':>8s}")
    print("-" * 78)
    out_rows = []
    for r in rows:
        rid = r.get("id", "?"); shape = r.get("shape", "")
        solvname = r.get("solvent") or a.solvent
        try:
            if r.get("q1") not in (None, ""):     # explicit-charge couple (our molecules, xyz geometry)
                states = ((int(r["q1"]), int(r["m1"])), (int(r["q2"]), int(r["m2"])))
                d = compute_couple(states, solv_of(solvname), smiles=r.get("smiles") or None,
                                   xyz=r.get("xyz") or None, xyz2=r.get("xyz2") or None,
                                   xc=a.xc, basis=a.basis,
                                   n_body=a.n_body, geom_method=a.geom, backend=a.backend)
                d["couple"] = r.get("couple", f"{d['q1']:+d}/{d['q2']:+d}")
            else:                                  # neutral-singlet reduction/oxidation
                d = compute(r["smiles"], r.get("couple", a.couple), solv_of(solvname),
                            xc=a.xc, basis=a.basis, n_body=a.n_body, geom_method=a.geom, backend=a.backend)
            ratio = d["pcm_lambda_o_eV"] / d["born_lambda_o_eV"] if d["born_lambda_o_eV"] else float("nan")
            print(f"{rid:18s} {str(d.get('couple',''))[:12]:12s} {d['solvent'][:13]:13s} "
                  f"{d['sasa_radius_A']:6.2f} {d['born_lambda_o_eV']:7.3f} {d['pcm_lambda_o_eV']:7.3f} "
                  f"{ratio:8.2f}", flush=True)
            out_rows.append(dict(id=rid, shape=shape, **d, pcm_over_born=round(ratio, 3)))
        except Exception as e:
            print(f"{rid:18s} ERROR: {type(e).__name__}: {e}", flush=True)
    if a.out and out_rows:
        cols = ["id", "shape", "smiles", "couple", "solvent", "eps_s", "eps_op",
                "sasa_radius_A", "born_lambda_o_eV", "pcm_lambda_o_eV", "pcm_over_born", "n_body",
                "pcm_lambda_o_geom1_eV", "pcm_lambda_o_geom2_eV", "xc", "basis", "backend",
                "cavity", "scf_converged", "geom_source"]
        with open(a.out, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore"); w.writeheader()
            w.writerows(out_rows)
        print(f"\nwrote {a.out}")


if __name__ == "__main__":
    main()
