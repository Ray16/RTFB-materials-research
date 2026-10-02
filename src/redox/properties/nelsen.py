"""Shared Nelsen 4-point inner-sphere reorganization-energy helpers (fresh-geometry variant).

Used by the D3TaLES reorg validation worker and the SMILES diagnostics, which each need to
embed a molecule, gas-optimize a charge state, and assemble the 4-point lambda_i from single
points at the two optimized geometries. Consolidated here so the DFT machinery lives in one
place (previously copy-pasted across ~5 scripts).

    lambda_i = [E_O(q_R) - E_O(q_O)] + [E_R(q_O) - E_R(q_R)]      (all in eV)

Energies are returned in eV. Density fitting (RI-J) is ON by default (~5x faster; the DF error
cancels in these same-molecule energy differences). Note: reorg.py has its OWN 4-point over the
production pipeline's cached SMD-optimized energies — this module is for on-the-fly gas-phase
recomputation (matched-protocol validation / diagnostics), not the production path.
"""
from __future__ import annotations
import numpy as np
from ase import Atoms

HARTREE_EV = 27.211386245988
BOHR = 0.52917721067


def embed(smiles):
    """RDKit ETKDG + MMFF starting geometry (ASE Atoms) for a neutral SMILES."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = AllChem.ETKDGv3(); p.randomSeed = 0xC0FFEE
    if AllChem.EmbedMolecule(m, p) != 0:
        if AllChem.EmbedMolecule(m, AllChem.ETKDGv2()) != 0:
            raise RuntimeError("embed failed")
    AllChem.MMFFOptimizeMolecule(m)
    c = m.GetConformer()
    return Atoms(symbols=[a.GetSymbol() for a in m.GetAtoms()],
                 positions=np.array([list(c.GetAtomPosition(i)) for i in range(m.GetNumAtoms())]))


def _dft(backend):
    if backend == "gpu":
        from gpu4pyscf import dft
        return dft
    from pyscf import dft
    return dft


def mol(atoms, charge, spin, basis):
    from pyscf import gto
    astr = "\n".join(f"{s} {p[0]} {p[1]} {p[2]}"
                     for s, p in zip(atoms.get_chemical_symbols(), atoms.positions))
    return gto.M(atom=astr, basis=basis, charge=charge, spin=spin, verbose=0)


def mean_field(m, xc, backend="gpu", density_fit=True, conv_tol=1e-9, max_cycle=300):
    """RKS/UKS mean-field with clean settings (RI-J by default; level_shift/damp crash gpu4pyscf)."""
    dft = _dft(backend)
    mf = (dft.RKS if m.spin == 0 else dft.UKS)(m)
    if density_fit:
        mf = mf.density_fit()
    mf.xc = xc
    mf.conv_tol = conv_tol
    mf.max_cycle = max_cycle
    return mf


def atoms_of(m):
    return Atoms(symbols=[m.atom_symbol(i) for i in range(m.natm)],
                 positions=m.atom_coords() * BOHR)


def energy(atoms, charge, spin, xc, basis, backend="gpu", density_fit=True):
    """Single-point energy (eV) of (charge, spin) at the given geometry."""
    return float(mean_field(mol(atoms, charge, spin, basis), xc, backend,
                            density_fit).kernel()) * HARTREE_EV


def gas_opt(start, charge, spin, xc, basis, backend="gpu", density_fit=True, maxsteps=100):
    """Gas-phase geometry optimization; returns (optimized Atoms, energy_eV)."""
    from pyscf.geomopt.geometric_solver import optimize
    m_opt = optimize(mean_field(mol(start, charge, spin, basis), xc, backend, density_fit),
                     maxsteps=maxsteps)
    at = atoms_of(m_opt)
    return at, energy(at, charge, spin, xc, basis, backend, density_fit)


def four_point_lambda(neu_at, E_O_at_O, ion_charge, ion_spin, xc, basis,
                      backend="gpu", density_fit=True):
    """Nelsen 4-point lambda_i (eV) for a neutral<->ion couple. neu_at/E_O_at_O are the
    optimized neutral geometry and its energy; the ion is gas-optimized here."""
    ion_at, E_R_at_R = gas_opt(neu_at, ion_charge, ion_spin, xc, basis, backend, density_fit)
    E_O_at_R = energy(ion_at, 0, 0, xc, basis, backend, density_fit)
    E_R_at_O = energy(neu_at, ion_charge, ion_spin, xc, basis, backend, density_fit)
    return (E_O_at_R - E_O_at_O) + (E_R_at_O - E_R_at_R)
