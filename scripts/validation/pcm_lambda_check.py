#!/usr/bin/env python
"""DFT-based outer-sphere validation: molecular-cavity PCM lambda_o for the lead candidates,
vs the shape-blind Born lambda_o. First runs the built-in sphere sanity check (PCM path must
reproduce analytic Born). Matches the existing pool PCM protocol (b3lyp/def2-svpd). GPU via
gpu_reserve; vertical (two SCFs at one geometry), no geometry opt, local-scratch only."""
import sys
sys.path.insert(0, "/nfs/lambda_stor_01/homes/rzhu/0_redox_refactor/src")
import numpy as np
from ase.io import read as aseread
from pyscf import gto
from redox.backends.dft import _build_smd_mf
from redox.properties.lambda_outer import (pcm_lambda_outer, sphere_sanity_check,
                                           born_lambda_outer, sasa_radius_A, solvent_constants)
from rdkit import Chem

sc = solvent_constants(); eps_s, eps_op = sc["eps_s"], sc["eps_op"]
print("sphere sanity check (PCM path must ~= analytic Born):", sphere_sanity_check(eps_s, eps_op), flush=True)

ARCH = "/nfs/lambda_stor_01/homes/rzhu/archive/0_redox/calcs/dft"
XC, BAS = "b3lyp", "def2-svpd"   # match the existing lambda_o_pcm protocol

def dm_at(geom_xyz, charge, mult):
    at = aseread(geom_xyz); syms = at.get_chemical_symbols()
    atom_str = "\n".join(f"{s} {p[0]} {p[1]} {p[2]}" for s, p in zip(syms, at.positions))
    mol = gto.M(atom=atom_str, basis=BAS, charge=charge, spin=mult - 1, verbose=0)
    mf = _build_smd_mf(mol, XC, "acetonitrile", "gpu", None, None)
    mf.kernel()
    dm = mf.make_rdm1()
    dm = dm.get() if hasattr(dm, "get") else np.asarray(dm)   # cupy -> numpy for gpu4pyscf
    return mol, dm

# lead catholyte: N'-methyl-DHP, couple neu/ox1 (vertical at the neutral geometry)
CASES = [("exp_dmpz_mb_methyl", "N'-methyl-DHP (catholyte)", "neu", (0, 1), (1, 2), 0.700)]
for gid, desc, geomstate, (qr, mr), (qo, mo), born in CASES:
    geom = f"{ARCH}/{gid}/{geomstate}/opt.xyz"
    mol_r, dm_r = dm_at(geom, qr, mr)
    _, dm_o = dm_at(geom, qo, mo)
    lam_pcm = pcm_lambda_outer(mol_r, dm_r, dm_o, eps_s, eps_op)
    print(f"\n{desc}: Born lambda_o = {born:.3f} eV | PCM lambda_o = {lam_pcm:.3f} eV | "
          f"PCM/Born = {lam_pcm/born:.3f}", flush=True)
