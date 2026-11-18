#!/usr/bin/env python
"""Per-couple catholyte anchor: compute the EXACT BMEPZ parent (5,10-bis(2-methoxyethyl)-
5,10-dihydrophenazine) at the same DFT+SMD / legacy-xTB-RRHO level as the pool, and compare
its two oxidation potentials to the measured -0.29 / +0.50 V vs Fc (ScienceDirect S2451929419303158).
This turns the catholyte bias from 'corroborated' into a measured per-couple error for the
EXACT same chemistry as the proposed DHP catholytes.

Writes only a tiny optimized-geometry xyz to LOCAL scratch (TMPDIR) — negligible NFS footprint.
Run on a GPU via gpu_reserve.  Measured couples from the BMEPZ parent paper.
"""
import os, sys, tempfile
from pathlib import Path
os.environ.setdefault('XTB_BIN', '/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/xtb')
from ase.io import read as aseread, write as asewrite
from rdkit import Chem
from rdkit.Chem import AllChem
from redox.backends.dft import dft_smd, _thermal_correction
from redox.core.common import read_result, DFT

SCRATCH = Path(os.environ.get("TMPDIR", "/tmp")) / "bmepz_anchor"
SCRATCH.mkdir(parents=True, exist_ok=True)
BMEPZ = "COCCN1c2ccccc2N(CCOC)c2ccccc21"   # 5,10-bis(2-methoxyethyl)-5,10-dihydrophenazine
STATES = [("neu", 0, 1), ("ox1", 1, 2), ("ox2", 2, 1)]
MEAS = {"E1": -0.29, "E2": 0.50}  # V vs Fc, MeCN

def start_xyz():
    m = Chem.AddHs(Chem.MolFromSmiles(BMEPZ))
    AllChem.EmbedMolecule(m, randomSeed=0xC0FFEE); AllChem.MMFFOptimizeMolecule(m)
    p = SCRATCH / "bmepz_start.xyz"
    conf = m.GetConformer()
    syms = [a.GetSymbol() for a in m.GetAtoms()]
    with open(p, "w") as f:
        f.write(f"{m.GetNumAtoms()}\n\n")
        for s, i in zip(syms, range(m.GetNumAtoms())):
            pos = conf.GetAtomPosition(i); f.write(f"{s} {pos.x} {pos.y} {pos.z}\n")
    return p

def fc_ref():
    fg = {}
    for s in ("neu", "ox"):
        r = read_result("ferrocene", s)
        at = aseread(str(DFT / "ferrocene" / s / "opt.xyz"))
        th = _thermal_correction(at, int(r["charge"]), int(r["mult"]))["g_thermal_eV"]
        fg[s] = r["e_smd_eV"] + th
    return fg["ox"] - fg["neu"]

def main():
    start = start_xyz()
    fcref = fc_ref()
    print(f"ferrocene legacy fcref = {fcref:.4f} eV", flush=True)
    G = {}
    for name, q, mult in STATES:
        opt_out = SCRATCH / f"{name}_opt.xyz"
        res = dft_smd(start, q, mult, backend="gpu", do_opt=True, do_freq=False, opt_out=opt_out)
        e_smd = res["e_smd_eV"]
        at = aseread(str(opt_out))
        th = _thermal_correction(at, q, mult)["g_thermal_eV"]
        G[name] = e_smd + th
        print(f"[{name}] q={q} mult={mult}  e_smd={e_smd:.3f}  g_thermal={th:.3f}  G={G[name]:.3f} eV", flush=True)
    E1 = (G["ox1"] - G["neu"]) - fcref
    E2 = (G["ox2"] - G["ox1"]) - fcref
    print("\n=== BMEPZ parent, computed vs measured (V vs Fc) ===")
    print(f"  E1: calc {E1:+.3f}  meas {MEAS['E1']:+.3f}  err {E1-MEAS['E1']:+.3f}")
    print(f"  E2: calc {E2:+.3f}  meas {MEAS['E2']:+.3f}  err {E2-MEAS['E2']:+.3f}")
    mae = (abs(E1-MEAS['E1']) + abs(E2-MEAS['E2'])) / 2
    bias = ((E1-MEAS['E1']) + (E2-MEAS['E2'])) / 2
    print(f"  per-couple catholyte MAE = {mae:.3f} V   bias = {bias:+.3f} V")

if __name__ == "__main__":
    main()
