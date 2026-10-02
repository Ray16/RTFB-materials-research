#!/usr/bin/env python
"""Resumable worker: recompute D3TaLES reorganization energies at THEIR TRUE level of theory —
IP-tuned LC-wHPBE / Def2SVP, gas phase — so the comparison is genuinely apples-to-apples.

Match details (from Duke et al. 2023 + the D3TaLES README):
  - functional: LC-wHPBE  ==  libxc `lc_wpbe` (fully long-range-corrected wPBE, HJS hole:
    0% short-range / 100% long-range exact exchange).
  - omega:      D3TaLES's OWN per-molecule IP-tuned value (the `omega` column of the dump),
    injected via mf.omega — so we do NOT re-derive the tuning, we reuse theirs.
  - basis:      def2-svp (Def2SVP, no diffuse — exactly as they ran it).
  - phase:      GAS (no solvent), geometries optimized at this level.
  - reorg:      standard 4-point (Nelsen) for BOTH couples they report:
                  hole     = neutral <-> +1 cation
                  electron = neutral <-> -1 anion

RI-J density fitting + gpu4pyscf (the range-separated exchange is far too slow otherwise).
Writes one JSON per molecule under results/d3tales_reorg_validation_d3level/calc/. Resumable.

  CUDA_VISIBLE_DEVICES=<i> PYTHONPATH=src python scripts/validation/reorg_d3tales/validate_reorg_worker_d3tales.py \
      --source results/d3tales_reorg_validation_d3level/molecules.csv --only 80JNKV
"""
from __future__ import annotations
import argparse
import json
import traceback
from pathlib import Path

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
BASE = ROOT / "results" / "d3tales_reorg_validation_d3level"
CALC = BASE / "calc"
GEO = BASE / "geom"
BASIS = "def2-svp"
XC = "lc_wpbe"                      # = Gaussian LC-wHPBE
HARTREE_EV = 27.211386245988


def _embed(smiles):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    p = AllChem.ETKDGv3(); p.randomSeed = 0xC0FFEE
    if AllChem.EmbedMolecule(m, p) != 0:
        p.useRandomCoords = True
        AllChem.EmbedMolecule(m, p)
    AllChem.MMFFOptimizeMolecule(m)
    conf = m.GetConformer()
    return [(a.GetSymbol(), tuple(conf.GetAtomPosition(a.GetIdx()))) for a in m.GetAtoms()]


def _mol(atoms, charge, spin):
    from pyscf import gto
    astr = "\n".join(f"{s} {p[0]} {p[1]} {p[2]}" for s, p in atoms)
    return gto.M(atom=astr, basis=BASIS, charge=charge, spin=spin, verbose=0)


def _mf(mol, omega, backend):
    if backend == "gpu":
        from gpu4pyscf import dft as gdft
        KS = gdft.RKS if mol.spin == 0 else gdft.UKS
    else:
        from pyscf import dft as cdft
        KS = cdft.RKS if mol.spin == 0 else cdft.UKS
    mf = KS(mol).density_fit()
    mf.xc = XC
    mf.omega = float(omega)        # inject D3TaLES's tuned omega for THIS molecule
    mf.conv_tol = 1e-9
    mf.max_cycle = 200
    return mf


def _opt(atoms, charge, spin, omega, backend):
    from pyscf.geomopt.geometric_solver import optimize
    molopt = optimize(_mf(_mol(atoms, charge, spin), omega, backend), maxsteps=100)
    BOHR = 0.52917721067
    optatoms = [(molopt.atom_symbol(i), tuple(molopt.atom_coords()[i] * BOHR))
                for i in range(molopt.natm)]
    return optatoms


def _E(atoms, charge, spin, omega, backend):
    return float(_mf(_mol(atoms, charge, spin), omega, backend).kernel())


def _couple(g_neu, E_neu_neu, q, spin, omega, backend):
    """4-point reorg for neutral<->(q,spin) ion (q=+1 hole, q=-1 electron). eV."""
    g_ion = _opt(g_neu, q, spin, omega, backend)
    E_ion_ion = _E(g_ion, q, spin, omega, backend)
    E_ion_at_neu = _E(g_neu, q, spin, omega, backend)     # ion charge, neutral geom
    E_neu_at_ion = _E(g_ion, 0, 0, omega, backend)        # neutral, ion geom
    lam = ((E_ion_at_neu - E_ion_ion) + (E_neu_at_ion - E_neu_neu)) * HARTREE_EV
    return round(lam, 4)


def process(row, backend):
    gid = row["id"]
    out = CALC / f"{gid}.json"
    if out.exists():
        try:
            if json.loads(out.read_text()).get("status") in ("ok", "partial"):
                return "skip"
        except Exception:
            pass
    rec = dict(id=gid, smiles=row["smiles"], omega=float(row["omega"]),
               d3_hole=row.get("d3_hole"), d3_electron=row.get("d3_electron"),
               n_atoms=row.get("n_atoms"),
               protocol="D3TaLES-matched: lc_wpbe(omega_tuned)/def2-svp gas 4-point",
               our_hole=None, our_electron=None, status="ok", error="")
    GEO.mkdir(parents=True, exist_ok=True)
    try:
        om = float(row["omega"])
        atoms0 = _embed(row["smiles"])
        g_neu = _opt(atoms0, 0, 0, om, backend)
        E_neu_neu = _E(g_neu, 0, 0, om, backend)
        for label, q, spin, key in (("hole", 1, 1, "our_hole"), ("electron", -1, 1, "our_electron")):
            d3 = rec.get(f"d3_{label}")
            if d3 is None or (isinstance(d3, float) and d3 != d3):
                continue
            try:
                rec[key] = _couple(g_neu, E_neu_neu, q, spin, om, backend)
            except Exception as e:
                rec["status"] = "partial"; rec["error"] += f"{label}:{type(e).__name__} "
    except Exception as e:
        rec["status"] = "error"; rec["error"] = f"{type(e).__name__}: {e}"
        rec["trace"] = traceback.format_exc()[-800:]
    CALC.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=1))
    return rec["status"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=str(BASE / "molecules.csv"))
    ap.add_argument("--only", default=None)
    ap.add_argument("--offset", type=int, default=0)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--backend", default="gpu")
    a = ap.parse_args()
    df = pd.read_csv(a.source)
    if a.only:
        df = df[df["id"] == a.only].reset_index(drop=True)
    else:
        df = df.iloc[a.offset::a.stride].reset_index(drop=True)
    if a.limit:
        df = df.head(a.limit)
    for _, row in df.iterrows():
        st = process(row, a.backend)
        r = json.loads((CALC / f"{row['id']}.json").read_text()) if (CALC / f"{row['id']}.json").exists() else {}
        print(f"[{row['id']}] {st}  hole={r.get('our_hole')}(d3 {r.get('d3_hole')})  "
              f"elec={r.get('our_electron')}(d3 {r.get('d3_electron')})", flush=True)
        if st in ("ok", "partial", "skip"):
            print("DONE_MARKER", flush=True)


if __name__ == "__main__":
    main()
