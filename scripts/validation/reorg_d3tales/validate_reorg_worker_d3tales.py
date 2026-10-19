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
Writes one JSON per molecule under results/d3tales_reorg_validation_d3level/calc/. Resumable only for matching solver-qualified records; historical files are preserved. Choose --out for a fresh run.

  python -m redox.workflows.gpu --gpu <i> --out logs/runs/d3level -- python scripts/validation/reorg_d3tales/validate_reorg_worker_d3tales.py \
      --source results/d3tales_reorg_validation_d3level/molecules.csv --only 80JNKV
"""
from __future__ import annotations
import argparse
import json
import traceback
from pathlib import Path

from redox.core.io import atomic_write_json
from redox.core.common import to_float
from redox.qm.gas_reference import (energy_Ha, optimized_molecule,
                                     verification, existing_result, record_guard)

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
    molopt = optimized_molecule(_mf(_mol(atoms, charge, spin), omega, backend), maxsteps=100)
    BOHR = 0.52917721067
    optatoms = [(molopt.atom_symbol(i), tuple(molopt.atom_coords()[i] * BOHR))
                for i in range(molopt.natm)]
    return optatoms


def _E(atoms, charge, spin, omega, backend):
    return energy_Ha(_mf(_mol(atoms, charge, spin), omega, backend))


def _couple(g_neu, E_neu_neu, q, spin, omega, backend, *, return_evidence=False):
    """4-point reorg for neutral<->(q,spin) ion (q=+1 hole, q=-1 electron). eV."""
    g_ion = _opt(g_neu, q, spin, omega, backend)
    E_ion_ion = _E(g_ion, q, spin, omega, backend)
    E_ion_at_neu = _E(g_neu, q, spin, omega, backend)     # ion charge, neutral geom
    E_neu_at_ion = _E(g_ion, 0, 0, omega, backend)        # neutral, ion geom
    lam = ((E_ion_at_neu - E_ion_ion) + (E_neu_at_ion - E_neu_neu)) * HARTREE_EV
    evidence = dict(energies_Ha=dict(neutral_at_neutral=E_neu_neu, ion_at_ion=E_ion_ion,
                                   ion_at_neutral=E_ion_at_neu, neutral_at_ion=E_neu_at_ion),
                    neutral_geometry_A=g_neu, ion_geometry_A=g_ion,
                    ion_charge=q, ion_spin=spin, lambda_unrounded_eV=lam)
    return (round(lam, 4), evidence) if return_evidence else round(lam, 4)


def process(row, backend):
    with record_guard(CALC / f"{row['id']}.json"):
        return _process(row, backend)


def _process(row, backend):
    gid = row["id"]
    out = CALC / f"{gid}.json"
    omega = to_float(row.get("omega"))
    if omega is None:
        raise ValueError("Finite tuned omega is required; no calculation started")
    identity = dict(id=str(gid),
                    smiles=str(row["smiles"]),
                    xc=XC, basis=BASIS, omega=omega,
                    phase="gas", density_fit=True, conv_tol=1e-9, max_cycle=200,
                    hole_reference_eV=to_float(row.get("d3_hole")),
                    electron_reference_eV=to_float(row.get("d3_electron")))
    required = tuple(name for name, ref in (("our_hole", identity["hole_reference_eV"]),
                                             ("our_electron", identity["electron_reference_eV"]))
                     if ref is not None)
    if existing_result(out, identity, __file__, required_fields=required) is not None:
        return "skip"
    rec = dict(id=gid, smiles=row["smiles"], omega=omega,
               d3_hole=to_float(row.get("d3_hole")), d3_electron=to_float(row.get("d3_electron")),
               n_atoms=to_float(row.get("n_atoms")),
               protocol=f"D3TaLES-matched: {XC}(omega_tuned)/{BASIS} gas 4-point",
               our_hole=None, our_electron=None, status="ok", error="")
    if rec["d3_hole"] is None and rec["d3_electron"] is None:
        raise ValueError("No reference couples are requested; no calculation started")
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
                rec[key], points = _couple(g_neu, E_neu_neu, q, spin, om, backend, return_evidence=True)
                rec.setdefault("four_point_evidence", {})[label] = points
                rec.setdefault("completed_couple_solver_verification", {})[label] = verification(
                    dict(identity, couple=label), __file__)
            except Exception as e:
                rec["status"] = "partial"; rec["error"] += f"{label}:{type(e).__name__} "
    except Exception as e:
        rec["status"] = "error"; rec["error"] = f"{type(e).__name__}: {e}"
        rec["trace"] = traceback.format_exc()[-800:]
    CALC.mkdir(parents=True, exist_ok=True)
    if rec["status"] == "ok":
        rec["solver_verification"] = verification(identity, __file__)
    atomic_write_json(out, rec)
    return rec["status"]


def main():
    global BASE, CALC, GEO
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=str(BASE / "molecules.csv"))
    ap.add_argument("--only", default=None)
    ap.add_argument("--offset", type=int, default=0)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--backend", default="gpu", choices=["cpu", "gpu"])
    ap.add_argument("--out", type=Path, help="fresh output base; existing unqualified records are preserved")
    a = ap.parse_args()
    if a.out is not None:
        BASE = a.out.resolve(); CALC = BASE / "calc"; GEO = BASE / "geom"
    df = pd.read_csv(a.source)
    if a.only:
        df = df[df["id"] == a.only].reset_index(drop=True)
    else:
        df = df.iloc[a.offset::a.stride].reset_index(drop=True)
    if a.limit:
        df = df.head(a.limit)
    problems = 0
    for _, row in df.iterrows():
        st = process(row, a.backend)
        problems += st not in ("ok", "skip")
        r = json.loads((CALC / f"{row['id']}.json").read_text()) if (CALC / f"{row['id']}.json").exists() else {}
        print(f"[{row['id']}] {st}  hole={r.get('our_hole')}(d3 {r.get('d3_hole')})  "
              f"elec={r.get('our_electron')}(d3 {r.get('d3_electron')})", flush=True)
        if st in ("ok", "skip"):
            print("DONE_MARKER", flush=True)

    if problems:
        raise SystemExit(f"{problems} reference records failed or remain unqualified")


if __name__ == "__main__":
    main()
