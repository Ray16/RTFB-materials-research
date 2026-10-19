#!/usr/bin/env python
"""Generic inner-sphere reorganization-energy validation worker (4-point Nelsen), for either
couple, at an arbitrary functional/basis, gas phase. Used to reproduce literature lambda values
at THEIR level of theory.

  couple = hole      : neutral singlet  <-> +1 cation doublet
  couple = electron  : neutral singlet  <-> -1 anion  doublet   (BOUND-ANION guarded)

For electron: a valid 4-point lambda needs a BOUND radical anion. We record the vertical-anion
HOMO (the added electron's orbital, at the neutral geometry); if it is > 0 the anion is unbound
and the lambda is ill-defined -> flagged `anion_unbound` (do not trust the number). This is the
same guard our production pipeline uses.

RI-J density fitting + gpu4pyscf (DF error cancels in the 4-point energy differences).
Resumable for matching solver-qualified records; choose a fresh --tag for historical outputs.
One JSON per molecule under <base>/<tag>/calc/.

  python -m redox.workflows.gpu --gpu <i> --out logs/runs/gas_reference --threads 2 -- python \
     scripts/validation/reorg_literature/validate_reorg_worker_generic.py --source <csv> --couple electron \
     --basis 6-31g(d,p) --tag b3lyp_631gdp_electron --backend gpu
"""
from __future__ import annotations
import argparse
import json
import traceback
from pathlib import Path

from redox.core.io import atomic_write_json
from redox.core.common import to_float
from redox.qm.gas_reference import (energy_Ha, optimized_molecule, homo_eV,
                                     verification, existing_result, record_guard)

import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

ROOT = Path(__file__).resolve().parents[3]
HARTREE_EV = 27.211386245988
XC = "b3lyp"
BASIS = "6-31g(d,p)"
CALC = None


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


def _mf(mol, backend):
    if backend == "gpu":
        from gpu4pyscf import dft as gdft
        KS = gdft.RKS if mol.spin == 0 else gdft.UKS
    else:
        from pyscf import dft as cdft
        KS = cdft.RKS if mol.spin == 0 else cdft.UKS
    mf = KS(mol).density_fit()
    mf.xc = XC; mf.conv_tol = 1e-9; mf.max_cycle = 200
    return mf


def _homo_eV(mf):
    return homo_eV(mf)


def _E(atoms, charge, spin, backend, want_homo=False):
    mf = _mf(_mol(atoms, charge, spin), backend)
    e = energy_Ha(mf)
    if want_homo:
        try:
            mf.mo_energy  # ensure populated
            h = _homo_eV(mf)
        except Exception:
            h = None
        return e, h
    return e


def _opt(atoms, charge, spin, backend):
    m = optimized_molecule(_mf(_mol(atoms, charge, spin), backend), maxsteps=100)
    B = 0.52917721067
    return [(m.atom_symbol(i), tuple(m.atom_coords()[i] * B)) for i in range(m.natm)]


def _four_point(g_neu, E_neu_neu, q, spin, backend, *, return_evidence=False):
    """4-point lambda for neutral(0,0) <-> ion(q,spin). Returns (lam, relax_ion, relax_neu, ion_homo)."""
    g_ion = _opt(g_neu, q, spin, backend)
    E_ion_ion = _E(g_ion, q, spin, backend)
    E_ion_at_neu, ion_homo = _E(g_neu, q, spin, backend, want_homo=True)   # vertical ion @ neutral geom
    E_neu_at_ion = _E(g_ion, 0, 0, backend)
    relax_ion = (E_ion_at_neu - E_ion_ion) * HARTREE_EV
    relax_neu = (E_neu_at_ion - E_neu_neu) * HARTREE_EV
    values = (round(relax_ion + relax_neu, 4), round(relax_ion, 4), round(relax_neu, 4), ion_homo)
    evidence = dict(energies_Ha=dict(neutral_at_neutral=E_neu_neu, ion_at_ion=E_ion_ion,
                                   ion_at_neutral=E_ion_at_neu, neutral_at_ion=E_neu_at_ion),
                    neutral_geometry_A=g_neu, ion_geometry_A=g_ion,
                    ion_charge=q, ion_spin=spin, lambda_unrounded_eV=relax_ion+relax_neu)
    return (*values, evidence) if return_evidence else values


def process(row, couple, backend):
    with record_guard(CALC / f"{row['id']}.json"):
        return _process(row, couple, backend)


def _process(row, couple, backend):
    gid = row["id"]
    out = CALC / f"{gid}.json"
    identity = dict(id=str(gid),
                    smiles=str(row["smiles"]),
                    couple=couple,
                    xc=XC, basis=BASIS,
                    phase="gas", density_fit=True, conv_tol=1e-9, max_cycle=200,
                    reference_eV=to_float(row.get(f"lit_{couple}_eV", row.get("lit_eV"))))
    if existing_result(out, identity, __file__, required_fields=("our_lambda", "relax_ion", "relax_neu")) is not None:
        return "skip"
    q, spin = (1, 1) if couple == "hole" else (-1, 1)
    rec = dict(id=gid, smiles=row["smiles"], couple=couple,
               lit_eV=to_float(row.get(f"lit_{couple}_eV", row.get("lit_eV"))),
               protocol=f"{XC}/{BASIS} gas 4-point ({couple})",
               our_lambda=None, relax_ion=None, relax_neu=None,
               ion_homo_eV=None, anion_unbound=None, status="ok", error="")
    try:
        atoms0 = _embed(row["smiles"])
        g_neu = _opt(atoms0, 0, 0, backend)
        E_neu_neu = _E(g_neu, 0, 0, backend)
        lam, ri, rn, ion_homo, points = _four_point(g_neu, E_neu_neu, q, spin, backend, return_evidence=True)
        rec["four_point_evidence"] = points
        rec.update(our_lambda=lam, relax_ion=ri, relax_neu=rn,
                   ion_homo_eV=ion_homo)
        if couple == "electron":
            rec["anion_unbound"] = None if ion_homo is None else bool(ion_homo > 0)
            if ion_homo is None or rec["anion_unbound"]:
                rec["diagnostic_lambda_eV"] = rec["our_lambda"]
                rec["our_lambda"] = None
                rec["status"] = "binding_unverified" if ion_homo is None else "anion_unbound"
                rec["error"] = "Radical-anion binding is unavailable or failed"
    except Exception as e:
        rec["status"] = "error"; rec["error"] = f"{type(e).__name__}: {e}"
        rec["trace"] = traceback.format_exc()[-800:]
    CALC.mkdir(parents=True, exist_ok=True)
    if rec["status"] == "ok":
        rec["solver_verification"] = verification(identity, __file__)
    atomic_write_json(out, rec)
    return rec["status"]


def main():
    global XC, BASIS, CALC
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True)
    ap.add_argument("--couple", choices=["hole", "electron"], required=True)
    ap.add_argument("--xc", default="b3lyp")
    ap.add_argument("--basis", default="6-31g(d,p)")
    ap.add_argument("--tag", required=True, help="output subdir under results/reorg_anchors/electron/")
    ap.add_argument("--only", default=None)
    ap.add_argument("--backend", default="gpu", choices=["cpu", "gpu"])
    a = ap.parse_args()
    XC, BASIS = a.xc, a.basis
    CALC = ROOT / "results" / "reorg_anchors" / "electron" / a.tag / "calc"
    print(f"[cfg] XC={XC} basis={BASIS} couple={a.couple} out={CALC}", flush=True)
    df = pd.read_csv(a.source)
    if a.only:
        df = df[df["id"] == a.only].reset_index(drop=True)
    problems = 0
    for _, row in df.iterrows():
        st = process(row, a.couple, a.backend)
        problems += st not in ("ok", "skip")
        r = json.loads((CALC / f"{row['id']}.json").read_text())
        ub = " UNBOUND-ANION!" if r.get("anion_unbound") else ""
        print(f"[{row['id']:16s}] {st:5s} our={r.get('our_lambda')} eV  lit={r.get('lit_eV')}  "
              f"ion_homo={r.get('ion_homo_eV')}eV{ub}  {r.get('error','')}", flush=True)

    if problems:
        raise SystemExit(f"{problems} reference records failed or remain unqualified")


if __name__ == "__main__":
    main()
