#!/usr/bin/env python
"""Reproduce the classic oligoacene INNER-SPHERE hole reorganization energies at THEIR
level of theory -- B3LYP/6-31G(d,p), gas phase -- with the standard 4-point (Nelsen) scheme.

Why this validation: the oligoacene hole lambda_i (Coropceanu/Bredas Chem. Rev. 2007) is the
cleanest possible cross-check for our reorg machinery -- rigid, closed-shell neutral <-> bound
doublet CATION (no diffuse-function / unbound-anion pathology), computed by exactly the 4-point
formula our pipeline uses. It is INDEPENDENT of D3TaLES (different functional/basis/family), so
agreement here confirms the method, not just the D3TaLES level.

Level matched to the literature:
  - functional: B3LYP (global hybrid; no range separation, no VV10 NLC).
  - basis:      6-31G(d,p)  (== 6-31G** ; the Coropceanu/Sanchez-Carrera acene level).
  - phase:      GAS. Geometries optimized at this level; single points at this level.
  - couple:     HOLE only  = neutral singlet <-> +1 cation doublet.
                (electron/anion is skipped: 6-31G(d,p) has no diffuse functions, so acene
                 anions are ill-conditioned -- the exact pathology we already flag elsewhere.)

RI-J density fitting + gpu4pyscf for speed (DF error cancels in the 4-point energy differences;
verified elsewhere to leave lambda unchanged to <1 meV). Writes one JSON per molecule under
results/reorg_anchors/inner_acene/calc/. Resumable for matching solver-qualified records; use a fresh --tag for historical outputs.

  python -m redox.workflows.gpu --gpu <idx> --out logs/runs/acene_reference --threads 2 -- python scripts/validation/reorg_literature/validate_reorg_worker_inner_acene.py \
      --source results/reorg_anchors/oligoacene_inner.csv --backend gpu
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
BASE = ROOT / "results" / "reorg_anchors" / "inner_acene"
# BASIS / CALC are set in main() from --basis/--tag so we can run 6-31G** and 6-311G**
# side-by-side without clobbering each other.
BASIS = "6-31g(d,p)"       # default == 6-31G**
CALC = BASE / "calc"
XC = "b3lyp"
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


def _mf(mol, backend):
    if backend == "gpu":
        from gpu4pyscf import dft as gdft
        KS = gdft.RKS if mol.spin == 0 else gdft.UKS
    else:
        from pyscf import dft as cdft
        KS = cdft.RKS if mol.spin == 0 else cdft.UKS
    mf = KS(mol).density_fit()
    mf.xc = XC
    mf.conv_tol = 1e-9
    mf.max_cycle = 200
    return mf


def _opt(atoms, charge, spin, backend):
    molopt = optimized_molecule(_mf(_mol(atoms, charge, spin), backend), maxsteps=100)
    BOHR = 0.52917721067
    return [(molopt.atom_symbol(i), tuple(molopt.atom_coords()[i] * BOHR))
            for i in range(molopt.natm)]


def _E(atoms, charge, spin, backend):
    return energy_Ha(_mf(_mol(atoms, charge, spin), backend))


def _hole_lambda(g_neu, E_neu_neu, backend, *, return_evidence=False):
    """4-point inner-sphere hole lambda: neutral singlet <-> +1 cation doublet. eV."""
    g_cat = _opt(g_neu, 1, 1, backend)                 # optimize cation (doublet)
    E_cat_cat = _E(g_cat, 1, 1, backend)               # cation @ cation geom
    E_cat_at_neu = _E(g_neu, 1, 1, backend)            # cation @ neutral geom
    E_neu_at_cat = _E(g_cat, 0, 0, backend)            # neutral @ cation geom
    relax_cat = (E_cat_at_neu - E_cat_cat) * HARTREE_EV
    relax_neu = (E_neu_at_cat - E_neu_neu) * HARTREE_EV
    lam = relax_cat + relax_neu
    values = (round(lam, 4), round(relax_cat, 4), round(relax_neu, 4))
    evidence = dict(energies_Ha=dict(neutral_at_neutral=E_neu_neu, ion_at_ion=E_cat_cat,
                                   ion_at_neutral=E_cat_at_neu, neutral_at_ion=E_neu_at_cat),
                    neutral_geometry_A=g_neu, ion_geometry_A=g_cat,
                    ion_charge=1, ion_spin=1, lambda_unrounded_eV=lam)
    return (*values, evidence) if return_evidence else values


def process(row, backend):
    with record_guard(CALC / f"{row['id']}.json"):
        return _process(row, backend)


def _process(row, backend):
    gid = row["id"]
    out = CALC / f"{gid}.json"
    identity = dict(id=str(gid),
                    smiles=str(row["smiles"]),
                    couple="hole",
                    xc=XC, basis=BASIS,
                    phase="gas", density_fit=True, conv_tol=1e-9, max_cycle=200,
                    reference_eV=to_float(row.get("lit_hole_eV")))
    if existing_result(out, identity, __file__, required_fields=("our_hole", "relax_cat", "relax_neu")) is not None:
        return "skip"
    rec = dict(id=gid, smiles=row["smiles"], lit_hole_eV=to_float(row.get("lit_hole_eV")),
               protocol=f"{XC}/{BASIS} gas 4-point (hole, neutral<->cation)",
               our_hole=None, relax_cat=None, relax_neu=None, status="ok", error="")
    try:
        atoms0 = _embed(row["smiles"])
        g_neu = _opt(atoms0, 0, 0, backend)
        E_neu_neu = _E(g_neu, 0, 0, backend)
        lam, rc, rn, points = _hole_lambda(g_neu, E_neu_neu, backend, return_evidence=True)
        rec["four_point_evidence"] = points
        rec["our_hole"], rec["relax_cat"], rec["relax_neu"] = lam, rc, rn
    except Exception as e:
        rec["status"] = "error"; rec["error"] = f"{type(e).__name__}: {e}"
        rec["trace"] = traceback.format_exc()[-800:]
    CALC.mkdir(parents=True, exist_ok=True)
    if rec["status"] == "ok":
        rec["solver_verification"] = verification(identity, __file__)
    atomic_write_json(out, rec)
    return rec["status"]


def main():
    global BASIS, CALC
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=str(ROOT / "results/reorg_anchors/oligoacene_inner.csv"))
    ap.add_argument("--only", default=None)
    ap.add_argument("--backend", default="gpu", choices=["cpu", "gpu"])
    ap.add_argument("--basis", default="6-31g(d,p)", help="pyscf basis, e.g. 6-31g(d,p) or 6-311g(d,p)")
    ap.add_argument("--tag", default=None, help="output subdir under inner_acene/ (default derived from basis)")
    a = ap.parse_args()
    BASIS = a.basis
    tag = a.tag or ("b3lyp_" + BASIS.replace("(", "").replace(")", "").replace("*", "s").replace(",", ""))
    CALC = BASE / tag / "calc"
    print(f"[cfg] XC={XC}  basis={BASIS}  out={CALC}", flush=True)
    df = pd.read_csv(a.source)
    if a.only:
        df = df[df["id"] == a.only].reset_index(drop=True)
    problems = 0
    for _, row in df.iterrows():
        st = process(row, a.backend)
        problems += st not in ("ok", "skip")
        r = json.loads((CALC / f"{row['id']}.json").read_text())
        lit = r.get("lit_hole_eV")
        diff = ""
        try:
            diff = f"{abs(float(r['our_hole']) - float(lit)):.3f}"
        except Exception:
            pass
        print(f"[{row['id']:11s}] {st:5s}  our_hole={r.get('our_hole')} eV  "
              f"lit={lit} eV  |diff|={diff}  {r.get('error','')}", flush=True)

    if problems:
        raise SystemExit(f"{problems} reference records failed or remain unqualified")


if __name__ == "__main__":
    main()
