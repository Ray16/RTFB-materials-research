#!/usr/bin/env python
"""Conformer-matched inner-sphere lambda for the merrifield candidates (generalizes
scripts/validation/reorg_d3tales/reorg_recompute_guarded.py from the single 0/-1 D3TaLES couple to arbitrary
multi-charge reduction/oxidation ladders).

Why: for the GRAFTED models the flexible 4-methylbenzyl tether rotates between charge states,
so the baseline 4-point lambda charges a whole tether isomerization to lambda (neutral<->ion
heavy-atom RMSD of 1-2.4 A; the standalone _sa refs, methyl-capped, stay <0.05 A). That soft
torsional mode is (a) not part of the fast inner-sphere ET and (b) physically frozen in the real
polymer (tether anchored to the polystyrene backbone). So the conformer-matched value is both
cleaner AND more representative of the grafted polymer.

Recipe (per molecule):
  reference conformer = the NEUTRAL 'neu' optimized geometry.
  For every non-neu state S: re-optimize S starting FROM the neu geometry with all rotatable-bond
  dihedrals FROZEN (geomeTRIC $freeze), giving S_cf in the neu conformer. neu itself is the ref.
  Then recompute each adjacent-charge couple's 4-point lambda using the cf geometries, in GAS and
  in SMD (SMD also rescues an unbound gas anion). Report raw vs cf vs cf+smd + RMSDs; pick best
  (cf+smd if available, else cf-gas, else raw).

  gpu_reserve run <idx> -- env OMP_NUM_THREADS=2 PYTHONPATH=src \
      python scripts/analysis/reorg_conformer_matched.py --only aq_benzylamino --backend gpu
"""
from __future__ import annotations
import argparse, json, traceback
from pathlib import Path
import numpy as np
from ase.io import read as ase_read
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

import redox.qm.dft as D
from redox.properties.reorg import _couples
from redox.core.common import DFT, read_manifest, read_result as _res

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "reorg_cf"
RMSD_MAX = 0.40   # heavy-atom RMSD (A) above which the baseline couple is conformer-contaminated

_CF_CONV = dict(convergence_energy=5e-6, convergence_grms=4e-4, convergence_gmax=1e-3,
                convergence_drms=2e-3, convergence_dmax=4e-3)


def heavy_rmsd(p1, p2):
    a1 = ase_read(str(p1)); a2 = ase_read(str(p2))
    m = np.array(a1.get_chemical_symbols()) != "H"
    A = a1.positions[m] - a1.positions[m].mean(0)
    B = a2.positions[m] - a2.positions[m].mean(0)
    H = A.T @ B; U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T)); R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return float(np.sqrt(((A @ R.T - B) ** 2).sum(1).mean()))


def rotatable_dihedral_constraints(smiles, out_path):
    """geomeTRIC $freeze locking every rotatable-bond dihedral (1-indexed, AddHs order == xyz)."""
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    patt = Chem.MolFromSmarts("[!$(*#*)&!D1]-!@[!$(*#*)&!D1]")
    lines, seen = ["$freeze"], set()
    for a, b in m.GetSubstructMatches(patt):
        na = [n.GetIdx() for n in m.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() != b]
        nb = [n.GetIdx() for n in m.GetAtomWithIdx(b).GetNeighbors() if n.GetIdx() != a]
        if not na or not nb:
            continue
        pick = lambda c: ([i for i in c if m.GetAtomWithIdx(i).GetAtomicNum() > 1] or c)[0]
        quad = (pick(na), a, b, pick(nb))
        if quad in seen:
            continue
        seen.add(quad)
        lines.append(f"dihedral {quad[0]+1} {quad[1]+1} {quad[2]+1} {quad[3]+1}")
    out_path.write_text("\n".join(lines) + "\n")
    return len(lines) - 1


# energy single points, memoized on (charge, mult, geom-path)
_CACHE = {}
def _e(kind, geom, q, m, bk):
    key = (kind, str(geom), q, m)
    if key in _CACHE:
        return _CACHE[key]
    if kind == "gas":
        r = D.dft_smd(geom, q, m, do_opt=False, do_gas=True, do_smd=False, do_freq=False, backend=bk)
        v = r.get("e_gas_eV"); diag = (r.get("gas_homo_eV"), bool(r.get("anion_unbound")))
    else:
        r = D.dft_smd(geom, q, m, do_opt=False, do_gas=False, do_smd=True, do_freq=False, backend=bk)
        v = r.get("e_smd_eV"); diag = (None, None)
    _CACHE[key] = (v, diag)
    return v, diag


def _smiles_for(gid):
    for r in read_manifest():
        if r["id"] == gid:
            return r["smiles"]
    raise SystemExit(f"no manifest smiles for {gid}")


def process(gid, backend):
    OUT.mkdir(parents=True, exist_ok=True)
    couples = _couples(gid)
    if not couples:
        return dict(id=gid, status="error", error="no couples (missing opt.xyz/result?)")
    neu_geom = DFT / gid / "neu" / "opt.xyz"
    if not neu_geom.exists():
        return dict(id=gid, status="error", error="no neu/opt.xyz reference")
    smiles = _smiles_for(gid)
    cfile = OUT / f"{gid}.freeze"
    n_frozen = rotatable_dihedral_constraints(smiles, cfile)

    # collect all distinct states across couples; neu is the frozen reference (geom as-is)
    states = {}
    for (sO, qO, mO), (sR, qR, mR) in couples:
        states[sO] = (qO, mO); states[sR] = (qR, mR)

    # conformer-matched geometry per state: neu -> its own opt; others -> re-opt from neu (frozen)
    cf_geom = {}
    rmsd_raw = {}; rmsd_cf = {}
    for st, (q, m) in states.items():
        raw = DFT / gid / st / "opt.xyz"
        if st == "neu":
            cf_geom[st] = neu_geom
            continue
        rmsd_raw[st] = round(heavy_rmsd(neu_geom, raw), 3) if raw.exists() else None
        outg = OUT / f"{gid}_{st}_cf.xyz"
        D.dft_smd(neu_geom, q, m, do_opt=True, do_gas=False, do_smd=False, do_freq=False,
                  opt_out=outg, constraints=str(cfile), conv_params=_CF_CONV,
                  assert_convergence=False, max_opt_steps=60, backend=backend)
        cf_geom[st] = outg
        rmsd_cf[st] = round(heavy_rmsd(neu_geom, outg), 3)

    rec = dict(id=gid, status="ok", smiles=smiles, n_frozen_dih=n_frozen,
               rmsd_neu_state_raw=rmsd_raw, rmsd_neu_state_cf=rmsd_cf, couples=[])
    for (sO, qO, mO), (sR, qR, mR) in couples:
        gO, gR = cf_geom[sO], cf_geom[sR]
        c = dict(couple=f"{sO}->{sR}", q_ox=qO, q_red=qR)
        try:
            eOO, _ = _e("gas", gO, qO, mO, backend); eOR, _ = _e("gas", gR, qO, mO, backend)
            eRO, dRO = _e("gas", gO, qR, mR, backend); eRR, _ = _e("gas", gR, qR, mR, backend)
            c["lambda_cf_gas_eV"] = round((eOR - eOO) + (eRO - eRR), 4)
            c["ion_unbound_gas"] = dRO[1]
        except Exception as e:
            c["lambda_cf_gas_eV"] = None; c["gas_err"] = f"{type(e).__name__}: {e}"
        try:
            sOO, _ = _e("smd", gO, qO, mO, backend); sOR, _ = _e("smd", gR, qO, mO, backend)
            sRO, _ = _e("smd", gO, qR, mR, backend); sRR, _ = _e("smd", gR, qR, mR, backend)
            c["lambda_cf_smd_eV"] = round((sOR - sOO) + (sRO - sRR), 4)
        except Exception as e:
            c["lambda_cf_smd_eV"] = None; c["smd_err"] = f"{type(e).__name__}: {e}"
        c["lambda_best_eV"], c["by"] = (
            (c["lambda_cf_smd_eV"], "cf+smd") if c.get("lambda_cf_smd_eV") is not None
            else (c.get("lambda_cf_gas_eV"), "cf_gas"))
        rec["couples"].append(c)
    (OUT / f"{gid}.json").write_text(json.dumps(rec, indent=1))
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", required=True, help="group id")
    ap.add_argument("--backend", default="gpu")
    a = ap.parse_args()
    try:
        r = process(a.only, a.backend)
    except Exception as e:
        r = dict(id=a.only, status="error", error=f"{type(e).__name__}: {e}",
                 trace=traceback.format_exc()[-800:])
        OUT.mkdir(parents=True, exist_ok=True)
        (OUT / f"{a.only}.json").write_text(json.dumps(r, indent=1))
    print(f"[{a.only}] {r.get('status')}")
    if r.get("status") == "ok":
        print(f"  frozen_dih={r['n_frozen_dih']}  rmsd(neu->ion) raw={r['rmsd_neu_state_raw']} cf={r['rmsd_neu_state_cf']}")
        for c in r["couples"]:
            print(f"  {c['couple']:<12} cf_gas={c.get('lambda_cf_gas_eV')} cf+smd={c.get('lambda_cf_smd_eV')} "
                  f"-> best={c.get('lambda_best_eV')} ({c.get('by')})")
    else:
        print("  ERROR:", r.get("error"))


if __name__ == "__main__":
    main()
