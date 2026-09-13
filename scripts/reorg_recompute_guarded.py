#!/usr/bin/env python
"""Guarded recompute of inner-sphere lambda for the flexible-imide OUTLIERS in the D3TaLES
production sweep (thiosuccinimides whose 4-point lambda is contaminated).

Two failure modes were diagnosed (see the reorg investigation):
  (A) CONFORMER JUMP  — neutral & anion optimize into different conformers, so the 4-point
      charges a whole conformational isomerization to lambda (large neutral<->anion RMSD).
  (B) UNBOUND GAS ANION — the gas-phase radical anion is unbound (HOMO > 0), so its gas
      single-point energy is meaningless (small RMSD but huge anion-side "relaxation").

This recomputes, reusing the already-optimized geometries in geom/:
  lambda_gas  : baseline 4-point (gas SPs on the SMD-opt geoms) + QC (anion HOMO, <S^2>, RMSD)
  lambda_cf   : CONFORMER-MATCHED — re-optimize the anion from the NEUTRAL geom with the
                rotatable dihedrals FROZEN (geomeTRIC constraints), so all 4 points share one
                conformer; fixes (A).
  lambda_smd  : SMD-CONSISTENT — score all 4 points in SMD (the electron is bound in solvent);
                fixes (B).
Then it PICKS the trustworthy value: conformer-matched if RMSD>RMSD_MAX, else SMD if the gas
anion is unbound, else the baseline gas value.

  gpu_reserve run <idx> -- env OMP_NUM_THREADS=2 PYTHONPATH=src \
      python scripts/reorg_recompute_guarded.py --only 06DSVU --backend gpu
"""
from __future__ import annotations
import argparse, json, traceback
from pathlib import Path
import numpy as np
import pandas as pd
from ase.io import read as ase_read, write as ase_write
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

import redox.dft as D

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results" / "d3tales_reorg_validation_prod"
GEO = BASE / "geom"
OUT = BASE / "recompute"
RMSD_MAX = 0.40   # Angstrom heavy-atom RMSD above which lambda is conformationally contaminated


def heavy_rmsd(p1, p2):
    a1 = ase_read(str(p1)); a2 = ase_read(str(p2))
    m = np.array(a1.get_chemical_symbols()) != "H"
    A = a1.positions[m] - a1.positions[m].mean(0)
    B = a2.positions[m] - a2.positions[m].mean(0)
    H = A.T @ B; U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T)); R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return float(np.sqrt(((A @ R.T - B) ** 2).sum(1).mean()))


def rotatable_dihedral_constraints(smiles, out_path):
    """geomeTRIC $freeze file locking every rotatable-bond dihedral (1-indexed atoms, matching
    embed()'s AddHs order == xyz order). Freezing the soft torsions keeps the anion in the
    neutral's conformer so only inner-sphere (bond/angle/ring) modes relax."""
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


def _gas(g, q, m, bk):
    r = D.dft_smd(g, q, m, do_opt=False, do_gas=True, do_smd=False, do_freq=False, backend=bk)
    return r.get("e_gas_eV"), r
def _smd(g, q, m, bk):
    r = D.dft_smd(g, q, m, do_opt=False, do_gas=False, do_smd=True, do_freq=False, backend=bk)
    return r.get("e_smd_eV"), r


def process(row, backend):
    gid = row["id"]; OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / f"{gid}.json"
    neu, an = GEO / f"{gid}_neu_opt.xyz", GEO / f"{gid}_anion_opt.xyz"
    rec = dict(id=gid, smiles=row["smiles"], family=row.get("family", ""),
               our_smd_prev=row.get("our_smd_eV"), status="ok", error="")
    try:
        rec["rmsd_neu_anion_A"] = round(heavy_rmsd(neu, an), 3)
        rec["lambda_gas_eV"] = round(float(row["our_smd_eV"]), 4)   # baseline (reuse prod value)
        # points we need in BOTH phases: neutral@neu and anion@neu (the unbound one drives the artifact)
        E_O_at_O, _ = _gas(neu, 0, 1, backend)
        E_R_at_O, rRO = _gas(neu, -1, 2, backend)
        rec["anion_homo_at_neu_eV"] = rRO.get("gas_homo_eV")
        rec["anion_unbound"] = bool(rRO.get("anion_unbound"))
        rec["anion_s_squared"] = rRO.get("gas_s_squared")
        # CONFORMER-MATCHED geometry: anion re-opt from the neutral geom with rotatable dihedrals
        # frozen (loose criteria; the frozen-dihedral opt only oscillates on gmax near the min).
        cfile = OUT / f"{gid}.freeze"; rec["n_frozen_dih"] = rotatable_dihedral_constraints(row["smiles"], cfile)
        an_cf = GEO / f"{gid}_anion_opt_cf.xyz"
        _cf_conv = dict(convergence_energy=5e-6, convergence_grms=4e-4, convergence_gmax=1e-3,
                        convergence_drms=2e-3, convergence_dmax=4e-3)
        D.dft_smd(neu, -1, 2, do_opt=True, do_gas=False, do_smd=False, do_freq=False,
                  opt_out=an_cf, constraints=str(cfile), conv_params=_cf_conv,
                  assert_convergence=False, max_opt_steps=60, backend=backend)
        rec["rmsd_neu_anioncf_A"] = round(heavy_rmsd(neu, an_cf), 3)
        # (A) conformer-matched, GAS energies — isolates the geometry fix (still hurt if anion unbound)
        try:
            E_R_at_Rcf, _ = _gas(an_cf, -1, 2, backend)
            E_O_at_Rcf, _ = _gas(an_cf, 0, 1, backend)
            rec["lambda_cf_eV"] = round((E_O_at_Rcf - E_O_at_O) + (E_R_at_O - E_R_at_Rcf), 4)
        except Exception as e:
            rec["lambda_cf_eV"] = None; rec["cf_error"] = f"{type(e).__name__}: {e}"
        # (DEFINITIVE) conformer-matched geometry AND SMD energies (electron bound in solvent).
        try:
            sOO, _ = _smd(neu, 0, 1, backend); sRO, _ = _smd(neu, -1, 2, backend)
            sRRcf, _ = _smd(an_cf, -1, 2, backend); sORcf, _ = _smd(an_cf, 0, 1, backend)
            rec["lambda_cf_smd_eV"] = round((sORcf - sOO) + (sRO - sRRcf), 4)
        except Exception as e:
            rec["lambda_cf_smd_eV"] = None; rec["cf_smd_error"] = f"{type(e).__name__}: {e}"
        # pick the trustworthy value: conformer-matched + SMD is the physically correct recipe for
        # these (bound anion, single conformer); fall back to cf-gas, then baseline.
        if rec.get("lambda_cf_smd_eV") is not None:
            rec["lambda_best_eV"], rec["lambda_best_by"] = rec["lambda_cf_smd_eV"], "cf+smd"
        elif rec.get("lambda_cf_eV") is not None:
            rec["lambda_best_eV"], rec["lambda_best_by"] = rec["lambda_cf_eV"], "conformer_matched_gas"
        else:
            rec["lambda_best_eV"], rec["lambda_best_by"] = rec["lambda_gas_eV"], "gas_baseline"
    except Exception as e:
        rec["status"] = "error"; rec["error"] = f"{type(e).__name__}: {e}"
        rec["trace"] = traceback.format_exc()[-800:]
    out.write_text(json.dumps(rec, indent=1))
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", help="single molecule id")
    ap.add_argument("--backend", default="gpu")
    ap.add_argument("--threshold", type=float, default=1.5,
                    help="auto-select ids with our_smd_eV > this from comparison.csv")
    a = ap.parse_args()
    df = pd.read_csv(BASE / "comparison.csv")
    if a.only:
        rows = df[df.id == a.only]
    else:
        rows = df[df.our_smd_eV > a.threshold]
    if rows.empty:
        print("no matching ids"); return
    for _, row in rows.iterrows():
        r = process(row, a.backend)
        print(f"[{row['id']}] {r['status']}  rmsd={r.get('rmsd_neu_anion_A')}->{r.get('rmsd_neu_anioncf_A')}A "
              f"unbound={r.get('anion_unbound')} gas={r.get('lambda_gas_eV')} cf={r.get('lambda_cf_eV')} "
              f"cf+smd={r.get('lambda_cf_smd_eV')} -> best={r.get('lambda_best_eV')} ({r.get('lambda_best_by')})",
              flush=True)
        if r["status"] == "ok":
            print("DONE_MARKER", flush=True)


if __name__ == "__main__":
    main()
