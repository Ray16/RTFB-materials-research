#!/usr/bin/env python
"""Resumable worker: recompute D3TaLES reorganization energies under THEIR protocol
(single RDKit conformer -> B3LYP/6-31G* gas-opt of neutral + ion -> Nelsen 4-point lambda_i)
for a stride-sharded slice of results/d3tales_reorg_validation/molecules.csv.

One JSON per molecule in .../calc/<id>.json (skipped if already present), so parallel workers
and reruns never redo finished molecules. Errors are recorded, never fatal.

  CUDA_VISIBLE_DEVICES=<idx> python scripts/validate_reorg_worker.py --offset 0 --stride 6
"""
from __future__ import annotations
import argparse, json, traceback
from pathlib import Path
import pandas as pd
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

from redox.nelsen import embed, gas_opt, four_point_lambda   # shared DFT 4-point machinery

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results" / "d3tales_reorg_validation"
CALC = BASE / "calc"
XC, BASIS = "b3lyp", "6-31g*"   # D3TaLES-matched protocol (RI-J on; ~5x faster, identical lambda)


def process(row):
    out = CALC / f"{row['id']}.json"
    if out.exists():
        try:
            if json.loads(out.read_text()).get("status") in ("ok", "partial"):
                return "skip"
        except Exception:
            pass
    rec = dict(id=row["id"], smiles=row["smiles"], family=row["family"],
               d3_electron=row["d3_electron"], d3_hole=row["d3_hole"],
               n_atoms=int(row["n_atoms"]), n_rot=int(row["n_rot"]),
               our_electron=None, our_hole=None, status="ok", error="")
    try:
        neu = embed(row["smiles"])
        neu_at, E_O_at_O = gas_opt(neu, 0, 0, XC, BASIS)
        if pd.notna(row["d3_electron"]):
            try:
                rec["our_electron"] = round(four_point_lambda(neu_at, E_O_at_O, -1, 1, XC, BASIS), 4)
            except Exception as e:
                rec["status"] = "partial"; rec["error"] += f"anion:{type(e).__name__} "
        if pd.notna(row["d3_hole"]):
            try:
                rec["our_hole"] = round(four_point_lambda(neu_at, E_O_at_O, +1, 1, XC, BASIS), 4)
            except Exception as e:
                rec["status"] = "partial"; rec["error"] += f"cation:{type(e).__name__} "
    except Exception as e:
        rec["status"] = "fail"; rec["error"] = f"{type(e).__name__}: {str(e)[:120]}"
        rec["trace"] = traceback.format_exc()[-400:]
    CALC.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=2))
    return rec["status"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--offset", type=int, default=0)
    ap.add_argument("--stride", type=int, default=1)
    a = ap.parse_args()
    df = pd.read_csv(BASE / "molecules.csv")
    mine = df.iloc[a.offset::a.stride]
    print(f"[worker {a.offset}/{a.stride}] {len(mine)} molecules", flush=True)
    for i, (_, row) in enumerate(mine.iterrows()):
        st = process(row)
        print(f"[w{a.offset}] {i+1}/{len(mine)} {row['id']} {row['family']} -> {st}", flush=True)


if __name__ == "__main__":
    main()
