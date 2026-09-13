#!/usr/bin/env python
"""Cheap screen: is the VERTICAL radical anion (anion charge on the NEUTRAL geometry) BOUND?

One gas SCF per molecule -> gas HOMO. HOMO > 0 == unbound vertical anion == the 4-point inner-
sphere lambda is ill-defined (the E_R_at_O cross point is a 'neutral + free electron' energy, not
a bound reduced state). Distinguishes 'unbound -> exclude/flag' from 'bound -> potentially rescuable'.

  gpu_reserve run <idx> -- env PYTHONPATH=src python scripts/reorg_screen_unbound.py \
      --offset 0 --stride 2 --backend gpu
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import pandas as pd
import redox.dft as D

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results" / "d3tales_reorg_validation_prod"
OUT = BASE / "screen"; OUT.mkdir(parents=True, exist_ok=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--threshold", type=float, default=1.5)
    ap.add_argument("--offset", type=int, default=0)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--backend", default="gpu")
    a = ap.parse_args()
    df = pd.read_csv(BASE / "comparison.csv")
    df = df[df.our_smd_eV > a.threshold].reset_index(drop=True).iloc[a.offset::a.stride]
    for _, row in df.iterrows():
        gid = row["id"]; f = OUT / f"{gid}.json"
        if f.exists():
            print(f"[{gid}] skip"); print("DONE_MARKER"); continue
        neu = BASE / "geom" / f"{gid}_neu_opt.xyz"
        rec = dict(id=gid, our_smd_eV=round(float(row["our_smd_eV"]), 3),
                   d3tales_eV=row.get("d3tales_eV"), status="ok")
        try:
            r = D.dft_smd(neu, -1, 2, do_opt=False, do_gas=True, do_smd=False, do_freq=False,
                          backend=a.backend)
            rec["anion_homo_at_neu_eV"] = r.get("gas_homo_eV")
            rec["anion_unbound"] = bool(r.get("anion_unbound"))
            rec["anion_s_squared"] = r.get("gas_s_squared")
        except Exception as e:
            rec["status"] = "error"; rec["error"] = f"{type(e).__name__}: {e}"
        f.write_text(json.dumps(rec, indent=1))
        print(f"[{gid}] homo={rec.get('anion_homo_at_neu_eV')} unbound={rec.get('anion_unbound')} "
              f"our={rec['our_smd_eV']} d3={rec['d3tales_eV']}", flush=True)
        if rec["status"] == "ok":
            print("DONE_MARKER", flush=True)


if __name__ == "__main__":
    main()
