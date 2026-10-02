#!/usr/bin/env python
"""Backfill the GFN2-xTB RRHO thermal correction for OROP benchmark states that lack one.

Most OROP benchmark states were computed with do_freq off, so their result.json has
g_thermal_eV = None — while the PRODUCTION pipeline always uses G = E_smd + g_thermal. A
benchmark that silently drops the thermal term does not measure the production protocol
(and a state pair with thermal on only ONE side is off by several eV). This computes the
missing term on the DFT+SMD-optimized geometry (opt.xyz) and writes it to a SEPARATE
`thermal.json` next to result.json — result.json is never modified.

CPU only (xtb). Threads are capped per worker so the shared node is not oversubscribed.

  python scripts/validation/orop/orop_backfill_thermal.py --workers 8 --threads 4
"""
from __future__ import annotations
import argparse
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
CALC = ROOT / "calcs" / "orop"


def _todo():
    out = []
    for rj in sorted(CALC.glob("*/*/result.json")):
        sd = rj.parent
        if (sd / "thermal.json").exists() or not (sd / "opt.xyz").exists():
            continue
        if json.loads(rj.read_text()).get("g_thermal_eV") is not None:
            continue
        out.append(sd)
    return out


def _one(sd: str):
    from ase.io import read
    from redox.qm.dft import _thermal_correction
    sd = Path(sd)
    r = json.loads((sd / "result.json").read_text())
    th = _thermal_correction(read(str(sd / "opt.xyz")), int(r["charge"]), int(r["mult"]))
    th.update(geometry=str(sd / "opt.xyz"), charge=int(r["charge"]), mult=int(r["mult"]),
              freq_level="gfn2-xtb (RRHO, 298.15K)")
    (sd / "thermal.json").write_text(json.dumps(th, indent=2))
    return str(sd), th["g_thermal_eV"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--threads", type=int, default=4, help="xtb OMP threads per worker")
    args = ap.parse_args()
    os.environ["XTB_THREADS"] = str(args.threads)
    for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
        os.environ[k] = str(args.threads)
    todo = _todo()
    print(f"{len(todo)} states need a thermal correction", flush=True)
    with ProcessPoolExecutor(args.workers) as ex:
        futs = {ex.submit(_one, str(sd)): sd for sd in todo}
        for i, f in enumerate(as_completed(futs), 1):
            try:
                sd, g = f.result()
                print(f"[{i}/{len(todo)}] {sd}  g_thermal={g:.4f} eV", flush=True)
            except Exception as exc:
                print(f"[{i}/{len(todo)}] FAIL {futs[f]}: {type(exc).__name__}: {exc}", flush=True)


if __name__ == "__main__":
    main()
