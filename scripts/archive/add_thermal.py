#!/usr/bin/env python
"""Thermal-only pass: add the gas-phase harmonic Gibbs correction (g_thermal_eV) to DFT
results that were computed WITHOUT it, reusing the already-optimized geometry.

Cheap by design: it does NOT redo the SMD opt or the wb97m-v single point. It reloads the
stored opt.xyz + the (charge, mult, opt_xc, opt_basis, opt_disp) from result.json and runs
ONLY the opt-level gas-phase Hessian (redox.dft._thermal_correction), then merges the result
back into result.json. Idempotent: skips states that already carry a numeric g_thermal_eV
unless --force.

  PYTHONPATH=src CUDA_VISIBLE_DEVICES=0 python scripts/add_thermal.py --only phenothiazine_parent --backend gpu
  PYTHONPATH=src python scripts/add_thermal.py --list        # just print what needs thermal
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DFT = ROOT / "calcs" / "dft"


def states_needing(force=False):
    todo = []
    for rj in sorted(DFT.glob("*/*/result.json")):
        gid, state = rj.parent.parent.name, rj.parent.name
        try:
            r = json.loads(rj.read_text())
        except Exception:
            continue
        if r.get("e_smd_eV") is None:
            continue  # energies not done; a full dft run must handle it, not us
        has = isinstance(r.get("g_thermal_eV"), (int, float))
        if has and not force:
            continue
        opt = rj.parent / "opt.xyz"
        if not opt.exists():
            print(f"[skip] {gid}/{state}: no opt.xyz (was it a --no-opt run?)")
            continue
        todo.append((gid, state))
    return todo


def add_one(gid, state, backend, force=False):
    from ase.io import read
    import redox.dft as D

    rj = DFT / gid / state / "result.json"
    r = json.loads(rj.read_text())
    if isinstance(r.get("g_thermal_eV"), (int, float)) and not force:
        print(f"[have] {gid}/{state} g_thermal={r['g_thermal_eV']:.4f}")
        return
    atoms = read(str(DFT / gid / state / "opt.xyz"))
    th = D._thermal_correction(atoms, r["charge"], r["mult"])   # GFN2-xTB RRHO (backend-agnostic)
    r.update(**th, freq_level="gfn2-xtb (RRHO, 298.15K)")
    rj.write_text(json.dumps(r, indent=2))
    print(f"[ok]   {gid}/{state} g_thermal={th['g_thermal_eV']:+.4f} eV "
          f"(n_imag={th['n_imag']}, freq_min={th['freq_min_cm']})", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None, help="'gid' or 'gid:state'")
    ap.add_argument("--shard", default=None, help="'n:i' round-robin over the todo list")
    ap.add_argument("--backend", default="gpu", choices=["cpu", "gpu"])
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--list", action="store_true", help="print todo and exit")
    args = ap.parse_args()

    todo = states_needing(force=args.force)
    if args.only:
        gid, _, st = args.only.partition(":")
        todo = [(g, s) for g, s in todo if g == gid and (not st or s == st)]
    if args.shard:
        n, i = (int(x) for x in args.shard.split(":"))
        todo = [t for k, t in enumerate(todo) if k % n == i]

    if args.list:
        for g, s in todo:
            print(f"{g}/{s}")
        print(f"# {len(todo)} states need thermal")
        return

    print(f"[start] {len(todo)} states: {todo}", flush=True)
    for gid, state in todo:
        try:
            add_one(gid, state, args.backend, force=args.force)
        except Exception as exc:
            print(f"[err]  {gid}/{state}: {type(exc).__name__}: {str(exc)[:160]}", flush=True)
    print("[done]", flush=True)


if __name__ == "__main__":
    main()
