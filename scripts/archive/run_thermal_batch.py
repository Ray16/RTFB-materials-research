#!/usr/bin/env python
"""Fan the thermal-only pass (scripts/add_thermal.py) across all idle LOCAL GPUs.

Picks free GPUs with scripts/free_gpus.py (respects the shared-node contention rule), then
runs one add_thermal worker per GPU, each taking a round-robin shard of the states that still
need g_thermal. Resumable: add_thermal skips states that already carry a numeric g_thermal_eV.

  nohup python scripts/run_thermal_batch.py > logs/thermal/driver.log 2>&1 &
"""
from __future__ import annotations
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
LOG = ROOT / "logs" / "thermal"
LOG.mkdir(parents=True, exist_ok=True)
PY = "/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python"


def free_local_gpus(nmax=8):
    out = subprocess.run([PY, "scripts/free_gpus.py", "-n", str(nmax)],
                         cwd=ROOT, capture_output=True, text=True)
    # local mode prints one comma-separated line, e.g. "1,5,6,7"
    txt = out.stdout.strip()
    if not txt:
        return []
    return [g.strip() for g in txt.replace("\n", ",").split(",") if g.strip()]


def n_todo():
    out = subprocess.run([PY, "scripts/add_thermal.py", "--list"],
                         cwd=ROOT, capture_output=True, text=True,
                         env={"PYTHONPATH": "src", "PATH": __import__("os").environ["PATH"]})
    return sum(1 for l in out.stdout.splitlines() if "/" in l and not l.startswith("#"))


def main():
    n = n_todo()
    print(f"[thermal] {n} states need g_thermal", flush=True)
    if n == 0:
        print("[thermal] nothing to do", flush=True)
        return
    gpus = free_local_gpus()
    if not gpus:
        print("[thermal] no idle GPU right now", flush=True)
        sys.exit(1)
    nsh = len(gpus)
    print(f"[thermal] fanning across {nsh} GPUs: {gpus}", flush=True)
    procs = []
    import os
    for i, gpu in enumerate(gpus):
        logf = open(LOG / f"worker_gpu{gpu}.log", "w")
        env = dict(os.environ, CUDA_VISIBLE_DEVICES=str(gpu), PYTHONPATH="src")
        p = subprocess.Popen(
            [PY, "scripts/add_thermal.py", "--shard", f"{nsh}:{i}", "--backend", "gpu"],
            cwd=ROOT, env=env, stdout=logf, stderr=subprocess.STDOUT)
        print(f"[thermal] worker shard {nsh}:{i} on GPU{gpu} (pid {p.pid})", flush=True)
        procs.append(p)
    for p in procs:
        p.wait()
    left = n_todo()
    print(f"[thermal] all workers done; {left} states still missing thermal", flush=True)
    print("RESULT thermal batch complete" if left == 0 else f"PARTIAL {left} remain", flush=True)


if __name__ == "__main__":
    main()
