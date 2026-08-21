#!/usr/bin/env python
"""Unattended driver for the viologen ion-pair fix.

Waits for the ion-pair UMA geometries, fans the DFT+SMD+Hessian jobs onto IDLE cluster GPUs
(never colliding — uses scripts/free_gpus.py), then runs redox.redox (main table, now with
thermal) and redox.ionpair (both viologen waves, released-counterion scheme). Fully resumable:
skips any species whose calcs/dft result already exists.

  nohup python scripts/run_viologen_fix.py > logs/viologen_fix/driver.log 2>&1 &

Emits single-line progress markers (LAUNCH/DONE/RESULT/ERROR) so a Monitor can watch it.
"""
from __future__ import annotations
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path("/nfs/lambda_stor_01/homes/rzhu/0_redox")
LOG = ROOT / "logs" / "viologen_fix"
LOG.mkdir(parents=True, exist_ok=True)
HOSTS = "lambda2,lambda4"
ACT = "source ~/miniforge3/etc/profile.d/conda.sh && conda activate redox"
PY = "/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python"

# id -> (states to expect, extra dft flags)
JOBS = {
    "pf6":             (["anion"], ""),
    "mv_ip2":          (["s0"], ""),
    "mv_ip1":          (["s0"], ""),
    "methyl_viologen": (["ox2", "ox1", "neu"], "--force"),   # recompute WITH thermal
}
# UMA geoms that DFT depends on (pf6 + methyl_viologen already relaxed earlier)
UMA_WAIT = {"mv_ip2": "s0", "mv_ip1": "s0"}


def result_ok(idc, state):
    p = ROOT / "calcs" / "dft" / idc / state / "result.json"
    if not p.exists():
        return False
    import json
    try:
        return json.loads(p.read_text()).get("e_smd_eV") is not None
    except Exception:
        return False


def species_done(idc):
    states, extra = JOBS[idc]
    # A --force species must always recompute (e.g. to add thermal to a pre-thermal result),
    # so never treat an existing result as "done" for it.
    if "--force" in extra:
        return False
    return all(result_ok(idc, s) for s in states)


def free_slots(n):
    out = subprocess.run([PY, "scripts/free_gpus.py", "--hosts", HOSTS, "-n", str(n)],
                         cwd=ROOT, capture_output=True, text=True)
    return [tuple(l.split()) for l in out.stdout.strip().splitlines() if l.strip()]


def launch(idc, host, gpu):
    _, extra = JOBS[idc]
    logf = open(LOG / f"dft_{idc}.log", "w")
    cmd = (f"cd {ROOT} && {ACT} && CUDA_VISIBLE_DEVICES={gpu} PYTHONPATH=src "
           f"python -m redox.dft --only {idc} --backend gpu {extra}")
    p = subprocess.Popen(["ssh", "-o", "BatchMode=yes", "-o", "ServerAliveInterval=60",
                          host, cmd], stdout=logf, stderr=subprocess.STDOUT)
    print(f"LAUNCH {idc} on {host} GPU{gpu} (force={'--force' in extra})", flush=True)
    return p


def main():
    t0 = time.time()

    # 1) wait for the ion-pair UMA geometries (pf6 + MV already done earlier)
    while True:
        missing = [i for i, s in UMA_WAIT.items()
                   if not (ROOT / "calcs" / "uma" / i / s / "relaxed.xyz").exists()]
        if not missing:
            break
        if time.time() - t0 > 7200:
            print(f"ERROR UMA geoms never appeared: {missing}", flush=True); sys.exit(1)
        time.sleep(20)
    print("DONE uma geometries ready", flush=True)

    # 2) launch DFT for species not already done, one per idle GPU (re-checked now)
    todo = [i for i in JOBS if not species_done(i)]
    print(f"DFT todo: {todo}", flush=True)
    procs = {}
    if todo:
        # try to get a slot per job; if fewer GPUs, launch what we can and loop
        pending = list(todo)
        while pending:
            slots = free_slots(len(pending))
            if not slots:
                print("WAIT no idle GPU, retrying in 60s", flush=True)
                time.sleep(60)
                continue
            for idc, (host, gpu) in zip(list(pending), slots):
                procs[idc] = launch(idc, host, gpu)
                pending.remove(idc)
            if pending:
                time.sleep(30)   # let launches settle before re-scanning
        for idc, p in procs.items():
            p.wait()
            print(f"DONE dft {idc} (exit {p.returncode})", flush=True)

    # 3) verify all done
    bad = [i for i in JOBS if not species_done(i)]
    if bad:
        print(f"ERROR species missing results: {bad}", flush=True); sys.exit(1)
    print("DONE all dft results present", flush=True)

    # 4) recompute main table (now with thermal for MV) + the ion-pair waves
    for mod in ("redox.redox", "redox.ionpair"):
        r = subprocess.run([PY, "-m", mod], cwd=ROOT,
                           env={"PYTHONPATH": "src", "PATH": __import__("os").environ["PATH"]},
                           capture_output=True, text=True)
        (LOG / f"{mod.split('.')[-1]}.out").write_text(r.stdout + "\n---STDERR---\n" + r.stderr)
        print(f"DONE ran {mod} (exit {r.returncode})", flush=True)

    print("RESULT viologen fix complete — see logs/viologen_fix/ionpair.out", flush=True)


if __name__ == "__main__":
    main()
