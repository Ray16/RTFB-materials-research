#!/usr/bin/env python
"""Wait for the 6 viologen-fix DFT results (3 ion-pair species + 3 MV states recomputed
WITH thermal), then regenerate the main E table and the ion-pair comparison table.

Resumable / non-blocking: run in background, watch logs/viologen_fix/finalize.log.
"""
from __future__ import annotations
import json
import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path("/nfs/lambda_stor_01/homes/rzhu/0_redox")
DFT = ROOT / "calcs" / "dft"
PY = "/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python"

# (id, state, must_have_thermal)
NEED = [
    ("pf6", "anion", True),
    ("mv_ip1", "s0", True),
    ("mv_ip2", "s0", True),
    ("methyl_viologen", "neu", True),
    ("methyl_viologen", "ox1", True),
    ("methyl_viologen", "ox2", True),
]


def ready(idc, state, need_thermal):
    p = DFT / idc / state / "result.json"
    if not p.exists():
        return False
    try:
        r = json.loads(p.read_text())
    except Exception:
        return False
    if r.get("e_smd_eV") is None:
        return False
    if need_thermal and not isinstance(r.get("g_thermal_eV"), (int, float)):
        return False
    return True


def main():
    t0 = time.time()
    while True:
        missing = [f"{i}/{s}" for i, s, th in NEED if not ready(i, s, th)]
        if not missing:
            break
        if time.time() - t0 > 14400:
            print(f"ERROR timeout, still missing: {missing}", flush=True)
            sys.exit(1)
        print(f"WAIT {len(NEED)-len(missing)}/{len(NEED)} ready; missing {missing}", flush=True)
        time.sleep(60)
    print("DONE all 6 results present with thermal", flush=True)

    env = {"PYTHONPATH": "src", "PATH": os.environ["PATH"]}
    for mod in ("redox.redox", "redox.ionpair"):
        r = subprocess.run([PY, "-m", mod], cwd=ROOT, env=env, capture_output=True, text=True)
        out = r.stdout + "\n---STDERR---\n" + r.stderr
        (ROOT / "logs" / "viologen_fix" / f"{mod.split('.')[-1]}.out").write_text(out)
        print(f"DONE ran {mod} (exit {r.returncode})", flush=True)
        print(r.stdout, flush=True)
    print("RESULT viologen fix complete", flush=True)


if __name__ == "__main__":
    main()
