#!/usr/bin/env python
"""Wait for the 3 microsolvated viologen DFT jobs, add xtb thermal, and compare the
explicit-solvation waves to bare continuum and experiment.

  nohup python scripts/finalize_microsolv.py > logs/microsolv/finalize.log 2>&1 &
"""
from __future__ import annotations
import json
import subprocess
import time
from pathlib import Path

ROOT = Path("/nfs/lambda_stor_01/homes/rzhu/0_redox")
DFT = ROOT / "calcs" / "dft"
PY = "/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python"
FC = 4.434749378960987
EXP = {"wave1 (MV2+/+.)": -0.45, "wave2 (MV+./0)": -0.88}
BARE = {"wave1 (MV2+/+.)": -0.757, "wave2 (MV+./0)": -1.318}   # bare ions + thermal
STATES = ["ox2", "ox1", "neu"]


def ready(st):
    p = DFT / "mv_solv" / st / "result.json"
    if not p.exists():
        return False
    try:
        return json.loads(p.read_text()).get("e_smd_eV") is not None
    except Exception:
        return False


def G(st):
    r = json.loads((DFT / "mv_solv" / st / "result.json").read_text())
    return r["e_smd_eV"] + (r.get("g_thermal_eV") or 0.0)


def main():
    t0 = time.time()
    while not all(ready(s) for s in STATES):
        if time.time() - t0 > 10800:
            print("ERROR timeout waiting for microsolv DFT", flush=True); return
        miss = [s for s in STATES if not ready(s)]
        print(f"WAIT missing {miss}", flush=True)
        time.sleep(60)
    print("DONE all 3 microsolvated states present", flush=True)

    # xtb thermal on the clusters (force)
    for s in STATES:
        subprocess.run([PY, "scripts/add_thermal.py", "--only", f"mv_solv:{s}", "--force"],
                       cwd=ROOT, env={"PYTHONPATH": "src",
                                      "PATH": __import__("os").environ["PATH"]})

    w1 = -(G("ox1") - G("ox2")) - FC       # MV2+ + e- -> MV+.
    w2 = -(G("neu") - G("ox1")) - FC       # MV+. + e- -> MV0
    calc = {"wave1 (MV2+/+.)": w1, "wave2 (MV+./0)": w2}

    print("\n" + "=" * 74)
    print("EXPLICIT MICROSOLVATION (MV . 4 MeCN, cluster-continuum) vs bare vs experiment")
    print("=" * 74)
    hdr = f"{'wave':18s} {'exp':>8s} {'bare+th':>9s} {'micro':>9s} {'|bare-exp|':>11s} {'|micro-exp|':>12s}"
    print(hdr); print("-" * len(hdr))
    maes = {"bare": [], "micro": []}
    for w in ["wave1 (MV2+/+.)", "wave2 (MV+./0)"]:
        eb = abs(BARE[w] - EXP[w]); em = abs(calc[w] - EXP[w])
        maes["bare"].append(eb); maes["micro"].append(em)
        print(f"{w:18s} {EXP[w]:+8.3f} {BARE[w]:+9.3f} {calc[w]:+9.3f} {eb:11.3f} {em:12.3f}")
    sp_exp = EXP['wave2 (MV+./0)'] - EXP['wave1 (MV2+/+.)']
    sp_bare = BARE['wave2 (MV+./0)'] - BARE['wave1 (MV2+/+.)']
    sp_micro = w2 - w1
    print(f"\n{'spacing':18s} {sp_exp:+8.3f} {sp_bare:+9.3f} {sp_micro:+9.3f}")
    print(f"\nMAE  bare+thermal = {sum(maes['bare'])/2:.3f} V   "
          f"microsolvated = {sum(maes['micro'])/2:.3f} V")
    print("RESULT microsolv viologen complete", flush=True)


if __name__ == "__main__":
    main()
