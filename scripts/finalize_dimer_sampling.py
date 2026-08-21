#!/usr/bin/env python
"""Wait for the top-K dimer-config DFT jobs, add xtb thermal, and report the viologen
dimerization dG = G(dimer) - 2 G(MV+.) as mean +/- spread over configs.
"""
from __future__ import annotations
import json, subprocess, time, statistics, os
from pathlib import Path

ROOT = Path("/nfs/lambda_stor_01/homes/rzhu/0_redox")
DFT = ROOT / "calcs" / "dft"
PY = "/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python"
CFGS = ["cfg00", "cfg01", "cfg03"]


def ready(c):
    p = DFT / "mv_dimer_cfg" / c / "result.json"
    try:
        return p.exists() and json.loads(p.read_text()).get("e_smd_eV") is not None
    except Exception:
        return False


def G(path):
    r = json.loads(Path(path).read_text())
    return r["e_smd_eV"] + (r.get("g_thermal_eV") or 0.0)


def main():
    t0 = time.time()
    while not all(ready(c) for c in CFGS):
        if time.time() - t0 > 7200:
            print("ERROR timeout", flush=True); return
        print(f"WAIT {[c for c in CFGS if not ready(c)]}", flush=True); time.sleep(60)
    print("DONE all dimer configs computed", flush=True)
    env = {"PYTHONPATH": "src", "PATH": os.environ["PATH"]}
    for c in CFGS:
        subprocess.run([PY, "-c",
            f"import sys;sys.path.insert(0,'src');import json;from redox.dft import _thermal_correction;"
            f"from ase.io import read;"
            f"p='calcs/dft/mv_dimer_cfg/{c}/result.json';r=json.load(open(p));"
            f"a=read('calcs/dft/mv_dimer_cfg/{c}/opt.xyz');"
            f"th=_thermal_correction(a,r['charge'],r['mult']);r.update(**th);"
            f"open(p,'w').write(json.dumps(r,indent=2));print('{c} thermal',th['g_thermal_eV'])"],
            cwd=ROOT, env=env)
    Gmono = G(DFT / "methyl_viologen" / "ox1" / "result.json")
    dGs = []
    for c in CFGS:
        dg = G(DFT / "mv_dimer_cfg" / c / "result.json") - 2 * Gmono
        dGs.append(dg)
        print(f"  {c}: dG_dim = {dg:+.3f} eV ({dg*96.485:+.1f} kJ/mol)", flush=True)
    mean = statistics.mean(dGs)
    sd = statistics.pstdev(dGs)
    print(f"\ndG_dim = {mean:+.3f} +/- {sd:.3f} eV  ({mean*96.485:+.1f} +/- {sd*96.485:.1f} kJ/mol)")
    print(f"config sigma = {sd*1000:.1f} meV (DFT). Dominant uncertainty is the +2 pimer continuum")
    print("desolvation, NOT config -> use dimerization as a RELATIVE ranking vs a reference.")
    print("RESULT dimer sampling complete", flush=True)


if __name__ == "__main__":
    main()
