#!/usr/bin/env python
"""Measured legacy (xTB msRRHO) vs native (projected internal-mode qRRHO) thermal, per state,
and the net effect on redox couples. Runs xtb --hess once per state: parse G(RRHO) from stdout
(legacy) and the cartesian `hessian` file -> qm.qrrho.nuclear_thermals (native). CPU only."""
import os, sys, subprocess, tempfile, re
from pathlib import Path
import numpy as np
sys.path.insert(0, "/nfs/lambda_stor_01/homes/rzhu/0_redox_refactor/src")
from ase.io import read as aseread
from redox.qm import qrrho as Q

XTB = "/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/xtb"
ARCH = "/nfs/lambda_stor_01/homes/rzhu/archive/0_redox/calcs/dft"
HARTREE_EV = 27.211386245988
Z = {"H":1,"C":6,"N":7,"O":8,"F":9,"Fe":26,"S":16}

def run_state(gid, state, charge, mult):
    geom = f"{ARCH}/{gid}/{state}/opt.xyz"
    at = aseread(geom)
    with tempfile.TemporaryDirectory() as d:
        xyz = Path(d)/"mol.xyz"; at.write(str(xyz))
        env = dict(os.environ, OMP_NUM_THREADS="4", MKL_NUM_THREADS="4", XTB_THREADS="4")
        p = subprocess.run([XTB, "mol.xyz", "--gfn", "2", "--chrg", str(charge),
                            "--uhf", str(mult-1), "--hess", "--acc", "1.0"],
                           cwd=d, env=env, capture_output=True, text=True)
        m = re.search(r"G\(RRHO\) contrib\.\s+([-\d.]+)", p.stdout)
        if not m:  # alt label
            m = re.search(r"\|\s*G\(RRHO\)\s+([-\d.]+)", p.stdout)
        legacy_g = float(m.group(1)) * HARTREE_EV if m else None
        # parse cartesian hessian (Turbomole $hessian ... $end), Hartree/Bohr^2
        htxt = (Path(d)/"hessian").read_text()
        body = htxt.split("$hessian",1)[1].split("$end",1)[0]
        vals = [float(x) for x in re.findall(r"[-+]?\d+\.\d+(?:[eE][-+]?\d+)?", body)]
        n3 = 3*len(at); H = np.array(vals[:n3*n3]).reshape(n3, n3)
        numbers = [Z[s] for s in at.get_chemical_symbols()]
        nt = Q.nuclear_thermals(H, at.positions, numbers, 1)
        native_g = nt["nuclear_by_cutoff_eV"]["50"]   # 50 cm^-1 rotor cutoff (standard)
    return legacy_g, native_g

CASES = {
 "ferrocene": [("neu",0,1),("ox",1,2)],
 "exp_dmpz_mb_methyl": [("neu",0,1),("ox1",1,2),("ox2",2,1)],
}
for gid, states in CASES.items():
    print(f"\n=== {gid} ===", flush=True)
    G = {}
    for st, q, mult in states:
        lg, ng = run_state(gid, st, q, mult)
        G[st] = (lg, ng)
        print(f"  {st}: legacy G(RRHO)={lg:.4f} eV   native qRRHO(50)={ng:.4f} eV   diff={ng-lg:+.4f}", flush=True)
    order = [s for s,_,_ in states]
    for a, b in zip(order, order[1:]):
        dleg = G[b][0]-G[a][0]; dnat = G[b][1]-G[a][1]
        print(f"  couple {a}->{b}: thermal contrib legacy {dleg:+.4f} / native {dnat:+.4f} eV  "
              f"=> potential shift from basis change = {(dnat-dleg)*1000:+.1f} mV", flush=True)
