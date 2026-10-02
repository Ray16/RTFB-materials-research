#!/usr/bin/env python
"""Sensitivity of E° and dG_disp to the thermal (RRHO) model, for the screening candidates.

Production thermal = GFN2-xTB quasi-RRHO (Grimme rotor interpolation of low modes, sthr = 50
cm-1) evaluated on the DFT-SMD geometry, which is NOT an xTB stationary point — some states
carry 1-2 imaginary modes. Rather than trust one setting, recompute the thermal term with the
rotor cutoff at 25 / 50 / 100 cm-1 and report the spread of every candidate E° and dG_disp.
A property whose spread exceeds its sigma is thermal-model-limited.

Writes calcs/dft/<id>/<state>/thermal_sthr<k>.json (never touches result.json) and
results/thermal_sensitivity.csv.

  python scripts/analysis/thermal_sensitivity.py --workers 4 --threads 2
"""
from __future__ import annotations
import argparse
import csv
import json
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
HARTREE_EV = 27.211386245988
STHR = (25, 50, 100)


def _xtb_g(xyz: Path, q: int, mult: int, sthr: int, threads: int):
    with tempfile.TemporaryDirectory(prefix="xtbsens_") as d:
        (Path(d) / "mol.xyz").write_text(xyz.read_text())
        (Path(d) / "x.inp").write_text(f"$thermo\n   sthr={sthr}\n$end\n")
        cmd = ["xtb", "mol.xyz", "--gfn", "2", "--chrg", str(q), "--uhf", str(mult - 1),
               "--hess", "--acc", "1.0", "--input", "x.inp"]
        r = subprocess.run(cmd, cwd=d, capture_output=True, text=True,
                           env=dict(os.environ, OMP_NUM_THREADS=str(threads)))
    m = re.search(r"G\(RRHO\) contrib\.\s+(-?\d+\.\d+)\s+Eh", r.stdout)
    n = re.search(r"#\s*imaginary freq\.\s+(\d+)", r.stdout)
    if not m:
        raise RuntimeError(f"xtb failed for {xyz}")
    return float(m.group(1)) * HARTREE_EV, int(n.group(1)) if n else None


def _job(args):
    sd, q, mult, sthr, threads = args
    sd = Path(sd)
    out = sd / f"thermal_sthr{sthr}.json"
    if out.exists():
        return json.loads(out.read_text())
    g, nimag = _xtb_g(sd / "opt.xyz", q, mult, sthr, threads)
    rec = dict(sthr_cm=sthr, g_thermal_eV=g, n_imag=nimag, method="gfn2-xtb qRRHO")
    out.write_text(json.dumps(rec, indent=2))
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--threads", type=int, default=2)
    a = ap.parse_args()
    from redox.core.common import DFT, free_energy, read_result, state_names
    from redox.screening.scorecard import _candidate_ids
    ids = sorted(_candidate_ids())
    jobs = []
    for gid in ids:
        for st in state_names(gid):
            r = read_result(gid, st, raw=True)
            for k in STHR:
                jobs.append((str(DFT / gid / st), int(r["charge"]), int(r["mult"]), k, a.threads))
    with ProcessPoolExecutor(a.workers) as ex:
        list(ex.map(_job, jobs))

    # E° and dG_disp per sthr, from the active-protocol E_smd + the alternative thermal term
    rows = []
    for gid in ids:
        st = sorted(((s, read_result(gid, s)) for s in state_names(gid)), key=lambda x: -x[1]["charge"])
        G = {k: {} for k in STHR}
        for s, r in st:
            for k in STHR:
                th = json.loads((DFT / gid / s / f"thermal_sthr{k}.json").read_text())
                G[k][s] = (r["e_smd_eV"] + th["g_thermal_eV"]) if r.get("e_smd_eV") is not None else None
        for (sO, rO), (sR, rR) in zip(st, st[1:]):
            if rO["charge"] - rR["charge"] != 1:
                continue
            Es = [(-(G[k][sR] - G[k][sO]) if None not in (G[k][sR], G[k][sO]) else None) for k in STHR]
            rows.append(dict(id=gid, quantity=f"E({sO}->{sR})",
                             **{f"sthr{k}": (round(e, 4) if e is not None else "") for k, e in zip(STHR, Es)},
                             spread_V=(round(max(Es) - min(Es), 4) if None not in Es else "")))
        for i in range(1, len(st) - 1):
            (sH, _), (sM, _), (sL, _) = st[i - 1], st[i], st[i + 1]
            ds = [((G[k][sH] + G[k][sL] - 2 * G[k][sM]) if None not in (G[k][sH], G[k][sL], G[k][sM])
                   else None) for k in STHR]
            rows.append(dict(id=gid, quantity=f"dG_disp({sM}) eV",
                             **{f"sthr{k}": (round(d, 4) if d is not None else "") for k, d in zip(STHR, ds)},
                             spread_V=(round(max(ds) - min(ds), 4) if None not in ds else "")))
    out = ROOT / "results" / "thermal_sensitivity.csv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "quantity"] + [f"sthr{k}" for k in STHR] + ["spread_V"])
        w.writeheader(); w.writerows(rows)
    sp = [r["spread_V"] for r in rows if r["spread_V"] != ""]
    print(f"wrote {out}: {len(rows)} quantities; max spread {max(sp):.4f} eV, "
          f"median {sorted(sp)[len(sp)//2]:.4f} eV" if sp else f"wrote {out} (no complete rows)")


if __name__ == "__main__":
    main()
