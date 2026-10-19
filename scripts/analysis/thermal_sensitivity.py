#!/usr/bin/env python
"""Historical gas-GFN2 RRHO sensitivity diagnostic.

Current native production thermochemistry uses projected internal-mode DFT qRRHO;
its 25/50/100 cm-1 terms are stored in addressed projected_thermal records.
This script instead computes the legacy GFN2-xTB quasi-RRHO model (sthr = 50
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
import hashlib
import math
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


def thermal_cache_path(sd, q, mult, sthr):
    sd = Path(sd)
    key = dict(geom_sha256=hashlib.sha256((sd / "opt.xyz").read_bytes()).hexdigest(),
               charge=q, mult=mult, sthr_cm=sthr, method="gfn2-xtb qRRHO", temperature_K=298.15)
    digest = hashlib.sha256(json.dumps(key, sort_keys=True).encode()).hexdigest()[:16]
    return sd / f"thermal_sthr{sthr}__{digest}.json", key


def referenced_potential(GO, GR, Gfc_ox, Gfc_neu):
    if any(x is None or not math.isfinite(x) for x in (GO, GR, Gfc_ox, Gfc_neu)):
        return None
    return GO - GR - (Gfc_ox - Gfc_neu)


def _xtb_g(xyz: Path, q: int, mult: int, sthr: int, threads: int):
    with tempfile.TemporaryDirectory(prefix="xtbsens_") as d:
        (Path(d) / "mol.xyz").write_text(xyz.read_text())
        (Path(d) / "x.inp").write_text(f"$thermo\n   sthr={sthr}\n$end\n")
        cmd = ["xtb", "mol.xyz", "--gfn", "2", "--chrg", str(q), "--uhf", str(mult - 1),
               "--hess", "--acc", "1.0", "--input", "x.inp"]
        r = subprocess.run(cmd, cwd=d, capture_output=True, text=True,
                           env=dict(os.environ, OMP_NUM_THREADS=str(threads), OMP_THREAD_LIMIT=str(threads),
                                    MKL_NUM_THREADS=str(threads), OPENBLAS_NUM_THREADS=str(threads)))
    m = re.search(r"G\(RRHO\) contrib\.\s+(-?\d+\.\d+)\s+Eh", r.stdout)
    n = re.search(r"#\s*imaginary freq\.\s+(\d+)", r.stdout)
    if r.returncode != 0 or not m:
        raise RuntimeError(f"xtb failed for {xyz}")
    return float(m.group(1)) * HARTREE_EV, int(n.group(1)) if n else None


def _job(args):
    sd, q, mult, sthr, threads = args
    sd = Path(sd)
    out, key = thermal_cache_path(sd, q, mult, sthr)
    if out.exists():
        record = json.loads(out.read_text())
        if record.get("input") == key and math.isfinite(record.get("g_thermal_eV", float("nan"))):
            return record
        raise RuntimeError(f"invalid addressed thermal cache: {out}; preserve and inspect it")
    g, nimag = _xtb_g(sd / "opt.xyz", q, mult, sthr, threads)
    rec = dict(sthr_cm=sthr, g_thermal_eV=g, n_imag=nimag, method="gfn2-xtb qRRHO", input=key)
    temp = out.with_suffix(f".{os.getpid()}.tmp")
    temp.write_text(json.dumps(rec, indent=2)); temp.replace(out)
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--only", help="comma-separated candidate IDs (Fc is always included)")
    ap.add_argument("--out", type=Path, default=ROOT / "results/thermal_sensitivity.csv")
    a = ap.parse_args()
    from redox.core.common import DFT, free_energy, read_result, state_names
    from redox.screening.scorecard import _candidate_ids
    ids = sorted(_candidate_ids())
    if a.only:
        chosen = set(a.only.split(","))
        if chosen - set(ids):
            ap.error(f"unknown candidates: {sorted(chosen - set(ids))}")
        ids = sorted(chosen)
    if a.workers < 1 or a.threads < 1 or a.workers * a.threads > max(1, len(os.sched_getaffinity(0)) // 2):
        ap.error("workers × threads exceeds shared-node CPU allowance")
    if os.getloadavg()[0] + a.workers * a.threads > (os.cpu_count() or 1):
        ap.error("shared node lacks current CPU capacity")
    if a.out.exists():
        ap.error("choose a new output file; previous sensitivity results are preserved")
    jobs = []
    for gid in sorted(set(ids) | {"ferrocene"}):
        for st in state_names(gid):
            r = read_result(gid, st, raw=True)
            for k in STHR:
                jobs.append((str(DFT / gid / st), int(r["charge"]), int(r["mult"]), k, a.threads))
    with ProcessPoolExecutor(a.workers) as ex:
        list(ex.map(_job, jobs))

    # E° and dG_disp per sthr, from the active-protocol E_smd + the alternative thermal term
    rows = []
    fc = {}
    for k in STHR:
        fc[k] = []
        for s in ("ox", "neu"):
            raw = read_result("ferrocene", s)
            if raw is None:
                raise RuntimeError("missing Fc reference state")
            path, _ = thermal_cache_path(DFT / "ferrocene" / s, int(raw["charge"]), int(raw["mult"]), k)
            th = json.loads(path.read_text())
            fc[k].append(raw["e_smd_eV"] + th["g_thermal_eV"] if raw.get("e_smd_eV") is not None else None)
    for gid in ids:
        st = sorted(((s, read_result(gid, s)) for s in state_names(gid)), key=lambda x: -x[1]["charge"])
        G = {k: {} for k in STHR}
        for s, r in st:
            for k in STHR:
                path, _ = thermal_cache_path(DFT / gid / s, int(r["charge"]), int(r["mult"]), k)
                th = json.loads(path.read_text())
                G[k][s] = (r["e_smd_eV"] + th["g_thermal_eV"]) if r.get("e_smd_eV") is not None else None
        for (sO, rO), (sR, rR) in zip(st, st[1:]):
            if rO["charge"] - rR["charge"] != 1:
                continue
            Es = [referenced_potential(G[k][sO], G[k][sR], *fc[k]) for k in STHR]
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
    out = a.out
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["id", "quantity"] + [f"sthr{k}" for k in STHR] + ["spread_V"])
        w.writeheader(); w.writerows(rows)
    sp = [r["spread_V"] for r in rows if r["spread_V"] != ""]
    print(f"wrote {out}: {len(rows)} quantities; max spread {max(sp):.4f} eV, "
          f"median {sorted(sp)[len(sp)//2]:.4f} eV" if sp else f"wrote {out} (no complete rows)")


if __name__ == "__main__":
    main()
