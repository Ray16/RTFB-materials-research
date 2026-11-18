#!/usr/bin/env python
"""Per-stage throughput report for pipeline campaigns (CPU only; reads records, runs nothing).

Sources, in order of preference:
  * run_manifest.json (and a resume's previous_run_manifest.json): wall time of every stage
    (build, uma, dft per redox state, sp, coverage, cross_points, review_tables), the host,
    and the GPU model from the node preflight;
  * calcs/dft/<id>/<state>/result.json  "timing_s" (written by redox.backends.dft since the
    throughput instrumentation): per-seed pre-optimization / optimization / SMD and gas
    single points, the thermal step, and atom / basis-function counts;
  * for runs that predate the instrumentation, per-seed optimizer time is parsed from the
    geomeTRIC "Time elapsed since start of run_optimizer" lines of the per-state DFT log;
    the remainder of the state wall time is reported as non-optimizer time (xTB pre-opt,
    SMD + gas single points with their SCF guesses, thermal).

Cached / resumed stages (wall < 60 s for a DFT state) are excluded, so only stages that
actually computed are reported.

    python scripts/validation/throughput_report.py RUN_DIR [RUN_DIR ...] [--out CSV]
"""
from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
CACHED_S = 60.0
_OPT = re.compile(r"Time elapsed since start of run_optimizer: ([\d.]+)")


def _manifests(run_dir: Path):
    found = [run_dir / "run_manifest.json"]
    found += sorted(run_dir.glob("resume_*/previous_run_manifest.json"))
    return [p for p in found if p.exists()]


def _gpu_names(m):
    return {g["index"]: g["name"] for g in (m.get("node_preflight") or {}).get("gpus", [])}


def _n_atoms(gid, state):
    p = ROOT / "calcs" / "dft" / gid / state / "opt.xyz"
    try:
        return int(p.read_text().split()[0])
    except (OSError, ValueError, IndexError):
        return None


def _record_timing(gid, state):
    p = ROOT / "calcs" / "dft" / gid / state / "result.json"
    try:
        d = json.loads(p.read_text())
    except (OSError, ValueError):
        return {}, {}
    return d.get("timing_s") or {}, d


def rows_for(run_dir: Path):
    out = []
    for mpath in _manifests(run_dir):
        m = json.loads(mpath.read_text())
        gpus = _gpu_names(m)
        cand = None
        for s in m.get("steps", []):
            argv = s.get("argv", [])
            if "--only" in argv:
                cand = argv[argv.index("--only") + 1]
            base = dict(run=run_dir.name, manifest=mpath.name, host=m.get("host"),
                        threads=m.get("threads"), nconf=m.get("nconf"), id=cand,
                        stage=s["stage"])
            if s["stage"] == "dft":
                for x in s.get("states", []):
                    a = x["attempts"][-1]
                    if not a.get("finished"):
                        continue
                    wall = a["finished"] - a["started"]
                    if wall < CACHED_S:
                        continue
                    timing, rec = _record_timing(cand, x["state"])
                    seeds = timing.get("seeds") if isinstance(timing.get("seeds"), list) else None
                    if seeds:
                        opt = sum(t.get("opt", 0) for t in seeds)
                        sp = sum(t.get("sp_smd", 0) + t.get("sp_gas", 0) for t in seeds)
                        src = "timing_s"
                    else:
                        try:
                            opt = sum(float(v) for v in _OPT.findall(Path(a["log"]).read_text()))
                        except OSError:
                            opt = None
                        sp = None
                        src = "log"
                    out.append(dict(base, state=x["state"], gpu=gpus.get(x.get("gpu"), x.get("gpu")),
                                    mult=rec.get("mult"), n_atoms=rec.get("n_atoms") or _n_atoms(cand, x["state"]),
                                    n_ao_sp=rec.get("n_ao_sp"), wall_min=round(wall / 60, 2),
                                    opt_min=(round(opt / 60, 2) if opt is not None else None),
                                    sp_min=(round(sp / 60, 2) if sp is not None else None),
                                    non_opt_min=(round((wall - opt) / 60, 2) if opt is not None else None),
                                    seed_spread_eV=rec.get("seed_spread_eV"), timing_source=src))
            else:
                st, fi = s.get("started"), s.get("finished")
                if st is None or fi is None:
                    continue
                wall = fi - st
                if s["stage"] in ("uma", "cross_points") and wall < 30:
                    continue   # cached on resume
                gidx = argv[argv.index("run") + 1] if "run" in argv else None
                out.append(dict(base, state="", gpu=gpus.get(int(gidx)) if gidx and gidx.isdigit() else "",
                                wall_min=round(wall / 60, 2)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("runs", nargs="+", type=Path)
    ap.add_argument("--out", type=Path, default=ROOT / "manuscript/tables/throughput_measured.csv")
    a = ap.parse_args()
    rows = [r for d in a.runs for r in rows_for(d.resolve())]
    cols = ["run", "manifest", "host", "gpu", "threads", "nconf", "id", "stage", "state", "mult",
            "n_atoms", "n_ao_sp", "wall_min", "opt_min", "sp_min", "non_opt_min",
            "seed_spread_eV", "timing_source"]
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {len(rows)} rows -> {a.out}")


if __name__ == "__main__":
    main()
