#!/usr/bin/env python
"""Independently recompute the scored pool's potentials at the xTB ("legacy") thermal tier.

The scored table data/manuscript/full_pool_scored_legacyG_MeCN_TBAPF6.csv (the basis of
Figures 5-7) has no producing script in the repository. This script rebuilds every couple
from the stored records and compares it with the table:

    G(state)  = E_SMD(active production single point; redox.core.protocol.load_record)
              + G_thermal(GFN2-xTB RRHO on the stored DFT geometry; dft._thermal_correction)
    E(O/R)    = G(O) - G(R) - E_abs(Fc),  E_abs(Fc) = G(Fc+) - G(Fc)  (same two rules)

Couples follow each candidate's declared ladder in the manifest: from the resting state
(n_e = 0) through n_e = -1, -2 for anolytes or +1, +2 for catholytes, as many as the table's n. A recomputed thermal value must match any stored
xTB value to 1e-4 eV. Runs no DFT. Writes a JSON record and prints the largest differences.

  OMP_NUM_THREADS=2 PYTHONPATH=src python scripts/validation/score_pool_legacy_tier.py --workers 8
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
os.environ.setdefault("XTB_BIN", shutil.which("xtb") or str(Path(sys.executable).with_name("xtb")))
TABLE = ROOT / "data/manuscript/full_pool_scored_legacyG_MeCN_TBAPF6.csv"
OUT = ROOT / "manuscript/evidence/2026-11-18_cpu_checks/pool_legacy_tier_rescore.json"


def ladders():
    from redox.core.common import read_manifest
    by = {}
    for r in read_manifest():
        by.setdefault(r["id"], {})[int(r["n_e"])] = r["state"]
    by.setdefault("ferrocene", {0: "neu", 1: "ox"})
    return by


def ladder(lad, row):
    step = -1 if row["pool"] == "anolyte" else 1
    return [lad[row["id"]][step * k] for k in range(int(row["n"]) + 1)]


def state_G(task):
    gid, state = task
    from ase.io import read as aseread
    from redox.backends.dft import _thermal_correction
    from redox.core.common import DFT
    from redox.core.protocol import load_record
    sd = DFT / gid / state
    try:
        res = json.loads((sd / "result.json").read_text())
        q, m = int(res["charge"]), int(res["mult"])
        rec = load_record(sd, q, m)
        if rec is None:
            return task, dict(error="no valid active single point")
        g = _thermal_correction(aseread(str(sd / "opt.xyz")), q, m)["g_thermal_eV"]
        stored = res.get("g_thermal_eV")
        return task, dict(e_smd_eV=rec["e_smd_eV"], g_thermal_eV=g, g_stored_eV=stored,
                          G_eV=rec["e_smd_eV"] + g, charge=q)
    except Exception as exc:  # noqa: BLE001 - reported, never hidden
        return task, dict(error=f"{type(exc).__name__}: {exc}"[:200])


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    rows = list(csv.DictReader(open(TABLE)))
    lad = ladders()
    tasks = {("ferrocene", "neu"), ("ferrocene", "ox")}
    for r in rows:
        tasks |= {(r["id"], s) for s in ladder(lad, r)}
    with ProcessPoolExecutor(max_workers=a.workers) as ex:
        G = dict(ex.map(state_G, sorted(tasks)))
    problems = [f"{k[0]}/{k[1]}: {v['error']}" for k, v in G.items() if "error" in v]
    problems += [f"{k[0]}/{k[1]}: recomputed thermal {v['g_thermal_eV']:.6f} != stored {v['g_stored_eV']:.6f}"
                 for k, v in G.items() if "error" not in v and isinstance(v.get("g_stored_eV"), (int, float))
                 and abs(v["g_thermal_eV"] - v["g_stored_eV"]) > 1e-4]
    fc = G[("ferrocene", "ox")]["G_eV"] - G[("ferrocene", "neu")]["G_eV"]
    out = []
    for r in rows:
        lst = ladder(lad, r)
        try:
            pairs = [(lst[i + 1], lst[i]) if r["pool"] == "catholyte" else (lst[i], lst[i + 1])
                     for i in range(len(lst) - 1)]          # (oxidized, reduced)
            Es = sorted(G[(r["id"], o)]["G_eV"] - G[(r["id"], rd)]["G_eV"] - fc for o, rd in pairs)
        except KeyError:
            out.append(dict(id=r["id"], error="missing state")); continue
        table = sorted(float(x) for x in r["E_couples"].split(";"))
        out.append(dict(id=r["id"], pool=r["pool"], recomputed_V=Es, table_V=table,
                        max_abs_diff_V=max(abs(x - y) for x, y in zip(Es, table)),
                        mean_recomputed_V=sum(Es) / len(Es), mean_table_V=float(r["E_mean"])))
    ok = [o for o in out if "max_abs_diff_V" in o]
    worst = sorted(ok, key=lambda o: -o["max_abs_diff_V"])[:8]
    summary = dict(n_candidates=len(rows), n_compared=len(ok), fc_legacy_V=fc,
                   max_abs_diff_V=max(o["max_abs_diff_V"] for o in ok),
                   n_within_1mV=sum(o["max_abs_diff_V"] <= 0.001 for o in ok),
                   n_within_5mV=sum(o["max_abs_diff_V"] <= 0.005 for o in ok), problems=problems)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    from redox.core.io import atomic_write_json
    atomic_write_json(OUT, dict(purpose="independent re-scoring of the pool table at the xTB thermal tier",
                                summary=summary, candidates=out,
                                states={f"{k[0]}/{k[1]}": v for k, v in G.items()}))
    print(json.dumps(summary, indent=1))
    for o in worst:
        print(f"{o['id']:34s} diff {o['max_abs_diff_V']:.4f} V  recomputed {[round(x,3) for x in o['recomputed_V']]}  table {o['table_V']}")
    return 0 if not problems else 1


if __name__ == "__main__":
    raise SystemExit(main())
