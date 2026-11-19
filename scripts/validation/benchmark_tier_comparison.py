#!/usr/bin/env python
"""Re-score the frozen 39-endpoint benchmark at the POOL's thermal tier (CPU only).

The published benchmark (MAE 0.247 V) uses the lead-confirmation thermal tier: a native
solution-DFT Hessian with projected qRRHO ("native_solution_dft_qrrho_spin"), for every
molecule and for ferrocene. The 91-candidate pool is ranked with the cheaper tier: a GFN2-xTB
RRHO correction on the DFT geometry ("legacy"). This script asks how the benchmark changes
when the pool's tier is used instead, keeping every electronic energy fixed:

    G_legacy(state) = [G_native(state) - g_native(state)] + g_legacy(state)

G_native is the independently checked value in the frozen review (states.json). g_native is
the state's projected native thermal record (projected_thermal/*.json, the certificate-locked
projected qRRHO that defines this tier); the bracket is checked against the stored solvated
single-point energy before use (it must agree to 1e-6 eV). g_legacy is recomputed with the
pipeline's own GFN2-xTB RRHO routine (redox.backends.dft._thermal_correction) on the stored
optimized geometry; where an earlier g_legacy is stored it must agree to 1e-4 eV.
Ferrocene is treated the same way, so the reference moves with the tier.

No DFT is run. Writes one JSON record.

  OMP_NUM_THREADS=4 PYTHONPATH=src python scripts/validation/benchmark_tier_comparison.py
"""
from __future__ import annotations

import csv
import glob
import json
import os
import shutil
import statistics as st
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
os.environ.setdefault("XTB_BIN", shutil.which("xtb") or str(Path(sys.executable).with_name("xtb")))
from ase.io import read as aseread  # noqa: E402
from redox.backends.dft import _thermal_correction  # noqa: E402
from redox.core.common import DFT  # noqa: E402
from redox.core.io import atomic_write_json  # noqa: E402

# Frozen benchmark record (states.json/summary.json of the 2026-10-22 review are stored
# alongside it as review_states.json/benchmark_review.json). Raw QM records come from
# redox.core.common.DFT, i.e. $REDOX_DATA_ROOT/calcs/dft.
EVIDENCE = ROOT / "manuscript/evidence/2026-10-22_benchmark_complete"
OUT = ROOT / "manuscript/evidence/2026-11-18_cpu_checks/benchmark_legacy_tier.json"


def _native_g(gid, state):
    recs = sorted(glob.glob(str(DFT / gid / state / "projected_thermal" / "*.json")))
    vals = [json.loads(Path(r).read_text()).get("g_thermal_total_eV") for r in recs]
    return [v for v in vals if isinstance(v, (int, float))]


def _smd_energies(gid, state):
    es = []
    r = DFT / gid / state / "result.json"
    if r.exists():
        es.append(json.loads(r.read_text()).get("e_smd_eV"))
    for f in glob.glob(str(DFT / gid / state / "sp" / "*.json")):
        es.append(json.loads(Path(f).read_text()).get("e_smd_eV"))
    return [e for e in es if isinstance(e, (int, float))]


def main() -> int:
    states = json.loads((EVIDENCE / "review_states.json").read_text())
    fc_native = json.loads((EVIDENCE / "benchmark_review.json").read_text())["native_Fc_absolute_V"]
    points = list(csv.DictReader(open(EVIDENCE / "benchmark_points.csv")))
    needed = {f"{p['id']}:{s}" for p in points for s in p["event"].split("->")}
    needed |= {"ferrocene:neu", "ferrocene:ox"}

    per_state, problems = {}, []
    for key in sorted(needed):
        gid, state = key.split(":")
        G_nat = states[key]["G_eV"]
        cands = _native_g(gid, state)
        smd = _smd_energies(gid, state)
        # pick the native record whose implied electronic energy matches a stored SMD energy
        match = [(g, e) for g in cands for e in smd if abs((G_nat - g) - e) < 1e-6]
        if not match:
            problems.append(f"{key}: no native thermal record reproduces a stored E_smd")
            continue
        g_nat, e_smd = match[0]
        at = aseread(str(DFT / gid / state / "opt.xyz"))
        res = json.loads((DFT / gid / state / "result.json").read_text())
        g_leg = _thermal_correction(at, int(res["charge"]), int(res["mult"]))["g_thermal_eV"]
        stored = res.get("g_thermal_eV")
        if isinstance(stored, (int, float)) and abs(stored - g_leg) > 1e-4:
            problems.append(f"{key}: recomputed legacy thermal {g_leg:.6f} != stored {stored:.6f}")
        per_state[key] = dict(G_native_eV=G_nat, g_native_eV=g_nat, e_smd_eV=e_smd,
                              g_legacy_eV=g_leg, g_legacy_stored_eV=stored,
                              G_legacy_eV=e_smd + g_leg)
        print(f"{key:45s} g_nat {g_nat:8.4f}  g_leg {g_leg:8.4f}  diff {g_leg-g_nat:+.4f}", flush=True)

    if problems:
        print("\n".join(problems))
    fc_leg = per_state["ferrocene:ox"]["G_legacy_eV"] - per_state["ferrocene:neu"]["G_legacy_eV"]
    fc_nat_check = per_state["ferrocene:ox"]["G_native_eV"] - per_state["ferrocene:neu"]["G_native_eV"]

    rows = []
    for p in points:
        ox, red = p["event"].split("->")
        a, b = per_state.get(f"{p['id']}:{ox}"), per_state.get(f"{p['id']}:{red}")
        if not a or not b:
            continue
        e_nat = a["G_native_eV"] - b["G_native_eV"] - fc_native
        e_leg = a["G_legacy_eV"] - b["G_legacy_eV"] - fc_leg
        rows.append(dict(id=p["id"], event=p["event"], family=p["family"], wave=int(p["wave"]),
                         exp_V=float(p["exp_V"]), calc_native_V=e_nat,
                         calc_native_published_V=float(p["calc_V"]), calc_legacy_V=e_leg,
                         shift_V=e_leg - e_nat))

    def metrics(key, subset):
        err = [r[key] - r["exp_V"] for r in subset]
        return dict(n=len(err), MAE_V=st.mean(abs(e) for e in err),
                    RMSE_V=(st.mean(e * e for e in err)) ** 0.5, bias_V=st.mean(err),
                    max_abs_V=max(abs(e) for e in err))

    ox_fams = {"phenothiazine", "nitroxide"}
    summary = dict(
        n_endpoints=len(rows),
        native_recomputed=metrics("calc_native_V", rows),
        legacy=metrics("calc_legacy_V", rows),
        legacy_reduction=metrics("calc_legacy_V", [r for r in rows if r["family"] not in ox_fams]),
        legacy_oxidation=metrics("calc_legacy_V", [r for r in rows if r["family"] in ox_fams]),
        shift_V=dict(mean=st.mean(r["shift_V"] for r in rows),
                     min=min(r["shift_V"] for r in rows), max=max(r["shift_V"] for r in rows)),
        max_native_reproduction_error_V=max(abs(r["calc_native_V"] - r["calc_native_published_V"])
                                            for r in rows),
        fc_native_V=fc_native, fc_native_from_states_V=fc_nat_check, fc_legacy_V=fc_leg,
    )
    record = dict(purpose="frozen benchmark re-scored at the pool's GFN2-xTB RRHO thermal tier",
                  summary=summary, problems=problems, endpoints=rows, states=per_state,
                  inputs=dict(states=str(EVIDENCE / "review_states.json"),
                              summary=str(EVIDENCE / "benchmark_review.json"),
                              points=str((EVIDENCE / "benchmark_points.csv").relative_to(ROOT))))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_json(OUT, record)
    print(json.dumps(summary, indent=1))
    return 0 if not problems and len(rows) == len(points) else 1


if __name__ == "__main__":
    raise SystemExit(main())
