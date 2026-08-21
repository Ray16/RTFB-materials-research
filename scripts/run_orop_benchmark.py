#!/usr/bin/env python
"""Blind external benchmark of our DFT+SMD pipeline against OROP experimental MeCN potentials.

OROP (data/raw/validation/OROP) ships, per oxidation system s (1..193):
  - two geometries  <s>-1-MeCN.xyz (OXIDIZED)  and  <s>-2-MeCN.xyz (REDUCED),
  - charge/spin of the oxidized state (feature CSV row s-1; verified aligned by e- parity),
  - an experimental redox potential (V vs Fc, implicit_solvation_results.csv).

We re-optimize BOTH states in SMD(MeCN) at OUR level (r2scan/def2-svp(d) // wb97m-v/def2-tzvp),
add a GFN2-xTB RRHO thermal correction, reference to our own level-matched Fc/Fc+, and compare
E_vs_Fc to experiment. Nothing is fit. Geometries only seed our optimization.

  python scripts/run_orop_benchmark.py --only 103 --backend gpu        # one system, both states
  python scripts/run_orop_benchmark.py --fanout 103,96,...             # distribute across idle GPUs
  python scripts/run_orop_benchmark.py --aggregate 103,96,...          # build results/orop_benchmark.csv
"""
from __future__ import annotations
import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OROP = ROOT / "data" / "raw" / "validation" / "OROP"
XYZ = OROP / "implicit_optimized_xyz"
CALC = ROOT / "calcs" / "orop"
PY = "/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python"


def load_index():
    """sys -> dict(charge_ox, mult_ox, charge_red, mult_red, exp, xyz_ox, xyz_red)."""
    feat = np.loadtxt(OROP / "features" / "feature_input_implicit_OROP_B3LYP-D3.csv",
                      skiprows=2, delimiter=",")
    res = pd.read_csv(OROP / "implicit_solvation_results.csv")
    res.columns = [c.strip() for c in res.columns]
    exp = res["exp. redox pot. [V]"].values
    imp = res["implicit solvation redox pot. [V]"].values
    idx = {}
    for i in range(len(feat)):
        s = i + 1
        q_ox = int(feat[i, 3]); m_ox = int(feat[i, 2])
        q_red = q_ox - 1
        # reduced state has one more electron -> parity flips; default low-spin
        m_red = 1 if (m_ox == 2) else 2   # ox doublet->red singlet ; ox singlet->red doublet
        fo, fr = XYZ / f"{s}-1-MeCN.xyz", XYZ / f"{s}-2-MeCN.xyz"
        if not (fo.exists() and fr.exists()):
            continue
        idx[s] = dict(charge_ox=q_ox, mult_ox=m_ox, charge_red=q_red, mult_red=m_red,
                      exp=float(exp[i]), imp_dft=float(imp[i]), xyz_ox=fo, xyz_red=fr)
    return idx


def _done(s, state):
    p = CALC / str(s) / state / "result.json"
    if not p.exists():
        return False
    try:
        return json.loads(p.read_text()).get("e_smd_eV") is not None
    except Exception:
        return False


def run_system(s, backend, force=False):
    import redox.dft as D
    idx = load_index()
    if s not in idx:
        print(f"[skip] sys {s}: no geometry/feature", flush=True); return
    rec = idx[s]
    for state, xyz, q, m in [("ox", rec["xyz_ox"], rec["charge_ox"], rec["mult_ox"]),
                             ("red", rec["xyz_red"], rec["charge_red"], rec["mult_red"])]:
        outdir = CALC / str(s) / state
        if _done(s, state) and not force:
            print(f"[have] sys {s}/{state}", flush=True); continue
        outdir.mkdir(parents=True, exist_ok=True)
        print(f"[run ] sys {s}/{state} q={q} m={m} ({xyz.name})", flush=True)
        try:
            res = D.dft_smd(xyz, q, m, do_opt=True, do_freq=True, backend=backend,
                            opt_out=outdir / "opt.xyz")
            (outdir / "result.json").write_text(json.dumps(res, indent=2))
            print(f"[ok  ] sys {s}/{state} e_smd={res['e_smd_eV']:.3f} "
                  f"gth={res.get('g_thermal_eV')}", flush=True)
        except Exception as exc:
            print(f"[err ] sys {s}/{state}: {type(exc).__name__}: {str(exc)[:160]}", flush=True)


def fanout(systems, force=False):
    """Distribute systems across idle local GPUs, one system (both states) per GPU worker."""
    gpus_txt = subprocess.run([PY, "scripts/free_gpus.py", "-n", "8"], cwd=ROOT,
                              capture_output=True, text=True).stdout.strip()
    gpus = [g for g in gpus_txt.replace("\n", ",").split(",") if g.strip()]
    if not gpus:
        print("no idle GPU", flush=True); sys.exit(1)
    todo = [s for s in systems if not (_done(s, "ox") and _done(s, "red")) or force]
    print(f"[orop] {len(todo)} systems to run across GPUs {gpus}", flush=True)
    LOG = ROOT / "logs" / "orop"; LOG.mkdir(parents=True, exist_ok=True)
    running = {}   # gpu -> proc
    queue = list(todo)
    import time
    while queue or running:
        for g in list(gpus):
            if g not in running and queue:
                s = queue.pop(0)
                env = dict(os.environ, CUDA_VISIBLE_DEVICES=str(g), PYTHONPATH="src")
                lf = open(LOG / f"sys{s}.log", "w")
                cmd = [PY, "scripts/run_orop_benchmark.py", "--only", str(s), "--backend", "gpu"]
                if force:
                    cmd.append("--force")
                running[g] = (subprocess.Popen(cmd, cwd=ROOT, env=env, stdout=lf,
                                               stderr=subprocess.STDOUT), s)
                print(f"[orop] launched sys {s} on GPU{g}", flush=True)
        for g in list(running):
            p, s = running[g]
            if p.poll() is not None:
                print(f"[orop] sys {s} on GPU{g} done (exit {p.returncode})", flush=True)
                del running[g]
        time.sleep(5)
    print("RESULT orop fanout complete", flush=True)


def aggregate(systems):
    import redox.dft  # noqa
    fc = 4.434749378960987   # our level-matched Fc/Fc+ absolute (config/electrolyte.py)
    idx = load_index()
    rows = []
    for s in systems:
        if not (_done(s, "ox") and _done(s, "red")):
            continue
        ro = json.loads((CALC / str(s) / "ox" / "result.json").read_text())
        rr = json.loads((CALC / str(s) / "red" / "result.json").read_text())
        Go = ro["e_smd_eV"] + (ro.get("g_thermal_eV") or 0.0)
        Gr = rr["e_smd_eV"] + (rr.get("g_thermal_eV") or 0.0)
        E_abs = -(Gr - Go)          # O + e- -> R
        E_vs_Fc = E_abs - fc
        rec = idx[s]
        rows.append(dict(sys=s, charge_ox=rec["charge_ox"], exp=rec["exp"],
                         calc=round(E_vs_Fc, 3), imp_dft=round(rec["imp_dft"], 3),
                         err=round(E_vs_Fc - rec["exp"], 3),
                         err_impdft=round(rec["imp_dft"] - rec["exp"], 3)))
    df = pd.DataFrame(rows)
    if df.empty:
        print("no completed systems yet"); return
    (ROOT / "results").mkdir(exist_ok=True)
    df.to_csv(ROOT / "results" / "orop_benchmark.csv", index=False)
    from scipy.stats import spearmanr, kendalltau
    mae = df.err.abs().mean(); rmse = (df.err**2).mean()**0.5
    rho = spearmanr(df.calc, df.exp).correlation if len(df) > 2 else float("nan")
    mae_imp = df.err_impdft.abs().mean()
    print(df.sort_values(["charge_ox", "exp"]).to_string(index=False))
    print(f"\n=== GLOBAL (n={len(df)}) ===")
    print(f"  OUR   MAE={mae:.3f}  RMSE={rmse:.3f}  signed={df.err.mean():+.3f}  Spearman(all)={rho:.3f}")
    print(f"  OROP raw implicit-DFT MAE (same systems)={mae_imp:.3f}")
    # Within-class stats: ranking is Tier-1's real job; a constant per-class offset does
    # not hurt ranking, so within-class Spearman/Kendall is the metric that matters.
    print(f"\n=== BY CHARGE CLASS (within-class ranking is the Tier-1 metric) ===")
    print(f"  {'q_ox':>4s} {'n':>3s} {'MAE':>6s} {'signed':>7s} {'Spearman':>9s} {'Kendall':>8s}")
    for q, sub in df.groupby("charge_ox"):
        if len(sub) >= 3:
            sr = spearmanr(sub.calc, sub.exp).correlation
            kt = kendalltau(sub.calc, sub.exp).correlation
            sr_s, kt_s = f"{sr:9.3f}", f"{kt:8.3f}"
        else:
            sr_s, kt_s = f"{'n<3':>9s}", f"{'n<3':>8s}"
        print(f"  {q:+4d} {len(sub):3d} {sub.err.abs().mean():6.3f} {sub.err.mean():+7.3f} {sr_s} {kt_s}")
    print(f"\nwrote results/orop_benchmark.csv")


DEFAULT = [103, 96, 101, 90, 16, 115,       # charge_ox=+1 (radical-cation / neutral)
           104, 105, 109, 108, 171, 107,    # charge_ox=0  (neutral / radical-anion)
           192, 191]                         # charge_ox=+2 and -1 (multiply charged)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", type=int)
    ap.add_argument("--fanout", type=str, help="comma list, or 'default'")
    ap.add_argument("--aggregate", type=str, help="comma list, or 'default'")
    ap.add_argument("--backend", default="gpu")
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    def parse(x):
        if x == "default":
            return DEFAULT
        if x == "all":   # every system with a calcs/orop/<s> dir
            return sorted(int(p.name) for p in CALC.iterdir()
                          if p.is_dir() and p.name.isdigit())
        return [int(v) for v in x.split(",")]

    if a.only is not None:
        run_system(a.only, a.backend, a.force)
    elif a.fanout:
        fanout(parse(a.fanout), a.force)
    elif a.aggregate:
        aggregate(parse(a.aggregate))
    else:
        print("nothing to do; use --only / --fanout / --aggregate")


if __name__ == "__main__":
    main()
