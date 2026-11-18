#!/usr/bin/env python
"""Recompute the main-text numbers that come from additional analyses of the stored records.

CPU only; runs no chemistry. Writes manuscript/evidence/2026-11-18_cpu_checks/supporting_stats.json.

  * conformer seeds: per-state seed counts and seed spreads (pool and all stored records)
  * outer-sphere check: IEF-PCM / Born ratio for pool molecules with both values
  * basis-function counts (def2-TZVPD) of the DMPZ molecules in the throughput table
  * UMA surrogate test: per-couple-type linear map from UMA gas-phase energy differences to
    DFT SMD electronic potentials (leave-one-out), rank correlations, and the anolyte
    pre-screen funnel against the Figure 5 front

  OMP_NUM_THREADS=2 PYTHONPATH=src python scripts/validation/manuscript_supporting_stats.py
"""
from __future__ import annotations

import csv
import glob
import json
import os
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts" / "plotting" / "manuscript"))
from redox.core.io import atomic_write_json  # noqa: E402

CALCS = ROOT / "calcs"
TABLE = ROOT / "results/validation/full_pool_scored_legacyG_MeCN_TBAPF6.csv"
OUT = ROOT / "manuscript/evidence/2026-11-18_cpu_checks/supporting_stats.json"


def load(p):
    return json.loads(Path(p).read_text()) if Path(p).exists() else None


def seeds():
    pool = [r["id"] for r in csv.DictReader(open(TABLE))]
    states, multi = [], []
    for g in pool:
        for f in glob.glob(str(CALCS / f"dft/{g}/*/result.json")):
            d = load(f)
            states.append(d.get("n_seeds"))
            if d.get("n_seeds") and d["n_seeds"] > 1 and d.get("seed_spread_eV") is not None:
                multi.append(d["seed_spread_eV"])
    allr, pen, failed = 0, [], 0
    for f in glob.glob(str(CALCS / "dft/*/*/result.json")):
        d = load(f)
        e = d.get("seed_e_smd_eV")
        if not e or len(e) < 2 or d.get("winning_seed") is None:
            continue
        allr += 1
        v = [x for x in e if x is not None]
        if e[0] is None:
            failed += 1
        else:
            pen.append(e[0] - min(v))
    return dict(pool_states=len(states), pool_three_seed=sum(s == 3 for s in states),
                pool_multi_seed=len(multi), spread_le_5meV=sum(s <= 0.005 for s in multi),
                spread_gt_25meV=sum(s > 0.025 for s in multi), spread_gt_100meV=sum(s > 0.1 for s in multi),
                all_records_multi_seed=allr, first_seed_gt_25meV=sum(p > 0.025 for p in pen),
                first_seed_failed=failed)


def pcm_ratio():
    pool = {r["id"] for r in csv.DictReader(open(TABLE))}
    src = Path("/nfs/lambda_stor_01/homes/rzhu/archive/0_redox/results/lambda_outer.csv")
    rows = [r for r in csv.DictReader(open(src)) if r["id"] in pool and r.get("lambda_o_pcm_eV")]
    rat = [float(r["pcm_over_born"]) for r in rows]
    return dict(source=str(src), couples=len(rows), molecules=len({r["id"] for r in rows}),
                ratio_min=min(rat), ratio_max=max(rat))


def basis_counts():
    from pyscf import gto
    from ase.io import read
    out = {}
    for gid in ["exp_dmpz_mb_methyl", "exp_dmpz_mb_ethyl", "hyp_dmpz_mb_trifluoroethyl",
                "hyp_dmpz_mb_methoxyethyl", "hyp_dmpz_27dimethyl_ethyl"]:
        a = read(str(CALCS / f"dft/{gid}/neu/opt.xyz"))
        s = "\n".join(f"{x} {p[0]} {p[1]} {p[2]}" for x, p in zip(a.get_chemical_symbols(), a.positions))
        out[gid] = dict(atoms=len(a), nao_def2_tzvpd=gto.M(atom=s, basis="def2-tzvpd", verbose=0).nao_nr())
    return out


def surrogate():
    import fig5_pareto as fp
    rows = list(csv.DictReader(open(TABLE)))
    pairs = [("ox1", "neu"), ("ox2", "ox1"), ("neu", "red1"), ("red1", "red2")]
    data = []
    for r in rows:
        gid = r["id"]
        for ox, red in pairs:
            dO, dR = load(CALCS / f"dft/{gid}/{ox}/result.json"), load(CALCS / f"dft/{gid}/{red}/result.json")
            uO, uR = load(CALCS / f"uma/{gid}/{ox}/result.json"), load(CALCS / f"uma/{gid}/{red}/result.json")
            if not all((dO, dR, uO, uR)) or dO.get("e_smd_eV") is None or dR.get("e_smd_eV") is None:
                continue
            data.append((gid, f"{ox}/{red}", -(dR["e_smd_eV"] - dO["e_smd_eV"]),
                         -(float(uR["energy_eV"]) - float(uO["energy_eV"]))))
    per_class, P, D, loo = {}, [], [], {}
    for cls in sorted({x[1] for x in data}):
        R = [x for x in data if x[1] == cls]
        d = np.array([x[2] for x in R]); u = np.array([x[3] for x in R])
        pred = np.array([np.polyval(np.polyfit(np.delete(u, i), np.delete(d, i), 1), u[i]) for i in range(len(R))])
        per_class[cls] = dict(n=len(R), loo_MAE_V=float(np.mean(abs(pred - d))),
                              spearman=float(spearmanr(u, d)[0]), dft_sd_V=float(d.std()))
        P += list(pred); D += list(d)
        for (gid, *_), pv in zip(R, pred):
            loo.setdefault(gid, {})[cls] = float(pv)
    P, D = np.array(P), np.array(D)
    # anolyte funnel vs the Figure 5 front (two-electron, non-control anolytes of the 91 pool)
    frozen, _ = fp.load_frozen()
    ano = [r for r in frozen if r["pool"] == "anolyte" and not r["ctrl"]]
    front = {r["id"] for r in fp.pareto(ano, True)}
    have = [r for r in ano if {"neu/red1", "red1/red2"} <= set(loo.get(r["id"], {}))]
    for r in have:
        r["Es"] = float(np.mean([loo[r["id"]]["neu/red1"], loo[r["id"]]["red1/red2"]]))
    funnel = {}
    for m in (0.15, 0.3, 0.45):
        keep = [a["id"] for a in have if not any(b is not a and b["cap"] >= a["cap"] and b["Es"] + m < a["Es"] for b in have)]
        funnel[str(m)] = dict(sent_to_dft=len(keep), front_recovered=len(front & set(keep)))
    return dict(couples=len(D), molecules=len({x[0] for x in data}), loo_MAE_V=float(np.mean(abs(P - D))),
                per_class=per_class, anolyte_two_electron=len(ano), anolyte_with_both_couples=len(have),
                anolyte_front=len(front), funnel=funnel)


def main():
    rec = dict(seeds=seeds(), outer_sphere_check=pcm_ratio(), basis_counts=basis_counts(),
               uma_surrogate=surrogate())
    OUT.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_json(OUT, rec)
    print(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
