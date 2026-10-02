#!/usr/bin/env python
"""SCF ground-state + level-crossing check for the E° reference frame.

Two questions, one set of single points on the EXISTING optimized geometries (opt.xyz):

  1) Did production SCF land on the lowest-energy SCF solution?  Each state is re-solved from
     several initial guesses (minao / atom / huckel, plus — for open shells — orbital-occupation
     guesses built from the closed-shell N+1 system minus one electron from HOMO-k, and from the
     N-1 system plus one electron in LUMO+k). The lowest converged energy found is taken as the
     ground-state estimate (variational; more guesses can only lower it); production is suspect
     if it sits measurably above it. Motivated by ferrocenium: production E(Fc+) - E(Fc) in gas,
     each at its SMD-optimized geometry, is 8.97 eV vs the NIST evaluated ionization energy
     6.71 +/- 0.08 eV (webbook.nist.gov, CAS 102-54-5); its SMD dG_solv is -4.92 eV vs a Born
     estimate of -1.2 to -2.1 eV (SASA-sphere 5.84 A / vdW-volume-sphere 3.27 A radius, eps 37.5).

  2) Is the OROP offset a method/reference effect?  The same states are also run at OROP's level
     (B3LYP-D3(0)/6-31G*, SMD and C-PCM MeCN), so E vs Fc can be compared level by level for the
     in-house anchors, ferrocene and a sample of OROP systems (whose own B3LYP numbers we have).

One JSON per (level, id, state): calcs/validation/scf_level/<level>/<id>/<state>.json, holding
every guess's gas + solvated energy, <S^2>, convergence. Resumable. Aggregate with --aggregate.

  gpu_reserve run <idx> -- env OMP_NUM_THREADS=2 PYTHONPATH=src \\
      python scripts/validation/reference/scf_level_check.py --shard 4:0 --backend gpu
  python scripts/validation/reference/scf_level_check.py --list          # task list
  python scripts/validation/reference/scf_level_check.py --aggregate     # -> results/validation/
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from redox.core.common import HARTREE_EV  # noqa: E402
from redox.qm.dft import _free_gpu_pool  # noqa: E402

OUT = ROOT / "calcs" / "validation" / "scf_level"
RES = ROOT / "results" / "validation"

# MeCN static permittivity from the project config (config/project.json solvent.eps_r), so
# the C-PCM level uses the same dielectric as every other continuum quantity in the pipeline.
EPS_MECN = float(json.loads((ROOT / "config" / "project.json").read_text())["solvent"]["eps_r"])

LEVELS = {
    # production energy level (redox.core.protocol.ACTIVE_SP)
    "prod": dict(xc="wb97m-v", nlc="vv10", disp=None, basis="def2-tzvpd", solv="smd"),
    # OROP's level. data/raw/validation/OROP/raw_ferrocene-ref-values.txt states "Unless
    # specified functional is b3lyp-D3 and basis is 6-31gs"; the D3 DAMPING is not stated there,
    # so zero-damping D3 here is an ASSUMPTION (dispersion largely cancels in a 1e energy
    # difference at fixed geometry, but this is not verified). Run with our SMD and with C-PCM.
    "b3lyp_smd": dict(xc="b3lyp", nlc=None, disp="d3zero", basis="6-31g*", solv="smd"),
    "b3lyp_cpcm": dict(xc="b3lyp", nlc=None, disp="d3zero", basis="6-31g*", solv="cpcm",
                       eps=EPS_MECN),
}

# in-house anchors + ferrocene (ids under calcs/dft/) and their states
ANCHORS = {
    "ferrocene": ["neu", "ox"],
    "tempo_parent": ["ox", "rad"],
    "phenothiazine_parent": ["ox", "neu"],
    "anthraquinone_parent": ["neu", "red1", "red2"],
    "methyl_viologen": ["ox2", "ox1", "neu"],
    "methylpyridinium": ["ox", "red"],
}
# OROP sample (calcs/orop/<sys>/{ox,red}): spans +1/0 and 0/-1 classes and the error range
OROP_SAMPLE = [1, 2, 10, 25, 40, 60, 80, 100, 120, 150, 170, 190]


def tasks():
    out = []
    for lv in LEVELS:
        for gid, sts in ANCHORS.items():
            for st in sts:
                out.append((lv, gid, ROOT / "calcs" / "dft" / gid / st, st))
        for s in OROP_SAMPLE:
            for st in ("ox", "red"):
                d = ROOT / "calcs" / "orop" / str(s) / st
                if (d / "opt.xyz").exists():
                    out.append((lv, f"orop_{s}", d, st))
    return out


def _mol(xyz: Path, charge: int, spin: int, basis: str):
    from pyscf import gto
    from redox.qm.dft import _ecp_for
    lines = xyz.read_text().splitlines()
    n = int(lines[0])
    atom = "\n".join(lines[2:2 + n])
    syms = [l.split()[0] for l in lines[2:2 + n]]
    return gto.M(atom=atom, basis=basis, ecp=_ecp_for(syms, basis), charge=charge,
                 spin=spin, verbose=0)


def _mf(mol, lv, solvated, backend, restricted=None):
    from redox.qm.dft import _dft_module
    dft = _dft_module(backend)
    closed = mol.spin == 0 if restricted is None else restricted
    mf = (dft.RKS if closed else dft.UKS)(mol).density_fit()
    if solvated:
        if lv["solv"] == "smd":
            mf = mf.SMD()
            mf.with_solvent.solvent = "acetonitrile"
        else:
            mf = mf.PCM()
            mf.with_solvent.method = "C-PCM"
            mf.with_solvent.eps = lv["eps"]
    mf.xc = lv["xc"]
    if lv["disp"]:
        mf.disp = lv["disp"]
    if lv["nlc"]:
        mf.nlc = lv["nlc"]
    mf.max_cycle = 200
    mf.conv_tol = 1e-9
    return mf


def _solve(mf, dm0=None):
    e = float(mf.kernel(dm0=dm0))
    if not mf.converged:
        try:
            m2 = mf.newton()
            e = float(m2.kernel(dm0=mf.make_rdm1()))
            mf = m2
        except Exception as exc:  # pragma: no cover
            print(f"   [warn] newton failed: {exc}", flush=True)
    s2 = None
    if getattr(mf, "mol", None) is not None and mf.mol.spin:
        try:
            ss = mf.spin_square()
            s2 = float(ss[0] if isinstance(ss, (tuple, list)) else ss)
        except Exception:
            pass
    return e, bool(mf.converged), s2


def _occ_guesses(xyz, q, mult, lv, solvated, backend, n_k=3):
    """UKS density guesses for a doublet from closed-shell N+1 (remove from HOMO-k) and N-1
    (add to LUMO+k) references. Returns [(label, dm)]."""
    import numpy as np
    out = []
    if mult != 2:
        return out
    for dq, kind in ((-1, "rm"), (+1, "add")):
        try:
            ref = _mf(_mol(xyz, q + dq, 0, lv["basis"]), lv, solvated, backend)
            ref.kernel()
            if not ref.converged:
                continue
            C, occ = ref.mo_coeff, ref.mo_occ
            try:
                occ_np = np.asarray(occ.get())
            except AttributeError:
                occ_np = np.asarray(occ)
            nocc = int((occ_np > 0).sum())
            target = _mf(_mol(xyz, q, 1, lv["basis"]), lv, solvated, backend)
            for k in range(n_k):
                oa = occ_np / 2.0
                ob = occ_np / 2.0
                if kind == "rm":
                    h = nocc - 1 - k
                    if h < 0:
                        break
                    ob = ob.copy(); ob[h] = 0.0
                else:
                    l = nocc + k
                    if l >= len(oa):
                        break
                    oa = oa.copy(); oa[l] = 1.0
                xp = type(C)
                try:
                    import cupy
                    oa_, ob_ = cupy.asarray(oa), cupy.asarray(ob)
                except ImportError:
                    oa_, ob_ = oa, ob
                dm = target.make_rdm1((C, C), (oa_, ob_))
                try:
                    dm = dm.get()        # keep guesses on the host; free GPU between solves
                except AttributeError:
                    pass
                out.append((f"{kind}{k}", dm))
            del ref, target, C
            _free_gpu_pool()
        except Exception as exc:
            print(f"   [warn] occupation guess {kind} failed: {type(exc).__name__}: {exc}",
                  flush=True)
    return out


def _multi_guess(xyz: Path, q: int, mult: int, lv: dict, phases, backend, tag=""):
    """Solve one (geometry, charge, mult) from every guess, per phase. Returns {phase: {...}}."""
    out = {}
    for phase in phases:
        solvated = phase == "solv"
        guesses = []
        for key in ("minao", "atom", "huckel"):
            mf = _mf(_mol(xyz, q, mult - 1, lv["basis"]), lv, solvated, backend)
            mf.init_guess = key
            try:
                e, c, s2 = _solve(mf)
                guesses.append(dict(guess=key, e_Ha=e, converged=c, s2=s2))
            except Exception as exc:
                guesses.append(dict(guess=key, error=f"{type(exc).__name__}: {exc}"[:200]))
            del mf
            _free_gpu_pool()
        for lab, dm in _occ_guesses(xyz, q, mult, lv, solvated, backend):
            mf = _mf(_mol(xyz, q, mult - 1, lv["basis"]), lv, solvated, backend)
            try:
                e, c, s2 = _solve(mf, dm0=dm)
                guesses.append(dict(guess=lab, e_Ha=e, converged=c, s2=s2))
            except Exception as exc:
                guesses.append(dict(guess=lab, error=f"{type(exc).__name__}: {exc}"[:200]))
            del mf, dm
            _free_gpu_pool()
        ok = [g for g in guesses if g.get("converged")]
        best = min(ok, key=lambda g: g["e_Ha"]) if ok else None
        out[phase] = dict(guesses=guesses, best=best)
        print(f"  {tag} {phase}: " + ", ".join(
            f"{g['guess']}={g['e_Ha']:.6f}{'' if g.get('converged') else '(nc)'}"
            for g in guesses if "e_Ha" in g), flush=True)
    return out


def run_task(lv_name, gid, sdir: Path, st, backend):
    out_p = OUT / lv_name / gid / f"{st}.json"
    if out_p.exists():
        return "skip"
    res = json.loads((sdir / "result.json").read_text())
    q, mult = int(res["charge"]), int(res["mult"])
    xyz = sdir / "opt.xyz"
    lv = LEVELS[lv_name]
    rec = dict(level=lv_name, level_def=lv, id=gid, state=st, charge=q, mult=mult,
               geom=str(xyz), g_thermal_eV=res.get("g_thermal_eV"))
    t0 = time.time()
    rec["phases"] = _multi_guess(xyz, q, mult, lv, ("gas", "solv"), backend,
                                 tag=f"{lv_name} {gid}/{st}")
    rec["wall_s"] = round(time.time() - t0, 1)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    out_p.write_text(json.dumps(rec, indent=1))
    return "done"


# ------------------------------------------------------- candidate set (decision-feeding energies)
CAND_OUT = OUT / "candidates"


def candidate_tasks():
    """Every energy that feeds a ranked candidate's E deg or lambda, at the production level:
    each relaxed state (gas + SMD, the active-protocol record) and each current lambda cross
    point (gas). Molecule ids = those in results/scorecard.csv."""
    import csv
    from redox.core.common import read_result, state_names
    from redox.properties.reorg import _cache_path
    ids = sorted({r["id"] for r in csv.DictReader(open(ROOT / "results" / "scorecard.csv"))})
    out = []
    for gid in ids:
        sts = state_names(gid)
        qm = {}
        for st in sts:
            r = read_result(gid, st)
            qm[st] = (int(r["charge"]), int(r["mult"]))
            out.append(dict(kind="state", id=gid, name=st, q=qm[st][0], mult=qm[st][1],
                            xyz=ROOT / "calcs" / "dft" / gid / st / "opt.xyz",
                            phases=("gas", "solv"),
                            prod=dict(gas=r.get("e_gas_eV"), solv=r.get("e_smd_eV"))))
        for sp in sts:
            for at in sts:
                if sp == at or abs(qm[sp][0] - qm[at][0]) != 1:
                    continue
                c = _cache_path(gid, sp, at, *qm[sp])
                if not c.exists():
                    continue
                d = json.loads(c.read_text())
                out.append(dict(kind="cross", id=gid, name=f"{sp}_at_{at}", q=qm[sp][0],
                                mult=qm[sp][1],
                                xyz=ROOT / "calcs" / "dft" / gid / at / "opt.xyz",
                                phases=("gas",), prod=dict(gas=d.get("e_gas_eV")),
                                cache=str(c)))
    return out


def run_candidate(t, backend):
    out_p = CAND_OUT / t["id"] / f"{t['name']}.json"
    if out_p.exists():
        return "skip"
    lv = LEVELS["prod"]
    t0 = time.time()
    rec = dict(level="prod", level_def=lv, **{k: (str(v) if isinstance(v, Path) else v)
                                              for k, v in t.items() if k != "phases"})
    rec["phases"] = _multi_guess(t["xyz"], t["q"], t["mult"], lv, t["phases"], backend,
                                 tag=f"cand {t['id']}/{t['name']}")
    for ph, b in rec["phases"].items():
        pe = t["prod"].get(ph)
        b["production_minus_best_eV"] = (pe - b["best"]["e_Ha"] * HARTREE_EV
                                         if (b["best"] and pe is not None) else None)
    rec["wall_s"] = round(time.time() - t0, 1)
    out_p.parent.mkdir(parents=True, exist_ok=True)
    out_p.write_text(json.dumps(rec, indent=1))
    return "done"


def _thermal_and_prod(rec: dict):
    """(g_thermal_eV, production gas eV, production SMD eV) for the state a record was run on:
    thermal from result.json, else the OROP backfill thermal.json next to it; production
    energies from the ACTIVE-protocol record (redox.core.protocol.load_record)."""
    from redox.core.protocol import load_record
    sd = Path(rec["geom"]).parent
    res = json.loads((sd / "result.json").read_text())
    g = res.get("g_thermal_eV")
    if g is None and (sd / "thermal.json").exists():
        g = json.loads((sd / "thermal.json").read_text()).get("g_thermal_eV")
    pr = load_record(sd, int(rec["charge"]), int(rec["mult"])) or {}
    return g, pr.get("e_gas_eV"), pr.get("e_smd_eV")


# ------------------------------------------------------------------------------- aggregate
def aggregate():
    import pandas as pd
    sys.path.insert(0, str(ROOT / "config"))
    import validation as V
    sys.path.insert(0, str(ROOT / "scripts" / "validation" / "orop"))
    from run_orop_benchmark import load_index
    exp = {(v["id"], e["event"]): e["exp_V_vs_Fc"] for v in V.VALIDATION
           for e in v.get("events", []) if "exp_V_vs_Fc" in e}
    couples = {"ferrocene": [("ox", "neu", "ox->neu")],
               "tempo_parent": [("ox", "rad", "ox->rad")],
               "phenothiazine_parent": [("ox", "neu", "ox->neu")],
               "anthraquinone_parent": [("neu", "red1", "neu->red1"), ("red1", "red2", "red1->red2")],
               "methyl_viologen": [("ox2", "ox1", "ox2->ox1"), ("ox1", "neu", "ox1->neu")],
               "methylpyridinium": [("ox", "red", "ox->red")]}
    orop = load_index()

    def load(lv, gid, st):
        p = OUT / lv / gid / f"{st}.json"
        if not p.exists():
            return None
        r = json.loads(p.read_text())
        return r if r["level_def"] == LEVELS[lv] else None

    # production energies as actually used by the pipeline, for the "is production on the
    # ground state" comparison
    scf_rows, e_rows = [], []
    for lv in LEVELS:
        for p in sorted((OUT / lv).glob("*/*.json")):
            r = json.loads(p.read_text())
            if r["level_def"] != LEVELS[lv]:
                print(f"[stale] {p} level_def {r['level_def']} != {LEVELS[lv]} — rerun it")
                continue
            row = dict(level=lv, id=r["id"], state=r["state"], charge=r["charge"], mult=r["mult"])
            for ph in ("gas", "solv"):
                b = r["phases"][ph]["best"]
                g0 = next((g for g in r["phases"][ph]["guesses"] if g["guess"] == "minao"), {})
                row[f"{ph}_best_guess"] = b and b["guess"]
                row[f"{ph}_minao_above_best_eV"] = (
                    (g0["e_Ha"] - b["e_Ha"]) * HARTREE_EV if b and g0.get("converged") else None)
                row[f"{ph}_best_s2"] = b and b.get("s2")
            if lv == "prod":
                _, pg, ps = _thermal_and_prod(r)
                for ph, pe in (("gas", pg), ("solv", ps)):
                    b = r["phases"][ph]["best"]
                    if b and pe is not None:
                        row[f"{ph}_production_above_best_eV"] = pe - b["e_Ha"] * HARTREE_EV
            scf_rows.append(row)

        def G(gid, st, ph="solv"):
            r = load(lv, gid, st)
            if not r or not r["phases"][ph]["best"]:
                return None
            g = _thermal_and_prod(r)[0]
            if g is None:
                return None
            return r["phases"][ph]["best"]["e_Ha"] * HARTREE_EV + g

        g_ox, g_neu = G("ferrocene", "ox"), G("ferrocene", "neu")
        fc_abs = (g_ox - g_neu) if None not in (g_ox, g_neu) else None
        r_ox, r_neu = load(lv, "ferrocene", "ox"), load(lv, "ferrocene", "neu")
        fc_ip = None
        if r_ox and r_neu and r_ox["phases"]["gas"]["best"] and r_neu["phases"]["gas"]["best"]:
            fc_ip = (r_ox["phases"]["gas"]["best"]["e_Ha"]
                     - r_neu["phases"]["gas"]["best"]["e_Ha"]) * HARTREE_EV
        e_rows.append(dict(level=lv, id="ferrocene", event="Fc_abs_V", value=fc_abs))
        e_rows.append(dict(level=lv, id="ferrocene", event="gas_dE_ox_minus_neu_eV (NIST IE 6.71+/-0.08)", value=fc_ip))
        for gid, cps in couples.items():
            for o, rd, ev in cps:
                go, gr = G(gid, o), G(gid, rd)
                if None in (go, gr, fc_abs):
                    continue
                E = go - gr - fc_abs
                e_rows.append(dict(level=lv, id=gid, event=ev, value=E, exp=exp.get((gid, ev)),
                                   err=(E - exp[(gid, ev)]) if (gid, ev) in exp else None))
        for s in OROP_SAMPLE:
            go, gr = G(f"orop_{s}", "ox"), G(f"orop_{s}", "red")
            if None in (go, gr, fc_abs) or s not in orop:
                continue
            E = go - gr - fc_abs
            e_rows.append(dict(level=lv, id=f"orop_{s}", event=f"q_ox={orop[s]['charge_ox']}",
                               value=E, exp=orop[s]["exp"], err=E - orop[s]["exp"],
                               orop_b3lyp=orop[s]["imp_dft"]))
    RES.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(scf_rows).to_csv(RES / "scf_ground_state_check.csv", index=False)
    # candidate set: the CURRENT production energy (active-protocol record / current lambda
    # cross-point cache — not the value stored when the check ran) vs the lowest solution found
    from redox.core.protocol import load_record
    cand = []
    for p in sorted(CAND_OUT.glob("*/*.json")):
        r = json.loads(p.read_text())
        if r["kind"] == "state":
            rec = load_record(Path(r["xyz"]).parent, r["q"], r["mult"]) or {}
            now = dict(gas=rec.get("e_gas_eV"), solv=rec.get("e_smd_eV"))
        else:
            c = Path(r["cache"])
            now = dict(gas=json.loads(c.read_text()).get("e_gas_eV") if c.exists() else None)
        for ph, b in r["phases"].items():
            pe = now.get(ph)
            cand.append(dict(id=r["id"], name=r["name"], kind=r["kind"], charge=r["q"],
                             mult=r["mult"], phase=ph, n_guesses=len(b["guesses"]),
                             n_converged=sum(1 for g in b["guesses"] if g.get("converged")),
                             best_guess=b["best"] and b["best"]["guess"],
                             production_minus_lowest_eV=(pe - b["best"]["e_Ha"] * HARTREE_EV
                                                         if (b["best"] and pe is not None)
                                                         else None),
                             at_check_time_eV=b.get("production_minus_best_eV")))
    pd.DataFrame(cand).to_csv(RES / "scf_candidates_check.csv", index=False)
    df = pd.DataFrame(e_rows)
    df.to_csv(RES / "level_crossing_Efc.csv", index=False)
    pd.set_option("display.width", 200)
    print(pd.DataFrame(scf_rows).round(3).to_string())
    print(df.round(3).to_string())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--shard", default="1:0", help="N:I — run every N-th task starting at I")
    ap.add_argument("--backend", default="gpu")
    ap.add_argument("--only", help="comma list of level:id:state")
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--aggregate", action="store_true")
    ap.add_argument("--reverse", action="store_true",
                    help="walk the task list backwards (to add workers alongside forward ones)")
    ap.add_argument("--claim", action="store_true",
                    help="atomic mkdir claim per task under <out>/_claims (skip claimed tasks)")
    ap.add_argument("--set", default="reference", choices=["reference", "candidates"],
                    help="reference = anchors/Fc/OROP at 3 levels; candidates = every energy "
                         "feeding a ranked candidate, production level")
    a = ap.parse_args()
    if a.aggregate:
        return aggregate()
    if a.set == "candidates":
        T = candidate_tasks()
        if a.list:
            for t in T:
                print(t["kind"], t["id"], t["name"], t["q"], t["mult"])
            print(len(T), "tasks")
            return
        n, i = map(int, a.shard.split(":"))
        order = list(enumerate(T))[::-1] if a.reverse else list(enumerate(T))
        for k, t in order:
            if k % n != i:
                continue
            if a.claim:
                try:
                    (CAND_OUT / "_claims").mkdir(parents=True, exist_ok=True)
                    (CAND_OUT / "_claims" / f"{t['id']}__{t['name']}").mkdir()
                except FileExistsError:
                    continue
            print(f"[{time.strftime('%H:%M:%S')}] cand {t['id']}/{t['name']}", flush=True)
            try:
                print("   ", run_candidate(t, a.backend), flush=True)
            except Exception as exc:
                print(f"   [FAIL] {type(exc).__name__}: {exc}", flush=True)
        return
    T = tasks()
    if a.only:
        want = {tuple(x.split(":")) for x in a.only.split(",")}
        T = [t for t in T if (t[0], t[1], t[3]) in want]
    if a.list:
        for t in T:
            print(t[0], t[1], t[3])
        print(len(T), "tasks")
        return
    n, i = map(int, a.shard.split(":"))
    order = list(enumerate(T))[::-1] if a.reverse else list(enumerate(T))
    for k, (lv, gid, sd, st) in order:
        if k % n != i:
            continue
        if (OUT / lv / gid / f"{st}.json").exists():
            continue
        if a.claim:
            try:
                (OUT / "_claims").mkdir(parents=True, exist_ok=True)
                (OUT / "_claims" / f"{lv}__{gid}__{st}").mkdir()
            except FileExistsError:
                continue
        print(f"[{time.strftime('%H:%M:%S')}] {lv} {gid}/{st}", flush=True)
        try:
            print("   ", run_task(lv, gid, sd, st, a.backend), flush=True)
        except Exception as exc:
            print(f"   [FAIL] {type(exc).__name__}: {exc}", flush=True)


if __name__ == "__main__":
    main()
