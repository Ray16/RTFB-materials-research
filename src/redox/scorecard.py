"""Unified per-candidate scorecard — the input contract for the selection engine.

Joins the five per-axis result tables (redox potentials, reorganization energy, reversibility,
disproportionation stability, capacity/proxies), applies the gating rules, rolls everything up
to ONE row per candidate molecule, and attaches a sigma (uncertainty) and trust level to every
axis so the Pareto engine can do sigma-aware domination.

GATING (config/electrolyte.py + config/scorecard_config.py):
  1. Reversibility HARD FILTER   : a couple counts only if verdict == reversible.
  2. Window gating               : a couple counts toward accessible n only if its E_vs_Fc lies
                                   inside WINDOW_V_VS_FC.  (A molecule with zero accessible
                                   reversible couples is rejected.)
  3. Anolyte/catholyte split     : by the accessible couples' potentials vs ANOLYTE_CATHOLYTE_
                                   DIVIDER_V (all below -> anolyte, all above -> catholyte, mixed
                                   -> ambipolar).

Data model: candidate = molecule; couple-level properties are rolled up:
  - E_operating : mean E_vs_Fc over accessible couples (the working potential)
  - n_accessible: # reversible in-window couples  ->  capacity = n*F/MW
  - lambda      : mean lambda over accessible couples (kinetics; lower = better)
  - dG_disp     : disproportionation dG of the accessible interior intermediate, if any
                  (worst/min over intermediates); n/a for a single-couple molecule
  - reversible  : all accessible couples reversible (True by construction of the filter)
  - SA, dGsolv  : molecule-level proxies

  PYTHONPATH=src python -m redox.scorecard
"""
from __future__ import annotations
import csv
import importlib.util
import statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"

# not real candidates: the Fc/Fc+ internal reference
REFERENCE_IDS = {"ferrocene"}


def _cfg(mod):
    spec = importlib.util.spec_from_file_location(mod, ROOT / "config" / f"{mod}.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def _load(name):
    p = RESULTS / f"{name}.csv"
    if not p.exists():
        return []
    with p.open() as f:
        return list(csv.DictReader(f))


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def build():
    ele = _cfg("electrolyte")
    sc = _cfg("scorecard_config")
    wlo, whi = ele.WINDOW_V_VS_FC
    divider = ele.ANOLYTE_CATHOLYTE_DIVIDER_V

    redox = _load("redox_potentials")
    lam = {(r["id"], r["couple"]): _f(r["lambda_i_eV"]) for r in _load("reorganization")}
    rev = {(r["id"], r["couple"]): r["verdict"] for r in _load("reversibility")}
    disp = {}
    for r in _load("stability_disproportionation"):
        disp.setdefault(r["id"], []).append(_f(r["dG_disp_kJmol"]))
    cap = {r["id"]: r for r in _load("capacity_and_proxies")}

    # group redox couples by molecule
    by_mol = {}
    for r in redox:
        if _f(r.get("E_vs_Fc_V")) is None:
            continue
        by_mol.setdefault(r["id"], {"name": r["name"], "family": r["family"], "couples": []})
        e = _f(r["E_vs_Fc_V"]); qo = int(r["q_ox"]); qr = int(r["q_red"])
        couple = r["event"]
        by_mol[r["id"]]["couples"].append(dict(
            couple=couple, E=e, q_ox=qo, q_red=qr,
            lam=lam.get((r["id"], couple)),
            reversible=(rev.get((r["id"], couple)) == "reversible"),
            in_window=(wlo <= e <= whi),
            sigma_E=sc.sigma_redox(qo, qr),
        ))

    rows = []
    for gid, m in by_mol.items():
        if gid in REFERENCE_IDS:        # Fc/Fc+ is the reference, not a candidate
            continue
        acc = [c for c in m["couples"] if c["reversible"] and c["in_window"]]
        rejected = [c for c in m["couples"] if not (c["reversible"] and c["in_window"])]
        if not acc:                     # GATE: no usable couple -> not a candidate
            rows.append(dict(id=gid, name=m["name"], family=m["family"],
                             status="REJECTED", reason="no reversible in-window couple",
                             n_accessible=0))
            continue
        Es = [c["E"] for c in acc]
        ano = [e for e in Es if e < divider]   # anolyte-side couples (low potential)
        cat = [e for e in Es if e >= divider]  # catholyte-side couples (high potential)
        role = "ambipolar" if (ano and cat) else ("anolyte" if ano else "catholyte")
        n_acc = len(acc)
        c = cap.get(gid, {})
        mw = _f(c.get("MW"))
        # capacity from ACCESSIBLE electrons (window+reversibility gated). NOTE: for an
        # ambipolar molecule this counts BOTH sides; on a single electrode only its own-side
        # electrons are usable (see per-side potentials).
        cspec = round(n_acc * ele.FARADAY / (mw * 3.6), 1) if mw else None
        lams = [c["lam"] for c in acc if c["lam"] is not None]
        disp_vals = disp.get(gid)          # interior intermediate(s), worst case = min
        rows.append(dict(
            id=gid, name=m["name"], family=m["family"], status="candidate", role=role,
            n_accessible=n_acc,
            E_anolyte_V=(round(statistics.mean(ano), 3) if ano else None),
            E_catholyte_V=(round(statistics.mean(cat), 3) if cat else None),
            sigma_E_V=round(max(c["sigma_E"] for c in acc), 3),
            specific_capacity_mAh_g=cspec, MW=(round(mw, 1) if mw else None),
            lambda_eV=(round(statistics.mean(lams), 3) if lams else None),
            sigma_lambda_eV=sc.SIGMA_LAMBDA_EV if lams else None,
            dG_disp_kJmol=(round(min(disp_vals), 1) if disp_vals else None),
            sigma_disp_eV=sc.SIGMA_DISP_EV if disp_vals else None,
            all_reversible=True,
            SA_score=_f(c.get("SA_score")),
            dGsolv_proxy_eV=_f(c.get("dGsolv_neutral_eV")),
            n_couples_rejected=len(rejected),
        ))

    # --- write ---
    RESULTS.mkdir(exist_ok=True)
    cols = ["id", "name", "family", "status", "role", "n_accessible", "E_anolyte_V",
            "E_catholyte_V", "sigma_E_V", "specific_capacity_mAh_g", "MW", "lambda_eV",
            "sigma_lambda_eV", "dG_disp_kJmol", "sigma_disp_eV", "all_reversible", "SA_score",
            "dGsolv_proxy_eV", "n_couples_rejected", "reason"]
    out = RESULTS / "scorecard.csv"
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k, "") for k in cols})

    cands = [r for r in rows if r["status"] == "candidate"]

    def _sortkey(r):
        return (r["role"], r["E_anolyte_V"] if r["E_anolyte_V"] is not None else r["E_catholyte_V"])

    print(f"{'id':22s} {'role':9s} {'n':>2s} {'E_an':>6s} {'E_cat':>6s} {'Cap':>5s} "
          f"{'lam':>5s} {'dGdisp':>7s} {'SA':>4s}")
    print("-" * 76)
    for r in sorted(cands, key=_sortkey):
        ea = f"{r['E_anolyte_V']:+.2f}" if r["E_anolyte_V"] is not None else "-"
        ec = f"{r['E_catholyte_V']:+.2f}" if r["E_catholyte_V"] is not None else "-"
        cap_s = f"{r['specific_capacity_mAh_g']:.0f}" if r["specific_capacity_mAh_g"] else "-"
        lam_s = f"{r['lambda_eV']:.2f}" if r["lambda_eV"] is not None else "-"
        dg_s = f"{r['dG_disp_kJmol']:.0f}" if r["dG_disp_kJmol"] is not None else "-"
        sa_s = f"{r['SA_score']:.1f}" if r["SA_score"] is not None else "-"
        print(f"{r['id']:22s} {r['role']:9s} {r['n_accessible']:2d} "
              f"{ea:>6s} {ec:>6s} {cap_s:>5s} {lam_s:>5s} {dg_s:>7s} {sa_s:>4s}")
    rej = [r for r in rows if r["status"] == "REJECTED"]
    if rej:
        print(f"\nREJECTED (gated out): {[r['id'] for r in rej]}")
    print(f"\n{len(cands)} candidates, {len(rej)} rejected. sigma + trust per axis attached "
          f"(config/scorecard_config.py). Window {wlo}..{whi} V vs Fc.")
    print(f"wrote {out}")
    return rows


if __name__ == "__main__":
    build()
