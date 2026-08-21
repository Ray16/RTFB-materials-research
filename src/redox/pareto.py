"""Selection engine: turn the per-candidate scorecard into a discovery shortlist.

Applies the design decided earlier:
  - split into ANOLYTE / CATHOLYTE pools (redox potential is not 'more is better' globally;
    an ambipolar molecule enters BOTH pools on its respective side);
  - objectives = the TRUSTWORTHY axes only, oriented so higher = better:
        anolyte  : voltage = -E_anolyte  | catholyte: voltage = +E_catholyte
        capacity = specific capacity      (both)
        stability = dG_disp               (both; missing for single-wave -> incomparable)
        kinetics = -lambda                (both)
  - SIGMA-AWARE domination: A dominates B only if A is not worse than B beyond combined noise
    on every shared objective AND strictly better (beyond noise) on at least one. sigma per
    axis from config/scorecard_config.py; capacity sigma = 0 (exact).
  - proxies (solubility, SA) are NOT objectives; carried as annotations/tiebreakers.
  - the Pareto front is the primary output; a transparent normalized figure-of-merit gives a
    secondary headline ranking.

  PYTHONPATH=src python -m redox.pareto
"""
from __future__ import annotations
import csv
import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results"

# figure-of-merit weights (transparent; document any change). Sum ~1.
FOM_WEIGHTS = dict(voltage=0.30, capacity=0.30, stability=0.20, kinetics=0.20)


def _cfg(mod):
    spec = importlib.util.spec_from_file_location(mod, ROOT / "config" / f"{mod}.py")
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _load_candidates():
    p = RESULTS / "scorecard.csv"
    rows = []
    with p.open() as f:
        for r in csv.DictReader(f):
            if r.get("status") == "candidate":
                rows.append(r)
    return rows


def _objectives(cand, pool, sc):
    """Return {name: (value_higher_is_better, sigma)} for a candidate in a pool.
    Missing objective -> omitted (incomparable on that axis)."""
    o = {}
    E = _f(cand["E_anolyte_V"]) if pool == "anolyte" else _f(cand["E_catholyte_V"])
    if E is not None:
        o["voltage"] = ((-E if pool == "anolyte" else E), _f(cand["sigma_E_V"]) or 0.0)
    cap = _f(cand["specific_capacity_mAh_g"])
    if cap is not None:
        o["capacity"] = (cap, sc.SIGMA_CAPACITY)               # exact
    dg = _f(cand["dG_disp_kJmol"])
    if dg is not None:
        o["stability"] = (dg, sc.SIGMA_DISP_EV * 96.485)       # eV sigma -> kJ/mol
    lam = _f(cand["lambda_eV"])
    if lam is not None:
        o["kinetics"] = (-lam, sc.SIGMA_LAMBDA_EV)
    return o


def _dominates(A, B, k=1.0):
    """sigma-aware: A dominates B iff, over objectives they SHARE, A is never worse beyond
    combined noise and is better beyond noise on >=1. Returns False if they share <2 axes."""
    shared = set(A) & set(B)
    if len(shared) < 2:
        return False
    better = False
    for o in shared:
        va, sa = A[o]; vb, sb = B[o]
        tol = k * (sa ** 2 + sb ** 2) ** 0.5
        if vb - va > tol:      # A worse than B beyond noise
            return False
        if va - vb > tol:      # A better than B beyond noise
            better = True
    return better


def _pareto_front(cands, pool, sc):
    objs = {c["id"]: _objectives(c, pool, sc) for c in cands}
    front = []
    for c in cands:
        cid = c["id"]
        dominated = any(_dominates(objs[o["id"]], objs[cid])
                        for o in cands if o["id"] != cid)
        if not dominated:
            front.append(cid)
    return front, objs


def _fom(cands, pool, sc):
    """Min-max normalized weighted figure of merit within the pool (secondary ranking).
    Missing stability -> neutral 0.5 (neither rewarded nor penalized)."""
    objs = {c["id"]: _objectives(c, pool, sc) for c in cands}
    axes = ["voltage", "capacity", "stability", "kinetics"]
    ranges = {}
    for a in axes:
        vals = [objs[c["id"]][a][0] for c in cands if a in objs[c["id"]]]
        ranges[a] = (min(vals), max(vals)) if vals else None
    scores = {}
    for c in cands:
        cid = c["id"]; s = 0.0
        for a in axes:
            lo_hi = ranges[a]
            if a in objs[cid] and lo_hi and lo_hi[1] > lo_hi[0]:
                v = objs[cid][a][0]
                norm = (v - lo_hi[0]) / (lo_hi[1] - lo_hi[0])
            else:
                norm = 0.5      # missing or degenerate -> neutral
            s += FOM_WEIGHTS[a] * norm
        scores[cid] = round(s, 3)
    return scores


def run_pool(pool, cands, sc):
    # pool membership = the molecule has an accessible couple on this side
    key = "E_anolyte_V" if pool == "anolyte" else "E_catholyte_V"
    pool_cands = [c for c in cands if _f(c[key]) is not None]
    if not pool_cands:
        return []
    front, _ = _pareto_front(pool_cands, pool, sc)
    fom = _fom(pool_cands, pool, sc)
    out = []
    for c in pool_cands:
        out.append(dict(pool=pool, id=c["id"], family=c["family"],
                        E_V=_f(c[key]), n=c["n_accessible"],
                        capacity=_f(c["specific_capacity_mAh_g"]),
                        lambda_eV=_f(c["lambda_eV"]),
                        dG_disp_kJmol=_f(c["dG_disp_kJmol"]),
                        SA=_f(c["SA_score"]), dGsolv=_f(c["dGsolv_proxy_eV"]),
                        pareto_optimal=(c["id"] in front), fom=fom[c["id"]]))
    return sorted(out, key=lambda x: -x["fom"])


def main():
    sc = _cfg("scorecard_config")
    cands = _load_candidates()
    all_rows = []
    for pool in ("anolyte", "catholyte"):
        rows = run_pool(pool, cands, sc)
        all_rows.extend(rows)
        print(f"\n=== {pool.upper()} pool ({len(rows)} candidates) ===")
        print(f"{'id':22s} {'E(V)':>6s} {'n':>2s} {'Cap':>5s} {'lam':>5s} {'dGdisp':>7s} "
              f"{'SA':>4s} {'FoM':>5s} {'Pareto':>7s}")
        print("-" * 78)
        for r in rows:
            cap = f"{r['capacity']:.0f}" if r['capacity'] else "-"
            lam = f"{r['lambda_eV']:.2f}" if r['lambda_eV'] is not None else "-"
            dg = f"{r['dG_disp_kJmol']:.0f}" if r['dG_disp_kJmol'] is not None else "-"
            sa = f"{r['SA']:.1f}" if r['SA'] is not None else "-"
            star = "  YES" if r["pareto_optimal"] else ""
            print(f"{r['id']:22s} {r['E_V']:+6.2f} {r['n']:>2s} {cap:>5s} {lam:>5s} {dg:>7s} "
                  f"{sa:>4s} {r['fom']:5.2f} {star:>7s}")
        top = [r for r in rows if r["pareto_optimal"]]
        print(f"Pareto-optimal ({pool}): {[r['id'] for r in top]}")

    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / "pareto_shortlist.csv"
    cols = ["pool", "id", "family", "E_V", "n", "capacity", "lambda_eV", "dG_disp_kJmol",
            "SA", "dGsolv", "pareto_optimal", "fom"]
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(all_rows)
    print(f"\nObjectives = trustworthy axes (voltage, capacity, stability, kinetics); proxies "
          f"(SA, dGsolv) are annotations. sigma-aware domination. Weights {FOM_WEIGHTS}.")
    print(f"wrote {out}")
    return all_rows


if __name__ == "__main__":
    main()
