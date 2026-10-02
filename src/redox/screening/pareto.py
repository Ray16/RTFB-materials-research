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
    axis = each row's OWN scorecard sigma column; capacity sigma = 0 (bookkeeping).
  - COMPLETENESS: a row missing any primary objective (voltage, capacity, kinetics, and
    stability when n>=2) is INCOMPLETE — excluded from the primary front and unable to
    dominate; it is reported on a separate partial-information list.
  - proxies (solubility, SA) are NOT objectives; carried as annotations/tiebreakers.
  - the sigma-aware Pareto front is the output; within a pool candidates are ordered by the
    primary axis (voltage). There is deliberately NO scalarized figure-of-merit: a weighted-sum
    FoM would (a) require arbitrary weights, (b) ignore sigma (so it disagrees with the
    sigma-aware front), and (c) be pool-relative via min-max normalization, so it silently
    rescales when the candidate set changes. Read the front + the raw axes instead.

  PYTHONPATH=src python -m redox.screening.pareto
"""
from __future__ import annotations
import csv

from redox.core.common import RESULTS
from redox.core.common import load_config as _cfg
from redox.core.common import to_float as _f


def _load_candidates():
    p = RESULTS / "scorecard.csv"
    rows = []
    with p.open() as f:
        for r in csv.DictReader(f):
            if r.get("status") == "candidate" and str(r.get("rankable", "True")) != "False":
                rows.append(r)
    return rows


PRIMARY = ("voltage", "capacity", "stability", "kinetics")


def _objectives(cand, pool, sc, cap_key="specific_capacity_mAh_g"):
    """Return ({name: (value_higher_is_better, sigma)}, missing) for one pool row.

    Every sigma is the candidate's OWN reported uncertainty (scorecard columns sigma_E_V,
    sigma_lambda_eV, sigma_disp_eV — derived from the benchmarks), not a config constant.
    `missing` lists PRIMARY objectives that are absent. Stability is NOT APPLICABLE (not
    missing) for a 1-electron path: there is no intermediate to disproportionate."""
    o, missing = {}, []
    E = _f(cand.get("E_V"))
    if E is not None:
        o["voltage"] = ((-E if pool == "anolyte" else E), _f(cand["sigma_E_V"]) or 0.0)
    else:
        missing.append("voltage")
    cap = _f(cand.get(cap_key))
    if cap is not None:
        o["capacity"] = (cap, sc.SIGMA_CAPACITY)               # bookkeeping on the path
    else:
        missing.append("capacity")
    dg = _f(cand.get("dG_disp_kJmol"))
    if dg is not None:
        o["stability"] = (dg, (_f(cand.get("sigma_disp_eV")) or sc.SIGMA_DISP_EV) * 96.485)
    elif str(cand.get("disp_applicable")) == "True":
        missing.append("stability")
    # Only QC-clean inner-sphere lambda is admitted. A row whose couples were all flagged or
    # not computed carries no kinetics value and is therefore INCOMPLETE, not "unranked".
    lam = _f(cand.get("lambda_i_eV"))
    if lam is not None:
        o["kinetics"] = (-lam, _f(cand.get("sigma_lambda_eV")) or sc.SIGMA_LAMBDA_EV)
    else:
        missing.append("kinetics")
    return o, missing


def _dominates(A, B, k=1.0):
    """sigma-aware: A dominates B iff, over objectives they SHARE, A is never worse beyond
    combined noise and is better beyond noise on >=1. Returns False if they share <2 axes.
    Callers guarantee A is COMPLETE (see _pareto_front), so a candidate can never gain
    dominance by lacking an axis."""
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


def _pareto_front(cands, pool, sc, cap_key="specific_capacity_mAh_g"):
    """Returns (front ids, objectives, dominated_by, missing).

    PRIMARY front: COMPLETE candidates only, dominated only by complete candidates.
    INCOMPLETE candidates (any primary objective missing) can be dominated by complete ones
    on the axes they share but can NEVER dominate anyone — they form a separate partial-
    information list (partial_front = not dominated by any complete candidate)."""
    objs, missing = {}, {}
    for c in cands:
        objs[c["id"]], missing[c["id"]] = _objectives(c, pool, sc, cap_key)
    complete = [c["id"] for c in cands if not missing[c["id"]]]
    front, dominated_by = [], {}
    for c in cands:
        cid = c["id"]
        dominated_by[cid] = [o for o in complete if o != cid and _dominates(objs[o], objs[cid])]
        if not missing[cid] and not dominated_by[cid]:
            front.append(cid)
    return front, objs, dominated_by, missing


def run_pool(pool, cands, sc, cap_key="specific_capacity_mAh_g", scenario=""):
    # pool membership = the molecule has an accessible couple on this side
    pool_cands = [c for c in cands if c.get("pool") == pool]
    if not pool_cands:
        return []
    front, _, dom_by, missing = _pareto_front(pool_cands, pool, sc, cap_key)
    out = []
    for c in pool_cands:
        out.append(dict(scenario=scenario, capacity_def=cap_key,
                        pool=pool, id=c["id"], family=c["family"],
                        E_V=_f(c["E_V"]), n=c["n_accessible"],
                        capacity=_f(c.get(cap_key)),
                        lambda_o_pcm_eV=_f(c.get("lambda_o_pcm_eV")),
                        lambda_het_eV=_f(c.get("lambda_het_eV")),
                        lambda_se_contact_eV=_f(c.get("lambda_se_contact_eV")),
                        dominated_by=";".join(dom_by[c["id"]]),
                        missing_objectives=";".join(missing[c["id"]]),
                        partial_front=(bool(missing[c["id"]]) and not dom_by[c["id"]]),
                        lambda_i_eV=_f(c["lambda_i_eV"]),
                        lambda_qc=c.get("lambda_qc", ""),
                        lambda_flags=c.get("lambda_flags", ""),
                        dG_disp_kJmol=_f(c["dG_disp_kJmol"]),
                        SA=_f(c["SA_score"]), dGsolv=_f(c["dGsolv_proxy_eV"]),
                        pareto_optimal=(c["id"] in front)))
    # order by the pool's primary axis (best voltage first): anolyte = most negative E,
    # catholyte = most positive E. The Pareto flag marks the shortlist; no scalarized score.
    return sorted(out, key=lambda x: (x["E_V"] if pool == "anolyte" else -x["E_V"]))


def main():
    sc = _cfg("scorecard_config")
    cands = _load_candidates()
    scenarios = getattr(sc, "CAPACITY_SCENARIOS", {"monomer": "specific_capacity_mAh_g"})
    all_rows = []
    fronts = {}
    for scen, cap_key in scenarios.items():
        print(f"\n######## capacity scenario: {scen}  ({cap_key}) ########")
        for pool in ("anolyte", "catholyte"):
            rows = run_pool(pool, cands, sc, cap_key, scen)
            all_rows.extend(rows)
            fronts[(scen, pool)] = {r["id"] for r in rows if r["pareto_optimal"]}
            _print_pool(pool, rows)
    # where does the shortlist depend on the supporting salt?
    names = list(scenarios)
    for pool in ("anolyte", "catholyte"):
        sets = [fronts.get((s, pool), set()) for s in names]
        if not any(sets):
            continue
        common = set.intersection(*sets)
        print(f"\n[{pool}] front in ALL salt scenarios: {sorted(common)}")
        for s, fs in zip(names, sets):
            print(f"[{pool}] only under {s}: {sorted(fs - common)}")
    _write(all_rows)
    return all_rows


def _print_pool(pool, rows):
    print(f"\n=== {pool.upper()} pool ({len(rows)} candidates) ===")
    print(f"{'id':22s} {'E(V)':>6s} {'n':>2s} {'Cap':>5s} {'lam':>5s} {'dGdisp':>7s} "
          f"{'SA':>4s} {'Pareto':>7s}  dominated by")
    print("-" * 90)
    for r in rows:
        cap = f"{r['capacity']:.0f}" if r['capacity'] else "-"
        lam = f"{r['lambda_i_eV']:.2f}" if r['lambda_i_eV'] is not None else "-"
        if r.get("lambda_qc") not in ("", "ok"):
            lam += "!"
        dg = f"{r['dG_disp_kJmol']:.0f}" if r['dG_disp_kJmol'] is not None else "-"
        sa = f"{r['SA']:.1f}" if r['SA'] is not None else "-"
        star = ("  YES" if r["pareto_optimal"] else
                (" part." if r["partial_front"] else ""))
        miss = f"  [INCOMPLETE: missing {r['missing_objectives']}]" if r["missing_objectives"] else ""
        print(f"{r['id']:22s} {r['E_V']:+6.2f} {r['n']:>2s} {cap:>5s} {lam:>5s} {dg:>7s} "
              f"{sa:>4s} {star:>7s}  {r['dominated_by']}{miss}")
    top = [r for r in rows if r["pareto_optimal"]]
    part = [r for r in rows if r["partial_front"]]
    print(f"Pareto-optimal ({pool}, complete data): {[r['id'] for r in top]}")
    if part:
        print(f"partial-information, not dominated ({pool}): {[r['id'] for r in part]}")


def _write(all_rows):
    RESULTS.mkdir(exist_ok=True)
    out = RESULTS / "pareto_shortlist.csv"
    cols = ["scenario", "capacity_def", "pool", "id", "family", "E_V", "n", "capacity",
            "lambda_i_eV", "lambda_qc", "lambda_flags", "lambda_o_pcm_eV", "lambda_het_eV",
            "lambda_se_contact_eV", "dG_disp_kJmol", "SA", "dGsolv", "pareto_optimal",
            "partial_front", "missing_objectives", "dominated_by"]
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(all_rows)
    print(f"\nObjectives = trustworthy axes (voltage, capacity, stability, kinetics=lambda_i); "
          f"lambda_o/lambda_het/lambda_se, SA, dGsolv are annotations. Rows missing a primary "
          f"objective are INCOMPLETE: never on the primary front, never dominate. Capacity run under each salt "
          f"scenario (max-load counter-ions). sigma-aware domination; no scalarized FoM.")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
