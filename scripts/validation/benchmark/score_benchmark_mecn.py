#!/usr/bin/env python
"""Score the production pipeline against the sourced MeCN benchmark (config/benchmark_mecn.py).

Reads results/redox_potentials.csv (computed E vs a live level-matched Fc) and writes
results/validation/benchmark_mecn.csv — one row per (molecule, couple, source) with the
experimental value, computed value, residual, tier and citation — plus summary statistics per
tier and per family. Wave spacings (E1 - E2, reference-free) are scored from the same table.
Nothing is fitted or corrected; excluded waves are listed with their reason, never scored.

  PYTHONPATH=src python scripts/validation/benchmark/score_benchmark_mecn.py
"""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "src"))
from redox.core.common import RESULTS, load_config  # noqa: E402
from redox.screening.scorecard import _bench_family  # noqa: E402

OUT = RESULTS / "validation" / "benchmark_mecn.csv"


def _f(x):
    try:
        v = float(x)
        return None if math.isnan(v) else v
    except (TypeError, ValueError):
        return None


def stats(res):
    n = len(res)
    if not n:
        return dict(n=0)
    mae = sum(abs(r) for r in res) / n
    rmse = math.sqrt(sum(r * r for r in res) / n)
    bias = sum(res) / n
    sd = math.sqrt(sum((r - bias) ** 2 for r in res) / (n - 1)) if n > 1 else float("nan")
    return dict(n=n, MAE=mae, RMSE=rmse, bias=bias, SD=sd)


def main():
    comp = {(r["id"], r["event"]): r for r in csv.DictReader(open(RESULTS / "redox_potentials.csv"))}
    rows = []
    for b in load_config("benchmark_mecn").BENCHMARK:
        fam = _bench_family(b["bench_family"])
        for ev in b["events"]:
            base = dict(id=b["id"], name=b["name"], family=fam, tier=ev["tier"],
                        tier_reason=ev["tier_reason"], source=ev["source"], doi=ev["doi"],
                        where=ev["where"], excluded_molecule=b.get("excluded") or "")
            if ev["kind"] == "E":
                c = comp.get((b["id"], ev["event"]), {})
                calc = _f(c.get("E_vs_Fc_V"))
                rows.append(dict(base, quantity=f"E[{ev['event']}]", exp_V=ev["exp_V_vs_Fc"],
                                 calc_V=calc, status=c.get("status", "not computed"),
                                 err_V=(calc - ev["exp_V_vs_Fc"]) if calc is not None else None))
            else:   # wave spacing E1 - E2 from the two computed couples of the ladder
                evs = [e["event"] for e in b["events"] if e["kind"] == "E" and e["source"] == ev["source"]]
                e1 = _f(comp.get((b["id"], evs[0]), {}).get("E_vs_Fc_V")) if len(evs) > 0 else None
                e2 = _f(comp.get((b["id"], evs[1]), {}).get("E_vs_Fc_V")) if len(evs) > 1 else None
                calc = (e1 - e2) if None not in (e1, e2) else None
                rows.append(dict(base, quantity="dE12", exp_V=ev["exp_V"], calc_V=calc,
                                 status="ok" if calc is not None else "INCOMPLETE",
                                 err_V=(calc - ev["exp_V"]) if calc is not None else None))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    cols = ["id", "name", "family", "quantity", "exp_V", "calc_V", "err_V", "status", "tier",
            "tier_reason", "source", "doi", "where", "excluded_molecule"]
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader()
        for r in rows:
            w.writerow({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})

    scored = [r for r in rows if r["err_V"] is not None and not r["excluded_molecule"]]
    print(f"{len(rows)} benchmark rows; {len(scored)} scored "
          f"({sum(1 for r in rows if r['err_V'] is None)} without a computed value)")
    for qty, sel in (("E (absolute, vs Fc)", lambda r: r["quantity"].startswith("E[")),
                     ("dE12 (reference-free spacing)", lambda r: r["quantity"] == "dE12")):
        for tier in ("A", "B"):
            s = stats([r["err_V"] for r in scored if sel(r) and r["tier"] == tier])
            if s["n"]:
                print(f"  {qty:30s} tier {tier}: n={s['n']:2d}  MAE={s['MAE']:.3f}  "
                      f"RMSE={s['RMSE']:.3f}  bias={s['bias']:+.3f}  SD={s['SD']:.3f} V")
        for fam in sorted({r["family"] for r in scored if sel(r) and r["tier"] == "A"} - {None}):
            s = stats([r["err_V"] for r in scored if sel(r) and r["tier"] == "A" and r["family"] == fam])
            print(f"      tier A {fam:22s} n={s['n']:2d}  MAE={s['MAE']:.3f}  bias={s['bias']:+.3f} V")
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
