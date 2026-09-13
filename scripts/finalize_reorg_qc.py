#!/usr/bin/env python
"""Fold the reorg QC verdict into the production comparison table.

Merges the unbound-anion screen (screen/*.json) and the guarded recompute (recompute/*.json)
into results/d3tales_reorg_validation_prod/comparison.csv, adding:
  anion_homo_at_neu_eV : gas HOMO of the VERTICAL anion (>0 => unbound)
  anion_unbound        : True if the vertical radical anion is unbound
  lambda_cf_smd_eV     : conformer-matched + SMD recompute (only for the outliers we recomputed)
  reliable             : False if the inner-sphere lambda is untrustworthy

An inner-sphere lambda is UNRELIABLE when the vertical anion is unbound (QC flag) or the value is
implausibly large (>1.5 eV) — the D3TaLES thiosuccinimide outliers, which no recompute recipe
(conformer-matched, SMD, or both) rescues. These are excluded from the lambda population/figures.
Idempotent: re-running just refreshes the merged columns.
"""
from __future__ import annotations
import csv, json, statistics
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results" / "d3tales_reorg_validation_prod"
CSV = BASE / "comparison.csv"
UNRELIABLE_ABOVE = 1.5   # eV; an inner-sphere lambda this large is unphysical (QC-confirmed unbound)


def _load(dir_, key):
    out = {}
    d = BASE / dir_
    if d.exists():
        for f in d.glob("*.json"):
            try:
                j = json.loads(f.read_text())
                out[j["id"]] = j.get(key)
            except Exception:
                pass
    return out


def _f(x):
    try:
        v = float(x); return v if v == v else None
    except (TypeError, ValueError):
        return None


def main():
    rows = list(csv.DictReader(CSV.open()))
    homo = _load("screen", "anion_homo_at_neu_eV")
    unbound = _load("screen", "anion_unbound")
    cfsmd = _load("recompute", "lambda_cf_smd_eV")
    n_unbound = n_excl = 0
    for r in rows:
        gid = r["id"]; our = _f(r.get("our_smd_eV"))
        r["anion_homo_at_neu_eV"] = round(homo[gid], 3) if homo.get(gid) is not None else ""
        r["anion_unbound"] = unbound.get(gid, "")
        r["lambda_cf_smd_eV"] = cfsmd.get(gid, "")
        unb = bool(unbound.get(gid)) if gid in unbound else False
        rel = not (unb or (our is not None and our > UNRELIABLE_ABOVE))
        r["reliable"] = rel
        n_unbound += unb
        n_excl += (not rel)
    cols = list(rows[0].keys())
    with CSV.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)

    def stats(xk, yk):
        ds = [abs(_f(r[xk]) - _f(r[yk])) for r in rows
              if r["reliable"] in (True, "True") and _f(r.get(xk)) is not None and _f(r.get(yk)) is not None]
        return len(ds), (statistics.mean(ds) if ds else float("nan"))
    print(f"total={len(rows)}  unbound_anion={n_unbound}  excluded(unreliable)={n_excl}  "
          f"reliable={len(rows)-n_excl}")
    for xk, yk, lbl in [("our_smd_eV", "d3tales_eV", "our SMD vs D3TaLES"),
                        ("our_smd_eV", "our_d3level_eV", "our SMD vs our-level")]:
        n, mad = stats(xk, yk)
        print(f"  reliable-only {lbl}: n={n}  MAD={mad:.3f} eV")
    print(f"wrote QC columns to {CSV}")


if __name__ == "__main__":
    main()
