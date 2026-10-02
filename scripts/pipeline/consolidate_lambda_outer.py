#!/usr/bin/env python
"""Merge per-chunk outer-sphere lambda runs (scripts/pipeline/compute_lambda_outer.py) into the canonical
results/lambda_outer.csv read by redox.screening.scorecard. One row per (molecule, couple).

  python scripts/pipeline/consolidate_lambda_outer.py
"""
import csv
import glob
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
# v2 (diffuse basis, SCF-checked, both geometries, metadata) is read FIRST so it wins over a
# legacy row for the same (molecule, couple); legacy rows remain only where v2 has no value.
SRC = (sorted(glob.glob(str(ROOT / "results/reorg_anchors/our_chunks_v2/chunk*_results.csv")))
       + sorted(glob.glob(str(ROOT / "results/reorg_anchors/our_chunks/chunk*_results.csv"))))
OUT = ROOT / "results/lambda_outer.csv"


def main():
    rows, seen = [], set()
    for f in SRC:
        for r in csv.DictReader(open(f)):
            couple = r.get("couple", "")
            suffix = "_" + couple.replace("/", "_")
            mol = r["id"][: -len(suffix)] if r["id"].endswith(suffix) else r["id"]
            key = (mol, couple)
            if key in seen:
                continue
            seen.add(key)
            rows.append(dict(id=mol, couple=couple,
                             lambda_o_born_eV=r["born_lambda_o_eV"],
                             lambda_o_pcm_eV=r["pcm_lambda_o_eV"],
                             pcm_over_born=r["pcm_over_born"],
                             sasa_radius_A=r["sasa_radius_A"], solvent=r["solvent"],
                             # provenance: chunks written before the metadata columns existed
                             # used the old defaults (B3LYP/6-31G(d,p), no diffuse, SCF
                             # convergence NOT checked) — labelled so, never silently merged.
                             method=(f"{r['xc']}/{r['basis']}" if r.get("xc")
                                     else "b3lyp/6-31g(d,p) (legacy)"),
                             scf_checked=bool(r.get("scf_converged")),
                             geom_source=r.get("geom_source", "")))
    cols = ["id", "couple", "lambda_o_born_eV", "lambda_o_pcm_eV",
            "pcm_over_born", "sasa_radius_A", "solvent", "method", "scf_checked", "geom_source"]
    with OUT.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols); w.writeheader(); w.writerows(rows)
    print(f"wrote {OUT} ({len(rows)} couples from {len(SRC)} chunk files)")


if __name__ == "__main__":
    main()
