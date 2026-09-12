#!/usr/bin/env python
"""Resumable worker: recompute D3TaLES inner-sphere reorganization energies with OUR EXACT
PRODUCTION pipeline, so the population is directly comparable to the candidates' lambda_i.

Drives redox.dft.dft_smd (the same function the pipeline uses), reproducing reorg.py:
  - geometries OPTIMIZED IN SMD(MeCN) at r2SCAN/def2-SVP (def2-SVPD for the anion),
  - 4-point energies = GAS single points at wB97M-V/def2-TZVP (def2-TZVPD for the anion),
    evaluated on those SMD-optimized geometries (do_smd=False -> the gas e_gas_eV that
    reorg.py's lambda_i is built from).
  - ELECTRON couple (neutral singlet -> radical-anion doublet), matching the reductive
    (anolyte) D3TaLES quinone/imide population and our reductive candidates.

Writes one JSON per molecule (+ opt geometries) under results/d3tales_reorg_validation_prod/
-> a NEW tree; it NEVER touches the B3LYP results in results/d3tales_reorg_validation/.
Resumable: a molecule with a finished JSON is skipped, so fleet shards/reruns never redo work.

  CUDA_VISIBLE_DEVICES=<idx> PYTHONPATH=src python scripts/validate_reorg_worker_prod.py \
      --source results/d3tales_reorg_validation/molecules.csv --offset 0 --stride 24 --backend gpu
"""
from __future__ import annotations
import argparse, json, traceback
from pathlib import Path
import numpy as np
import pandas as pd
from ase.io import write as ase_write
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

import redox.dft as D
from redox.nelsen import embed   # RDKit ETKDG+MMFF neutral starting geometry

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "results" / "d3tales_reorg_validation_prod"
CALC = BASE / "calc"
GEO = BASE / "geom"


def _gas(res):
    return res.get("e_gas_eV")


def _free_gpu(backend):
    """Release the cupy memory pool between DFT stages. The pipeline runs one process per
    state, so GPU memory is always fresh; this worker chains 4 dft_smd calls per molecule, so
    without freeing, gpu4pyscf/cupy accumulate allocations and OOM even on a 32 GB V100."""
    if backend == "gpu":
        try:
            import cupy
            cupy.get_default_memory_pool().free_all_blocks()
            cupy.get_default_pinned_memory_pool().free_all_blocks()
        except Exception:
            pass


def process(row, backend):
    gid = row["id"]
    out = CALC / f"{gid}.json"
    if out.exists():
        try:
            if json.loads(out.read_text()).get("status") in ("ok", "partial"):
                return "skip"
        except Exception:
            pass
    rec = dict(id=gid, smiles=row["smiles"], family=row.get("family", ""),
               d3_electron=row.get("d3_electron", ""), d3_hole=row.get("d3_hole", ""),
               n_atoms=row.get("n_atoms", ""), n_rot=row.get("n_rot", ""),
               protocol="pipeline-exact: SMD-opt r2SCAN/def2-SVP(D) + gas wB97M-V/def2-TZVP(D)",
               our_electron=None, relax_neu_meV=None, relax_anion_meV=None, status="ok", error="")
    GEO.mkdir(parents=True, exist_ok=True)
    try:
        # neutral starting geometry (RDKit) -> xyz
        neu0 = GEO / f"{gid}_neu_start.xyz"; ase_write(str(neu0), embed(row["smiles"]))
        neu_opt = GEO / f"{gid}_neu_opt.xyz"; an_opt = GEO / f"{gid}_anion_opt.xyz"

        # SMD-opt neutral (0, singlet); gas SP -> E_O_at_O
        rNeu = D.dft_smd(neu0, 0, 1, do_opt=True, do_gas=True, do_smd=False, do_freq=False,
                         opt_out=neu_opt, backend=backend)
        E_O_at_O = _gas(rNeu); _free_gpu(backend)
        # SMD-opt anion (-1, doublet) warm-started from the neutral SMD geom; gas SP -> E_R_at_R
        rAn = D.dft_smd(neu_opt, -1, 2, do_opt=True, do_gas=True, do_smd=False, do_freq=False,
                        opt_out=an_opt, backend=backend)
        E_R_at_R = _gas(rAn); _free_gpu(backend)
        # cross points (single points, no opt): E_O_at_R (neutral charge on anion geom),
        # E_R_at_O (anion charge on neutral geom)
        E_O_at_R = _gas(D.dft_smd(an_opt, 0, 1, do_opt=False, do_gas=True, do_smd=False,
                                  do_freq=False, backend=backend)); _free_gpu(backend)
        E_R_at_O = _gas(D.dft_smd(neu_opt, -1, 2, do_opt=False, do_gas=True, do_smd=False,
                                  do_freq=False, backend=backend)); _free_gpu(backend)
        lam = (E_O_at_R - E_O_at_O) + (E_R_at_O - E_R_at_R)
        rec.update(our_electron=lam,
                   relax_neu_meV=round((E_R_at_O - E_R_at_R) * 1000, 1),   # anion relaxing to its geom
                   relax_anion_meV=round((E_O_at_R - E_O_at_O) * 1000, 1),
                   converged=all(bool(r.get("converged_gas", True)) for r in (rNeu, rAn)))
    except Exception as e:
        rec["status"] = "error"; rec["error"] = f"{type(e).__name__}: {e}"
        rec["trace"] = traceback.format_exc()[-800:]
    CALC.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(rec, indent=1))
    return rec["status"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", default=str(ROOT / "results" / "d3tales_reorg_validation" / "molecules.csv"))
    ap.add_argument("--offset", type=int, default=0)
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--backend", default="gpu")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--only", default=None, help="process a single molecule id (for claim-based fleet)")
    a = ap.parse_args()
    df = pd.read_csv(a.source)
    if a.only:
        df = df[df["id"] == a.only].reset_index(drop=True)
        if df.empty:
            print(f"[only] id {a.only!r} not in {a.source}"); return
    else:
        df = df.iloc[a.offset::a.stride].reset_index(drop=True)
    if a.limit:
        df = df.head(a.limit)
    n_ok = n_skip = n_err = 0
    for _, row in df.iterrows():
        st = process(row, a.backend)
        n_ok += st == "ok"; n_skip += st == "skip"; n_err += st == "error"
        tag = ""
        if st == "ok":
            tag = f"  lam_e={json.loads((CALC/(row['id']+'.json')).read_text()).get('our_electron')}"
        print(f"[{row['id']}] {st}{tag}", flush=True)
        if st in ("ok", "skip"):
            print("DONE_MARKER", flush=True)   # fleet worker success signal
    print(f"shard offset={a.offset} stride={a.stride}: ok={n_ok} skip={n_skip} err={n_err}", flush=True)


if __name__ == "__main__":
    main()
