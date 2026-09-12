#!/usr/bin/env python
"""Sample viologen pimer stacking configurations to pin down the config-dependence of the
dimerization dG (addresses the single-configuration error; gives a real sigma).

Grid over (interplanar distance, longitudinal offset, parallel/antiparallel). UMA-relax every
config (cheap), rank by UMA energy, and write the lowest-K relaxed geometries for DFT+SMD. The
spread of the low-lying UMA energies is itself the config-sensitivity estimate.

  CUDA_VISIBLE_DEVICES=0 PYTHONPATH=src python scripts/sample_dimer.py \
      --monomer calcs/dft/methyl_viologen/ox1/opt.xyz --topk 3
"""
from __future__ import annotations
import argparse
from pathlib import Path

from ase.io import read
import redox.uma as U
import redox.dimerize as D

ROOT = Path(__file__).resolve().parents[1]

GRID = [(ip, off, ap)
        for ip in (3.2, 3.5)
        for off in (0.0, 1.6, 3.2)
        for ap in (True, False)]      # 12 configs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--monomer", required=True)
    ap.add_argument("--outdir", default="calcs/uma/mv_dimer")
    ap.add_argument("--topk", type=int, default=3)
    a = ap.parse_args()

    mono = read(a.monomer)
    calc = U.make_calculator("uma-s-1p2p1", "cuda")
    outdir = ROOT / a.outdir
    results = []
    for k, (ip, off, apar) in enumerate(GRID):
        seed = D.build_pi_dimer(mono, interplanar=ip, offset=off, antiparallel=apar)
        # [MV2]2+ pimer: 2 doublets -> closed-shell singlet, total charge +2
        e, fmax, ns, atoms = U.relax_one(calc, seed, 2, 1, fmax=0.05, steps=400)
        d = outdir / f"cfg{k:02d}"; d.mkdir(parents=True, exist_ok=True)
        U.write_xyz(atoms, d / "relaxed.xyz",
                    f"cfg{k} ip={ip} off={off} {'anti' if apar else 'par'} E={e:.4f} fmax={fmax:.3f}")
        results.append((e, k, ip, off, apar, fmax))
        print(f"[cfg{k:02d}] ip={ip} off={off} {'anti' if apar else 'par':4s} "
              f"E_uma={e:.4f} fmax={fmax:.3f}", flush=True)

    results.sort()
    e0 = results[0][0]
    print("\n[rank] lowest-energy configs (UMA, relative meV):")
    for e, k, ip, off, apar, fmax in results:
        print(f"  cfg{k:02d}  dE={ (e-e0)*1000:8.1f} meV  ip={ip} off={off} {'anti' if apar else 'par'}")
    spread = (results[min(a.topk, len(results))-1][0] - e0) * 1000
    print(f"\n[spread] top-{a.topk} UMA energy spread = {spread:.1f} meV (config sensitivity)")
    topk = [k for _, k, *_ in results[:a.topk]]
    (outdir / "topk.txt").write_text(",".join(f"cfg{k:02d}" for k in topk))
    print(f"[topk] {topk} -> {outdir}/topk.txt  (DFT these next)")


if __name__ == "__main__":
    main()
