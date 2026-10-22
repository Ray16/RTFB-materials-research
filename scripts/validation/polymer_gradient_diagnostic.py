#!/usr/bin/env python
"""Compare conservative SMD gradients with energy differences at a saved geometry.

Grid refinements are diagnostics; no geometry or production record is promoted.
Each SCF is checkpointed on shared storage before the next displacement.
"""

import argparse
import os
from pathlib import Path
import socket
import time

import numpy as np

from redox.core.io import atomic_write_json
from redox.workflows.chunked import file_hash


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xyz", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--grid-level", type=int, choices=(3, 5), required=True)
    args = parser.parse_args()
    from redox.workflows.resources import verify_gate_owner

    verify_gate_owner()
    if len(os.sched_getaffinity(0)) != 2:
        raise RuntimeError("Two disjoint CPU slots required")
    from ase.io import read
    from pyscf import gto
    from redox.qm.dft import _build_smd_mf, _free_gpu_pool
    from redox.qm.stable_geometry_scf import stabilize

    args.out.mkdir(parents=True, exist_ok=False)
    atoms = read(args.xyz)
    symbols = atoms.get_chemical_symbols()
    coordinates = atoms.positions / 0.52917721067
    settings = dict(
        xc="r2scan",
        disp="d4",
        basis="def2-svp",
        solvent="acetonitrile",
        solvent_model="SMD",
        grid_level=args.grid_level,
        grid_response=True,
        conv_tol=1e-12,
        conv_tol_grad=1e-8,
        charge=0,
        mult=1,
        steps_Bohr=[0.0005, 0.001],
    )
    state = dict(
        status="computing",
        host=socket.gethostname(),
        pid=os.getpid(),
        process_start=Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()[19],
        cpus=sorted(os.sched_getaffinity(0)),
        gpu=os.environ.get("CUDA_VISIBLE_DEVICES"),
        source_xyz=str(args.xyz.resolve()),
        source_geometry_sha256=file_hash(args.xyz),
        settings=settings,
        producer_sha256=file_hash(__file__),
        evaluations=[],
        started=time.time(),
    )
    atomic_write_json(args.out / "status.json", state)
    baseline = None

    def solve(position, name, gradient=False):
        mol = gto.M(
            atom=list(zip(symbols, position.tolist())),
            unit="Bohr",
            charge=0,
            spin=0,
            basis="def2-svp",
            verbose=0,
            max_memory=6000,
        )
        mf = _build_smd_mf(mol, "r2scan", "acetonitrile", "gpu", disp="d4")
        mf.grids.level = args.grid_level
        mf.conv_tol, mf.conv_tol_grad, mf.max_cycle = 1e-12, 1e-8, 300
        stabilize(mf, args.out / "scf_recovery.jsonl")
        mf.kernel(dm0=baseline)
        if not mf.converged:
            raise RuntimeError("Diagnostic SCF did not converge")
        record = dict(
            name=name,
            E_Ha=float(mf.e_tot),
            SCF_converged=True,
            coordinates_Bohr=position.tolist(),
        )
        raw_density = mf.make_rdm1()
        density = np.asarray(
            raw_density.get() if hasattr(raw_density, "get") else raw_density
        )
        grad = None
        if gradient:
            method = mf.nuc_grad_method()
            method.grid_response = True
            raw_gradient = method.kernel()
            grad = np.asarray(
                raw_gradient.get() if hasattr(raw_gradient, "get") else raw_gradient
            )
            record.update(
                gradient_Ha_Bohr=grad.tolist(),
                maximum_gradient_component=float(abs(grad).max()),
            )
        atomic_write_json(args.out / (name + ".json"), record)
        state["evaluations"].append(
            dict(name=name, file_sha256=file_hash(args.out / (name + ".json")))
        )
        atomic_write_json(args.out / "status.json", state)
        del mf, mol
        _free_gpu_pool()
        return record["E_Ha"], density, grad

    try:
        _, baseline, gradient = solve(coordinates, "baseline", True)
        worst = tuple(
            int(x) for x in np.unravel_index(abs(gradient).argmax(), gradient.shape)
        )
        axes = list(dict.fromkeys([worst, (0, 1)]))
        differences = []
        for atom, axis in axes:
            for step in settings["steps_Bohr"]:
                energies = []
                for sign in (-1, 1):
                    shifted = coordinates.copy()
                    shifted[atom, axis] += sign * step
                    energy, _, _ = solve(
                        shifted, f"atom{atom}_axis{axis}_h{step}_sign{sign}"
                    )
                    energies.append(energy)
                numerical = (energies[1] - energies[0]) / (2 * step)
                differences.append(
                    dict(
                        atom=atom,
                        axis=axis,
                        step_Bohr=step,
                        analytic=float(gradient[atom, axis]),
                        numerical=numerical,
                        absolute_difference=abs(numerical - gradient[atom, axis]),
                    )
                )
        state.update(
            status="diagnostic_complete_no_promotion",
            worst_component=worst,
            differences=differences,
            finished=time.time(),
        )
        atomic_write_json(args.out / "status.json", state)
    except BaseException as error:
        state.update(status="failed_preserved", error=repr(error), finished=time.time())
        atomic_write_json(args.out / "status.json", state)
        raise


if __name__ == "__main__":
    main()
