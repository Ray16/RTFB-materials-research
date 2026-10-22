#!/usr/bin/env python
"""Isolate finite-difference CDS asymmetry in a saved native Hessian.

Uses only the CPU legacy SMD CDS primitive, replayed against saved GPU control
components. Refined matrices are diagnostic until native and worst-column
derivative/minimum checks pass. No raw state or acceptance threshold changes.
"""

import argparse
from importlib import metadata
import json
import os
from pathlib import Path
import time

import numpy as np

from redox.core.io import atomic_write_json
from redox.workflows.chunked import file_hash


def finite_gradient_hessian(coordinates, gradient, step):
    if not np.isfinite(step) or step <= 0:
        raise ValueError("Positive finite displacement required")
    coordinates = np.asarray(coordinates, dtype=float)
    n = len(coordinates)
    result = np.empty((n, n, 3, 3))
    # Reproduce the installed library's +h, -2h coordinate arithmetic.
    working = coordinates.copy()
    for atom in range(n):
        for axis in range(3):
            working[atom, axis] += step
            plus = np.asarray(gradient(working.copy())).copy()
            working[atom, axis] -= 2 * step
            minus = np.asarray(gradient(working.copy())).copy()
            result[atom, :, axis] = (plus - minus) / (2 * step)
            working[atom, axis] = coordinates[atom, axis]
    return result.transpose(0, 2, 1, 3).reshape(3 * n, 3 * n)


def cds_for_state(state, step):
    from ase.io import read
    from pyscf import gto
    from pyscf.solvent import smd

    state = Path(state)
    raw = json.loads((state / "result.json").read_text())
    atoms = read(state / "opt.xyz")
    coordinates = atoms.positions / 0.52917721067
    mol = gto.M(
        atom=list(zip(atoms.get_chemical_symbols(), coordinates.tolist())),
        unit="Bohr",
        basis="def2-svp",
        charge=raw["charge"],
        spin=raw["mult"] - 1,
        verbose=0,
    )
    solvent = smd.SMD(mol, solvent="acetonitrile")

    def gradient(position):
        mol.set_geom_(position, unit="Bohr")
        solvent.reset(mol)
        return smd.get_cds_legacy(solvent)[1]

    return finite_gradient_hessian(coordinates, gradient, step)


def matrix_summary(matrix):
    if not np.isfinite(matrix).all():
        raise ValueError("Nonfinite Hessian component")
    skew = matrix - matrix.T
    i, j = np.unravel_index(abs(skew).argmax(), skew.shape)
    return dict(
        maximum_asymmetry=float(abs(skew).max()),
        relative_asymmetry=float(
            np.linalg.norm(skew) / max(np.linalg.norm(matrix), 1e-30)
        ),
        worst_indices=[int(i), int(j)],
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--hessian", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--control", type=Path, nargs=2, action="append", default=[])
    args = parser.parse_args()
    from redox.workflows.resources import (
        allocation_guard,
        cpu_allocation,
        host_capacity,
    )

    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    args.out.mkdir(parents=True, exist_ok=False)
    raw = np.load(args.hessian)
    report = dict(
        status="running_CPU_component_diagnostic",
        started=time.time(),
        protocol=dict(solvent="acetonitrile", CDS_steps_Bohr=[1e-4, 1e-5, 1e-6]),
        source_hessian_sha256=file_hash(args.hessian),
        source_geometry_sha256=file_hash(args.state / "opt.xyz"),
        producer_sha256=file_hash(__file__),
        software=dict(pyscf=metadata.version("pyscf")),
        original=matrix_summary(raw),
        controls=[],
        components={},
        replacements={},
        physical_QC_promoted=False,
    )
    atomic_write_json(args.out / "status.json", report)
    try:
        if len(args.control) != 2:
            raise ValueError(
                "Both saved singlet/radical component replays are required"
            )
        for index, (state, baseline) in enumerate(args.control):
            reproduced = cds_for_state(state, 1e-4)
            np.save(args.out / f"control_{index}.npy", reproduced)
            prior = np.load(baseline)
            maximum = float(abs(reproduced - prior).max())
            report["controls"].append(
                dict(
                    state=str(state),
                    baseline=str(baseline),
                    baseline_sha256=file_hash(baseline),
                    geometry_sha256=file_hash(state / "opt.xyz"),
                    maximum_difference=maximum,
                )
            )
            if maximum > 1e-9:
                raise ValueError("CPU CDS primitive differs from saved GPU baseline")
        matrices = {}
        for step in (1e-4, 1e-5, 1e-6):
            matrix = cds_for_state(args.state, step)
            matrices[step] = matrix
            name = f"CDS_hessian_{step}.npy"
            np.save(args.out / name, matrix)
            report["components"][str(step)] = dict(
                matrix_summary(matrix), sha256=file_hash(args.out / name)
            )
            if step != 1e-4:
                replaced = raw - matrices[1e-4] + matrix
                filename = f"diagnostic_replacement_{step}.npy"
                np.save(args.out / filename, replaced)
                report["replacements"][str(step)] = dict(
                    matrix_summary(replaced), sha256=file_hash(args.out / filename)
                )
            atomic_write_json(args.out / "status.json", report)
        if file_hash(args.hessian) != report["source_hessian_sha256"]:
            raise ValueError("Source Hessian changed during diagnostic")
        report.update(
            status="component_diagnostic_complete_no_admission", finished=time.time()
        )
        atomic_write_json(args.out / "status.json", report)
        print(json.dumps(report["replacements"], indent=2))
    except BaseException as error:
        report.update(
            status="failed_preserved", error=repr(error), finished=time.time()
        )
        atomic_write_json(args.out / "status.json", report)
        raise


if __name__ == "__main__":
    main()
