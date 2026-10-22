#!/usr/bin/env python
"""Check the legacy and conservative CDS forces without any electronic SCF."""

import argparse
from concurrent.futures import ProcessPoolExecutor
from importlib import metadata
import json
import os
from pathlib import Path
import time

import numpy as np

from redox.core.io import atomic_write_json
from redox.qm.cds_finite import energy_gradient
from redox.workflows.chunked import file_hash


def evaluate(spec):
    from redox.workflows.resources import (
        allocation_guard,
        cpu_allocation,
        host_capacity,
    )

    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    from ase.io import read
    from pyscf import gto
    from pyscf.solvent import smd

    atoms = read(spec["xyz"])
    coordinates = atoms.positions / 0.52917721067
    mol = gto.M(
        atom=list(zip(atoms.get_chemical_symbols(), coordinates.tolist())),
        unit="Bohr",
        basis="def2-svp",
        charge=spec["charge"],
        spin=spec["mult"] - 1,
        verbose=0,
    )
    solvent = smd.SMD(mol, solvent="acetonitrile")
    _, legacy = smd.get_cds_legacy(solvent)

    def energy(position):
        mol.set_geom_(position, unit="Bohr")
        solvent.reset(mol)
        return smd.get_cds_legacy(solvent)[0]

    numeric = energy_gradient(coordinates, energy, 1e-4)
    finer = energy_gradient(coordinates, energy, 5e-5)
    worst = np.unravel_index(abs(numeric - legacy).argmax(), legacy.shape)
    result = dict(
        spec,
        geometry_sha256=file_hash(spec["xyz"]),
        steps_Bohr=[1e-4, 5e-5],
        gradient_legacy=legacy.tolist(),
        gradient_energy_fd=numeric.tolist(),
        gradient_energy_fd_finer=finer.tolist(),
        maximum_legacy_disagreement=float(abs(numeric - legacy).max()),
        maximum_step_disagreement=float(abs(numeric - finer).max()),
        worst_component=[int(x) for x in worst],
        CDS_energy_function_unchanged=True,
    )
    comparisons = []
    for name in spec.get("full_gradient_diagnostics", []):
        root = Path(name)
        status = json.loads((root / "status.json").read_text())
        baseline = json.loads((root / "baseline.json").read_text())
        if status["source_geometry_sha256"] != result["geometry_sha256"]:
            raise ValueError("Full-SCF diagnostic uses a different geometry")
        corrected = np.asarray(baseline["gradient_Ha_Bohr"]) - legacy + numeric
        for check in status["differences"]:
            atom, axis = check["atom"], check["axis"]
            comparisons.append(
                dict(
                    grid=status["settings"]["grid_level"],
                    atom=atom,
                    axis=axis,
                    step_Bohr=check["step_Bohr"],
                    corrected=float(corrected[atom, axis]),
                    numerical=check["numerical"],
                    absolute_error=abs(corrected[atom, axis] - check["numerical"]),
                )
            )
    result["full_energy_difference_checks"] = comparisons
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    cases = json.loads(args.manifest.read_text())
    args.out.mkdir(parents=True, exist_ok=False)
    with ProcessPoolExecutor(max_workers=len(cases)) as pool:
        results = list(pool.map(evaluate, cases))
    for result in results:
        atomic_write_json(args.out / (result["id"] + ".json"), result)
    valid = all(r["maximum_step_disagreement"] <= 1e-6 for r in results) and all(
        c["absolute_error"] <= 1e-6
        for r in results
        for c in r["full_energy_difference_checks"]
    )
    summary = dict(
        status="energy_gradient_diagnostics_passed"
        if valid
        else "diagnostics_require_review",
        cases=results,
        producer_sha256=file_hash(__file__),
        derivative_module_sha256=file_hash(
            Path(__file__).resolve().parents[2] / "src/redox/qm/cds_finite.py"
        ),
        pyscf_version=metadata.version("pyscf"),
        checked_at=time.time(),
        acceptance_Ha_Bohr=1e-6,
        optimizer_production_admission=False,
    )
    atomic_write_json(args.out / "summary.json", summary)
    print(
        json.dumps(
            dict(
                status=summary["status"],
                cases=len(results),
                max_step_difference=max(
                    r["maximum_step_disagreement"] for r in results
                ),
            )
        )
    )
    if not valid:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
