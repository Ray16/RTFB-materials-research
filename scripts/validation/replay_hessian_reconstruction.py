#!/usr/bin/env python
"""Replay additive CDS replacement against two full independently recomputed Hessians."""

import argparse
import json
import os
from pathlib import Path

import numpy as np

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.qm import native_thermal as native, qrrho
from redox.qm.hessian_memory import derivative_certificate
from redox.workflows import polymer_campaign as api


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    from redox.workflows.resources import (
        allocation_guard,
        cpu_allocation,
        host_capacity,
    )

    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    certificate = derivative_certificate(native.NATIVE_THERMAL, 1e-5)
    args.out.mkdir(parents=True, exist_ok=False)
    helper = api.load_module(
        "CDS_component_replay", ROOT / "scripts/validation/diagnose_saved_hessian.py"
    )
    rows = []
    from ase.io import read

    for index, (state, mult) in enumerate([("neu", 1), ("red1", 2)]):
        source = ROOT / f"calcs/dft/exp_naphthalimide_control/{state}"
        address = next((source / "native_thermal").glob("*.json")).stem
        baseline = source / "native_thermal" / address / "hessian_cartesian.npy"
        full = (
            ROOT
            / f"results/validation/resume_campaign_2026-10-12/physical_campaign/cds_qualification/exp_naphthalimide_control__{state}/work/native_thermal/{address}/hessian_cartesian.npy"
        )
        coarse = (
            ROOT
            / f"results/validation/polymer_hessian_components_2026-10-21/control_{index}.npy"
        )
        fine = helper.cds_for_state(source, 1e-5)
        reconstructed = np.load(baseline) - np.load(coarse) + fine
        actual = np.load(full)
        maximum = float(abs(actual - reconstructed).max())
        relative = float(
            np.linalg.norm(actual - reconstructed) / np.linalg.norm(actual)
        )
        atoms = read(source / "opt.xyz")
        nuclear1 = qrrho.nuclear_thermals(
            reconstructed, atoms.positions, atoms.numbers, 1
        )
        nuclear2 = qrrho.nuclear_thermals(actual, atoms.positions, atoms.numbers, 1)
        thermal = max(
            abs(
                nuclear1["nuclear_by_cutoff_eV"][tag]
                - nuclear2["nuclear_by_cutoff_eV"][tag]
            )
            for tag in ("25", "50", "100")
        )
        np.save(args.out / f"reconstructed_{state}.npy", reconstructed)
        row = dict(
            mult=mult,
            state=state,
            maximum_matrix_difference=maximum,
            relative_matrix_error=relative,
            maximum_thermal_difference_eV=thermal,
            same_symmetry_number_used_for_difference_only=1,
            input_sha256={
                str(p): api.digest(p)
                for p in (baseline, full, coarse, source / "opt.xyz")
            },
        )
        rows.append(row)
        if max(maximum, relative, thermal) > 1e-6:
            raise ValueError(
                "Cached reconstruction differs from full native recomputation"
            )
    report = dict(
        status="two_full_recomputation_replays_passed",
        CDS_step_Bohr=1e-5,
        qualification_sha256=certificate,
        controls=rows,
        producer_sha256=api.digest(__file__),
        component_helper_sha256=api.digest(helper.__file__),
        scientific_model_changed=False,
        thresholds_changed=False,
    )
    atomic_write_json(args.out / "summary.json", report)
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
