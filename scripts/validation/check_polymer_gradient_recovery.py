#!/usr/bin/env python
"""Independently reconstruct sampled gradient checks from saved total energies."""
import argparse
import json
import math
import os
from pathlib import Path
import time

from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def reconstruct(campaign):
    plan = api.read_plan(campaign)
    api.verify_sources(plan["producer_sha256"])
    api.verify_sources(plan["input_sha256"])
    path = campaign / "gradient_qualification.json"
    report = json.loads(path.read_text())
    if (report["status"] != "fixed_geometry_gradient_and_sampled_derivatives_passed"
            or report["free_energy_computed"] or report["native_Hessian_qualified"]
            or report["full_large_geometry_workflow_qualified"] or report["synthesis_admission"]):
        raise ValueError("Fixed-point evidence has invalid status or was promoted beyond its scope")
    api.verify_sources(report["output_sha256"])
    folder = Path(plan["raw_root"]) / "energy_force_check"
    base = json.loads((folder / "baseline.json").read_text())
    force = base["gradient_Ha_Bohr"]
    shape = len(force)
    if not base["SCF_converged"] or any(len(row) != 3 or any(not math.isfinite(x) for x in row) for row in force):
        raise ValueError("Baseline SCF/gradient is not finite and converged")
    selected = sorted(((abs(force[i][j]), i, j) for i in range(shape) for j in range(3)))[-2:]
    checks, sources = [], {str(path): api.digest(path), str(folder / "baseline.json"): api.digest(folder / "baseline.json")}
    for _, atom, axis in selected:
        for step in (.0005, .001):
            energies = []
            for sign in (-1, 1):
                p = folder / f"atom{atom}_axis{axis}_h{step}_sign{sign}.json"
                data = json.loads(p.read_text())
                if not data["SCF_converged"] or not math.isfinite(data["E_Ha"]):
                    raise ValueError("Displaced SCF is incomplete or nonfinite")
                for i in range(shape):
                    for j in range(3):
                        expected = base["coordinates_Bohr"][i][j] + (sign * step if (i, j) == (atom, axis) else 0)
                        if abs(data["coordinates_Bohr"][i][j] - expected) > 1e-12:
                            raise ValueError("Displaced coordinates differ from the declared component/step")
                sources[str(p)] = api.digest(p)
                energies.append(data["E_Ha"])
            numerical = (energies[1] - energies[0]) / (2 * step)
            checks.append(dict(atom=atom, axis=axis, step_Bohr=step,
                analytic=force[atom][axis], numerical=numerical,
                absolute_error=abs(numerical - force[atom][axis])))
    if checks != report["checks"]["checks"] or report["checks"]["acceptance_Ha_Bohr"] != 1e-6:
        raise ValueError("Independent derivatives differ from the recorded gate")
    if any(c["absolute_error"] > 1e-6 for c in checks):
        raise ValueError("Independent derivative check exceeds the unchanged tolerance")
    maximum = max(abs(x) for row in force for x in row)
    previous = json.loads(Path(plan["previous_reference"]).read_text())
    difference = abs(maximum - previous["gradient_max_Ha_Bohr"])
    peak = report["observed_peak_GPU_MiB"]
    if (maximum != report["gradient_max_Ha_Bohr"] or difference != report["previous_gradient_max_difference"]
            or not 0 <= difference <= 1e-6 or not 0 < peak <= plan["maximum_memory_pilot_MiB"]):
        raise ValueError("Reference-gradient or resource-continuity bound failed")
    return dict(status="independent_fixed_geometry_derivatives_reconstructed", checked=time.time(),
        atoms=shape, sampled_components=2, displaced_SCF_points=8, checks=checks,
        maximum_error_Ha_Bohr=max(c["absolute_error"] for c in checks),
        gradient_max_reference_difference_Ha_Bohr=difference, observed_peak_GPU_MiB=peak,
        source_sha256=dict(sources, **{str(Path(__file__).resolve()): api.digest(__file__),
            str(campaign / "plan.json"): api.digest(campaign / "plan.json")}),
        scope="Fixed geometry, two components, two steps; no minimum, Hessian, ensemble G or candidate admission")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    if args.out.exists():
        raise FileExistsError("Preserve prior independent review")
    review = reconstruct(args.campaign.resolve())
    atomic_write_json(args.out, review)
    print(json.dumps({k: review[k] for k in ("status", "atoms", "maximum_error_Ha_Bohr", "observed_peak_GPU_MiB")}))


if __name__ == "__main__":
    main()
