#!/usr/bin/env python
"""Replay the large fixed-density solvent-response diagnostic from raw matrices."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import time


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify(mapping):
    for p, expected in mapping.items():
        if sha(p) != expected:
            raise ValueError("Response evidence changed: " + p)


def review(folder):
    import numpy as np
    plan_path, status_path = folder / "plan.json", folder / "status.json"
    plan, status = (json.loads(p.read_text()) for p in (plan_path, status_path))
    if status["status"] != "response_resource_and_sampled_derivatives_passed":
        raise ValueError("The response pilot has not passed")
    verify(plan["producer_sha256"])
    verify(plan["input_sha256"])
    verify(status["output_sha256"])
    raw = Path(plan["raw_root"])
    baseline = json.loads((raw / "baseline.json").read_text())
    if baseline["SCF_converged"] is not True:
        raise ValueError("The response baseline is not converged")
    expected_coords = np.array(baseline["coordinates_Bohr"])
    analytic = np.load(raw / "analytic_response_MO_occ.npy", mmap_mode="r")
    if analytic.ndim != 4 or analytic.shape[:2] != (len(expected_coords), 3) or not np.isfinite(analytic).all():
        raise ValueError("Analytic response dimensions/values are invalid")
    checks = json.loads((raw / "response_checks.json").read_text())
    expected = {(row["atom"], row["axis"], h) for row in checks for h in plan["settings"]["steps_Bohr"]}
    observed = {(row["atom"], row["axis"], row["step_Bohr"]) for row in checks}
    if len(checks) != 4 or len(expected) != 4 or observed != expected:
        raise ValueError("Two response components at both steps are required")
    output = []
    for row in checks:
        atom, axis, step = row["atom"], row["axis"], row["step_Bohr"]
        pair = []
        for sign in (-1, 1):
            label = f"atom{atom}_axis{axis}_h{step}_sign{sign}"
            info = json.loads((raw / (label+".json")).read_text())
            coords = expected_coords.copy()
            coords[atom, axis] += sign*step
            if (not info["fixed_AO_density"] or not info["fixed_MO_coefficients"] or info["SCF_performed"]
                    or not np.allclose(coords, info["coordinates_Bohr"], rtol=0, atol=1e-11)):
                raise ValueError("Displaced response context/geometry differs")
            value = np.load(raw / (label+".npy"))
            if value.shape != analytic.shape[2:] or not np.isfinite(value).all():
                raise ValueError("Displaced potential matrix has invalid dimensions/values")
            pair.append(value)
        numeric = (pair[1]-pair[0])/(2*step)
        reference = np.asarray(analytic[atom, axis])
        relative = float(np.linalg.norm(numeric-reference)/max(1e-14,np.linalg.norm(reference)))
        if abs(relative-row["relative_error"]) > 1e-12:
            raise ValueError("Reported response error does not reconstruct")
        output.append(dict(atom=atom, axis=axis, step_Bohr=step, relative_error=relative))
    limit = plan["acceptance"]["maximum_response_relative_error"]
    if limit != .01 or max(r["relative_error"] for r in output) > limit:
        raise ValueError("The original one-percent response criterion is not met")
    if (not 0 < status["observed_peak_GPU_MiB"] <= plan["acceptance"]["maximum_GPU_MiB"]
            or status["native_Hessian_qualified"] or status["G_computed"] or status["candidate_admission"]):
        raise ValueError("Resource bound or diagnostic scope is invalid")
    return dict(status="fixed_density_response_matrices_independently_reconstructed", checked=time.time(),
        atoms=len(expected_coords), response_shape=list(analytic.shape), checks=output,
        maximum_relative_error=max(r["relative_error"] for r in output),
        observed_peak_GPU_MiB=status["observed_peak_GPU_MiB"],
        source_sha256={str(p): sha(p) for p in (plan_path, status_path, Path(__file__).resolve())},
        qualified_native_Hessian=False, qualified_G=False, candidate_admission=False,
        scope="Fixed-density PCM response/resource qualification only; no full-Hessian minimum or free energy")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("folder", type=Path)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0,cpu_allocation(1))
    result = review(args.folder.resolve())
    with args.out.open("x") as output:
        output.write(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result))


if __name__ == "__main__":
    main()
