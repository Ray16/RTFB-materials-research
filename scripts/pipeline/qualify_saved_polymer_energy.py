#!/usr/bin/env python
"""Create a separate energy qualification after an explained coarse-step failure.

The failed source is immutable. A new force-validation policy requires both
component attribution and independently reconstructed fine-step quantum checks,
with the original 1e-6 Ha/Bohr tolerance. Native Hessian/G gates are unchanged.
"""
import argparse
import json
import math
import os
from pathlib import Path
import shutil
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

POLICY = dict(id="CDS_explained_step_converged_total_force_v1",
    absolute_force_error_Ha_Bohr=1e-6, CDS_attribution_residual_Ha_Bohr=1e-6,
    required_quantum_steps_Bohr=[0.0005, 0.0001, 0.00005, 0.00001],
    original_0_001_Bohr_failure_retained=True, baseline_repeats=3,
    original_native_Hessian_protocol_unchanged=True)


def qualification_checks(coarse, fine, surface):
    if coarse["passed"] or coarse["acceptance_Ha_Bohr"] != 1e-6:
        raise ValueError("This policy requires the preserved failed original check")
    if not fine["diagnostic_criteria_passed"] or not fine["derivative_criteria_passed"]:
        raise ValueError("The predeclared fine-step criteria did not pass")
    if (not all(math.isfinite(surface[k]) for k in ("maximum_remaining_difference_Ha_Bohr", "maximum_library_energy_difference_Ha"))
            or surface["maximum_remaining_difference_Ha_Bohr"] > 1e-6
            or surface["maximum_library_energy_difference_Ha"] > 1e-12):
        raise ValueError("The CDS term does not explain the coarse-step discrepancy")
    original = {(r["atom"], r["axis"]) for r in coarse["checks"]}
    fine_pairs = {(r["atom"], r["axis"]) for r in fine["checks"]}
    if len(original) < 2 or not original <= fine_pairs:
        raise ValueError("Original failing components are missing from the quantum refinement")
    checks = [dict(r) for r in coarse["checks"] if r["step_Bohr"] == 0.0005]
    checks += [dict(atom=r["atom"], axis=r["axis"], step_Bohr=r["step_Bohr"],
        analytic=r["analytic"], numerical=r["numeric"], absolute_error=r["error_Ha_Bohr"])
        for r in fine["checks"]]
    for a, x in original:
        steps = {r["step_Bohr"] for r in checks if (r["atom"], r["axis"]) == (a, x)}
        if steps != set(POLICY["required_quantum_steps_Bohr"]):
            raise ValueError("Required four-step coverage is incomplete")
    if not all(math.isfinite(r["absolute_error"]) and 0 <= r["absolute_error"] <= 1e-6 for r in checks):
        raise ValueError("A quantum force error exceeds the unchanged acceptance bound")
    return dict(passed=True, checks=checks, acceptance_Ha_Bohr=1e-6,
        validation_policy=POLICY, original_coarse_check_passed=False)


def prepare(args):
    source, fine_root, surface_root = args.source.resolve(), args.fine.resolve(), args.surface.resolve()
    old = api.read_plan(source)
    api.verify_sources(old["producer_sha256"])
    api.verify_sources(old["input_sha256"])
    outer_path = source / "workspace_energy_status.json"
    outer = json.loads(outer_path.read_text())
    if outer["status"] != "energy_attempt_failed_preserved" or "total force" not in outer["error"]:
        raise ValueError("The source is not the observed final-force failure")
    if len(old["tasks"]) != 1 or (old["tasks"][0]["charge"], old["tasks"][0]["mult"]) != (0, 1):
        raise ValueError("One existing neutral-singlet polymer is in scope")
    task = old["tasks"][0]
    original_path = source / "energy/tasks" / task["key"] / "status.json"
    report = json.loads(original_path.read_text())
    if report["status"] != "state_records_verified" or report["physical_QC"] != "passed":
        raise ValueError("Original geometry/active-energy audit is incomplete")
    api.verify_sources(report["output_sha256"])
    state = Path(report["state_dir"])
    geometry = state / "opt.xyz"
    coarse_path = state / "total_energy_force_check/summary.json"
    coarse = json.loads(coarse_path.read_text())
    fine_plan = api.read_plan(fine_root)
    if (fine_plan["source_campaign"] != str(source)
            or api.digest(fine_root / "input.xyz") != api.digest(geometry)
            or coarse["geometry_sha256"] != api.digest(geometry)):
        raise ValueError("Fine-step review does not belong to the exact source geometry")
    fine_script = ROOT / "scripts/validation/review_polymer_gradient_profile.py"
    surface_script = ROOT / "scripts/validation/review_cds_surface.py"
    fine = api.load_module("saved_energy_independent_fine_review", fine_script).review(fine_root)
    surface = api.load_module("saved_energy_independent_surface_review", surface_script).review(surface_root / "checked")
    surface_data = json.loads((surface_root / "checked/summary.json").read_text())
    if surface_data["source_sha256"].get(str(geometry)) != api.digest(geometry):
        raise ValueError("CDS attribution belongs to another geometry")
    checks = qualification_checks(coarse, fine, surface)
    baseline = json.loads((Path(fine_plan["raw_root"]) / "diagnostic/baseline.json").read_text())
    maximum = max(abs(v) for xyz in baseline["gradient_Ha_Bohr"] for v in xyz)
    if not math.isfinite(maximum) or maximum > 0.00045:
        raise ValueError("The independently checked geometry is not stationary")
    from redox.core.protocol import load_record
    active = load_record(state, 0, 1)
    if active is None or abs(active["e_smd_eV"]-report["E_smd_eV"]) > 1e-8:
        raise ValueError("Stored active energy no longer reconstructs")
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    shutil.copytree(source / "live_inputs", out / "live_inputs")
    shutil.copytree(source / "controls", out / "controls")
    for stage in ("energy", "thermal"):
        for d in ("tasks", "done", "claims", "failures", "workers", "pilot", "launches"):
            (out / stage / d).mkdir(parents=True)
    inputs = dict(old["input_sha256"])
    for p in (source / "plan.json", outer_path, original_path, geometry, coarse_path,
              fine_root / "plan.json", fine_root / "status.json", surface_root / "checked/summary.json",
              surface_root / "checked/CDS_energies.json"):
        inputs[str(p)] = api.digest(p)
    producer = dict(old["producer_sha256"])
    for p in (Path(__file__).resolve(), fine_script, surface_script):
        producer[str(p)] = api.digest(p)
    for task_file in (out / "live_inputs").glob("*/task.json"):
        inputs[str(task_file)] = api.digest(task_file)
    plan = dict(old, input_sha256=inputs, producer_sha256=producer,
        raw_root=str(ROOT / "calcs/validation" / out.name),
        source_campaign=str(source), source_state_read_only=str(state),
        force_validation_policy=POLICY, prepared=time.time(),
        recovery_scope="Unchanged saved geometry/active energy; separately reviewed CDS-explained four-step force validation",
        new_optimization_performed=False, new_hybrid_energy_calculation_performed=False)
    atomic_write_json(out / "plan.json", plan)
    atomic_write_json(out / "qualification_review.json", dict(fine=fine, CDS_attribution=surface,
        original_coarse_check=coarse, policy=POLICY, source_state=str(state),
        energy_only_qualification=True, Hessian_minimum_tested=False, native_G_qualified=False))
    new_report = dict(report, conservative_total_force_check=checks,
        gradient_max_Ha_Bohr=maximum, force_validation_policy=POLICY,
        original_failed_workflow=str(outer_path), original_failed_workflow_preserved=True,
        minimum_Hessian_tested=False, free_energy_computed=False)
    new_report["output_sha256"] = dict(report["output_sha256"], **inputs,
        **{str(out / "qualification_review.json"): api.digest(out / "qualification_review.json")})
    # Include measured cavity evidence for the normal native host forecast.
    execution = Path(old["raw_root"]) / "execution"
    witness = {str(p): api.digest(p) for p in execution.rglob("status.json")}
    new_report["output_sha256"].update(witness)
    path = out / "energy/tasks" / task["key"] / "status.json"
    atomic_write_json(path, new_report)
    atomic_write_json(out / "energy/done" / (task["key"]+".json"),
        dict(key=task["key"], status_path=str(path), status_sha256=api.digest(path),
             new_calculation_performed=False, validation_policy=POLICY["id"]))
    atomic_write_json(out / "workspace_energy_status.json", dict(
        status="stationary_electronic_case_checked_native_G_pending", energy_receipt=str(path),
        energy_receipt_sha256=api.digest(path), execution_sha256=witness,
        original_failed_workflow=str(outer_path), original_failed_workflow_preserved=True,
        validation_policy=POLICY, free_energy_computed=False))
    records, invalid = api.admitted_records(out, "energy")
    if invalid or len(records) != 1:
        raise ValueError("New scoped energy receipt failed normal identity/hash admission")
    print(json.dumps(dict(status="saved_energy_separately_qualified_native_G_pending",
        validation_policy=POLICY["id"], original_source_unchanged=True, E_smd_eV=report["E_smd_eV"])))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--fine", type=Path, required=True)
    p.add_argument("--surface", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    prepare(args)


if __name__ == "__main__":
    main()
