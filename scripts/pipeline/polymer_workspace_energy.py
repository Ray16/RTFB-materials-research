#!/usr/bin/env python
"""Qualified 32 GB execution branch for one existing polymer energy task.

Use the original optimizer, checkpointing, active energy and physical checks.
Native thermochemistry is deliberately a subsequent, separately qualified stage.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

ADAPTER = ROOT / "scripts/pipeline/smd_gradient_workspace.py"
QUALIFIER = ROOT / "scripts/validation/qualify_polymer_gradient_workspace.py"


def gradient_evidence(campaign):
    plan = api.read_plan(campaign)
    api.verify_sources(plan["producer_sha256"])
    api.verify_sources(plan["input_sha256"])
    path = campaign / "gradient_qualification.json"
    data = json.loads(path.read_text())
    if (data["status"] != "fixed_geometry_gradient_and_sampled_derivatives_passed"
            or data["free_energy_computed"] or data["native_Hessian_qualified"]
            or data["observed_peak_GPU_MiB"] > plan["maximum_memory_pilot_MiB"]
            or data["previous_gradient_max_difference"] > 1e-6
            or not data["checks"]["passed"]
            or data["checks"]["acceptance_Ha_Bohr"] != 1e-6):
        raise ValueError("Large-point gradient/resource evidence is not qualified")
    checks = data["checks"]["checks"]
    coordinates = {(r["atom"], r["axis"]) for r in checks}
    if len(checks) != 4 or len(coordinates) != 2:
        raise ValueError("Two components and two displacement sizes are required")
    for atom, axis in coordinates:
        rows = [r for r in checks if (r["atom"], r["axis"]) == (atom, axis)]
        if {r["step_Bohr"] for r in rows} != {.0005, .001}:
            raise ValueError("Derivative displacement steps differ")
    if any(not 0 <= r["absolute_error"] <= 1e-6 for r in checks):
        raise ValueError("Sampled derivative error exceeds the unchanged gate")
    api.verify_sources(data["output_sha256"])
    return plan, {str(campaign / "plan.json"): api.digest(campaign / "plan.json"),
                  str(path): api.digest(path)}


def prepare(args):
    point = args.gradient.resolve()
    proof, evidence = gradient_evidence(point)
    parent = Path(proof["parent_failed_point"])
    old = api.read_plan(parent)
    inputs = dict(old["input_sha256"], **proof["input_sha256"], **evidence)
    controls = api.load_module("gradient_control_verifier", QUALIFIER)
    witnesses = controls.checked_controls(Path(proof["gradient_controls"]), old["execution_workspace_GiB"])
    out = args.campaign.resolve()
    out.mkdir(parents=True, exist_ok=False)
    shutil.copytree(parent / "live_inputs", out / "live_inputs")
    shutil.copytree(parent / "controls", out / "controls")
    for stage in ("energy", "thermal"):
        for name in ("tasks", "done", "claims", "failures", "workers", "pilot", "launches"):
            (out / stage / name).mkdir(parents=True)
    for task in old["tasks"]:
        p = out / "live_inputs" / task["key"] / "task.json"
        ready = json.loads(p.read_text())
        ready["input_xyz"] = str(p.parent / "input.xyz")
        ready["input_sha256"][ready["input_xyz"]] = api.digest(ready["input_xyz"])
        atomic_write_json(p, ready)
        inputs[str(p)] = api.digest(p)
        inputs.update(ready["input_sha256"])
    inventory = json.loads((parent / "footprint_qualification.json").read_text())
    sources = dict(proof["producer_sha256"], **witnesses)
    for p in (Path(__file__).resolve(), ADAPTER, QUALIFIER):
        sources[str(p)] = api.digest(p)
    plan = dict(old, raw_root=str(ROOT / "calcs/polymer_models" / out.parent.name / out.name),
        input_sha256=inputs, producer_sha256=sources, gradient_qualification=str(point),
        SMD_memory_qualification_controls=[str(Path(proof["gradient_controls"]) / "controls" / str(spin) / "status.json") for spin in (1, 2)],
        host_required_bytes=inventory["host_required_bytes"],
        execution_driver=str(Path(__file__).resolve()), energy_only_execution=True,
        recovery_scope="Bounded SMD workspace and gradient contractions; original tasks/model/optimizer/QC retained",
        full_large_geometry_workflow_qualified=False, native_Hessian_qualified=False,
        created=time.time())
    atomic_write_json(out / "plan.json", plan)
    print(json.dumps(dict(status="qualified_workspace_energy_task_prepared", tasks=len(plan["tasks"]))))


def run(args):
    import psutil
    from redox.workflows.resources import verify_gate_owner, host_capacity, gpu_suitability
    verify_gate_owner()
    caps = ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "XTB_THREADS")
    if len(os.sched_getaffinity(0)) != 2 or any(os.environ.get(k) != "2" for k in caps):
        raise ValueError("Two disjoint CPU cores and native thread caps=2 required")
    campaign = args.campaign.resolve()
    plan = api.read_plan(campaign)
    api.verify_sources(plan["producer_sha256"])
    api.verify_sources(plan["input_sha256"])
    gradient_evidence(Path(plan["gradient_qualification"]))
    if plan["host_required_bytes"] > .75 * psutil.virtual_memory().available:
        raise RuntimeError("Measured density-fitting host-RAM requirement is not available")
    host_capacity(2, minimum_ram_gib=plan["host_required_bytes"] / 2**30)
    if gpu_suitability(int(os.environ["CUDA_VISIBLE_DEVICES"]))["memory_MiB"] < plan["minimum_large_GPU_MiB"]:
        raise RuntimeError("Wrong GPU memory class")
    driver = api.load_module("qualified_workspace_original_energy", ROOT / "scripts/pipeline/refine_polymer_conformers.py")
    task = driver.backend(campaign, "energy").task_input(campaign, plan["tasks"][0])
    # The backend atomically creates this task directory before computing.
    # Never overwrite a previous task or adopt another worker's output.
    if (campaign / "energy/tasks" / task["key"]).exists():
        raise FileExistsError("Preserve the existing energy attempt")
    raw = Path(plan["raw_root"])
    raw.mkdir(parents=True, exist_ok=False)
    adapter = api.load_module("qualified_workspace_energy_adapter", ADAPTER)
    report = dict(status="optimization_starting", host=os.uname().nodename, pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19],
        CPU_affinity=sorted(os.sched_getaffinity(0)), GPU=os.environ["CUDA_VISIBLE_DEVICES"],
        started=time.time(), free_energy_computed=False, native_Hessian_qualified=False,
        synthesis_admission=False, sampled_ensemble_convergence_claimed=False)
    status = campaign / "workspace_energy_status.json"
    atomic_write_json(status, report)
    try:
        with adapter.compact_smd_execution(raw / "execution", budget_GiB=plan["execution_workspace_GiB"]):
            driver.calculate_energy(campaign, task, pilot=False)
        p = campaign / "energy/tasks" / task["key"] / "status.json"
        value = json.loads(p.read_text())
        if value["status"] != "state_records_verified" or value.get("physical_QC") != "passed":
            raise ValueError("Original energy/physical checks did not pass")
        api.verify_sources(value["output_sha256"])
        api.verify_sources(plan["producer_sha256"])
        api.verify_sources(plan["input_sha256"])
        report.update(status="stationary_electronic_case_checked_native_G_pending", energy_receipt=str(p),
            energy_receipt_sha256=api.digest(p), observed_peak_GPU_MiB=value["observed_peak_GPU_MiB"],
            execution_sha256={str(p): api.digest(p) for p in (raw / "execution").rglob("*") if p.is_file()})
    except BaseException as error:
        report.update(status="energy_attempt_failed_preserved", error=repr(error))
        raise
    finally:
        report["finished"] = time.time()
        atomic_write_json(status, report)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "run"))
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--gradient", type=Path)
    args = parser.parse_args()
    if args.action == "prepare":
        if not args.gradient:
            parser.error("Preparation requires the checked large-gradient campaign")
        from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
        with allocation_guard():
            host_capacity(1)
            os.sched_setaffinity(0, cpu_allocation(1))
        prepare(args)
    else:
        run(args)


if __name__ == "__main__":
    main()
