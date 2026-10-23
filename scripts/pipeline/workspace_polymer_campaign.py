#!/usr/bin/env python
"""Resource-profile continuations using the existing polymer physics and queues.

Full calculations use fresh force storage. Checked energy parents can also be
imported into a separate thermal campaign without repeating their optimization.
All electronic, geometry and native-thermal acceptance checks remain unchanged.
"""
import argparse
import json
import math
import os
from pathlib import Path
import shutil
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

SCRIPT = Path(__file__).resolve()
ADAPTER = ROOT / "scripts/pipeline/smd_gradient_workspace.py"
CONTROL_CHECKER = ROOT / "scripts/validation/qualify_polymer_gradient_workspace.py"


def controls(plan):
    paths = [Path(p) for p in plan["SMD_memory_qualification_controls"]]
    roots = {p.parents[2] for p in paths}
    if len(paths) != 2 or len(roots) != 1:
        raise ValueError("One matched singlet/doublet control pair is required")
    checker = api.load_module("workspace_campaign_controls", CONTROL_CHECKER)
    return checker.checked_controls(roots.pop(), plan["execution_workspace_GiB"])


def profile_evidence(profile):
    plan = api.read_plan(profile)
    api.verify_sources(plan["producer_sha256"])
    api.verify_sources(plan["input_sha256"])
    controls(plan)
    path = profile / "footprint_qualification.json"
    data = json.loads(path.read_text())
    if (data["status"] not in ("fixed_geometry_workspace_and_numeric_continuity_qualified",
                              "fixed_geometry_workspace_qualified_larger_reference_unavailable")
            or "finished" not in data or data["free_energy_computed"]
            or data["optimized_geometry_qualified"] or data["native_Hessian_qualified"]
            or not 0 < data["observed_peak_GPU_MiB"] <= plan["maximum_memory_pilot_MiB"]):
        raise ValueError("Resource point is incomplete, over budget, or promoted beyond scope")
    pilot = profile / "energy/pilot/status.json"
    record = json.loads(pilot.read_text())
    if record["status"] != "memory_pilot_passed":
        raise ValueError("Complete energy/gradient resource evidence is required")
    api.verify_sources(record["output_sha256"])
    execution = Path(plan["raw_root"]) / "pilot_execution/status.json"
    used = json.loads(execution.read_text())
    if (used.get("gradient_solver_calls", 0) < 1 or used["budget_GiB"] != plan["execution_workspace_GiB"]
            or used.get("direct_LU_solver_retained") is not True
            or used["source_sha256"].get(str(ADAPTER)) != api.digest(ADAPTER)):
        raise ValueError("Qualified bounded-gradient profile was not exercised")
    api.verify_sources(used["source_sha256"])
    return plan, {str(p): api.digest(p) for p in (path, pilot, execution)}


def fresh_campaign(source, out, changes, witnesses):
    old = api.read_plan(source)
    api.verify_sources(old["producer_sha256"])
    api.verify_sources(old["input_sha256"])
    if len(old["tasks"]) != 1 or old["tasks"][0]["charge"] != 0 or old["tasks"][0]["mult"] != 1:
        raise ValueError("This execution branch requires one existing neutral-singlet task")
    out.mkdir(parents=True, exist_ok=False)
    shutil.copytree(source / "live_inputs", out / "live_inputs")
    shutil.copytree(source / "controls", out / "controls")
    for stage in ("energy", "thermal"):
        for name in ("tasks", "done", "claims", "failures", "workers", "pilot", "launches"):
            (out / stage / name).mkdir(parents=True)
    inputs = dict(old["input_sha256"], **witnesses)
    inputs[str(source / "plan.json")] = api.digest(source / "plan.json")
    for task in old["tasks"]:
        p = out / "live_inputs" / task["key"] / "task.json"
        ready = json.loads(p.read_text())
        ready["input_xyz"] = str(p.parent / "input.xyz")
        ready["input_sha256"][ready["input_xyz"]] = api.digest(ready["input_xyz"])
        atomic_write_json(p, ready)
        inputs[str(p)] = api.digest(p)
        inputs.update(ready["input_sha256"])
    sources = dict(old["producer_sha256"])
    for p in (SCRIPT, ADAPTER, CONTROL_CHECKER):
        sources[str(p)] = api.digest(p)
    new = dict(old, **changes)
    new.update(producer_sha256=sources, input_sha256=inputs,
        raw_root=str(ROOT / "calcs/polymer_models" / out.parent.name / out.name),
        execution_driver=str(SCRIPT), source_campaign=str(source), created=time.time(),
        recovery_scope="Qualified workspace execution with original task/model/QC and separately stored full calculations")
    new["producer_sha256"].update(controls(new))
    return new


def prepare_profile(profile, out):
    plan, evidence = profile_evidence(profile)
    new = fresh_campaign(profile, out, dict(host_required_bytes=json.loads(
        (profile / "footprint_qualification.json").read_text())["host_required_bytes"],
        purpose="Full optimization/active energy and separately checked native thermochemistry for the exact declared graph",
        energy_only_execution=False, native_pipeline_enabled=True), evidence)
    pilot = profile / "energy/pilot/status.json"
    copied = dict(json.loads(pilot.read_text()), resource_profile_reuse=dict(
        original_status=str(pilot), original_sha256=api.digest(pilot),
        new_calculation_performed=False, scope="Same task/model/adapter; fresh optimization storage"))
    atomic_write_json(out / "energy/pilot/status.json", copied)
    energy = profile / "energy/pilot/energy_and_gradient.json"
    shutil.copyfile(energy, out / "energy/pilot/energy_and_gradient.json")
    for p in (out / "energy/pilot/status.json", out / "energy/pilot/energy_and_gradient.json"):
        new["input_sha256"][str(p)] = api.digest(p)
    atomic_write_json(out / "plan.json", new)
    return new


def checked_energy(source):
    from redox.core.protocol import load_record
    old = api.read_plan(source)
    api.verify_sources(old["producer_sha256"])
    api.verify_sources(old["input_sha256"])
    task = old["tasks"][0]
    path = source / "energy/tasks" / task["key"] / "status.json"
    report = json.loads(path.read_text())
    if (report["status"] != "state_records_verified" or report.get("physical_QC") != "passed"
            or not report.get("conservative_total_force_check", {}).get("passed")
            or any(report["task"].get(k) != v for k, v in task.items())):
        raise ValueError("A checked stationary electronic/physical parent is required")
    api.verify_sources(report["output_sha256"])
    active = load_record(Path(report["state_dir"]), task["charge"], task["mult"])
    if active is None or not math.isfinite(report["E_smd_eV"]) or abs(active["e_smd_eV"] - report["E_smd_eV"]) > 1e-8:
        raise ValueError("Electronic parent does not reconstruct on its exact geometry")
    workflow = source / "workspace_energy_status.json"
    if workflow.exists():
        d = json.loads(workflow.read_text())
        if (d["status"] != "stationary_electronic_case_checked_native_G_pending"
                or d["energy_receipt_sha256"] != api.digest(path)):
            raise ValueError("Direct energy workflow has not completed its final checks")
    return report, path


def prepare_native(source, out, control_root, budget, minimum_gpu):
    if not math.isfinite(budget) or budget <= 0 or minimum_gpu < (budget + 1) * 1024:
        raise ValueError("Positive workspace budget and physical device headroom are required")
    report, path = checked_energy(source)
    changes = dict(execution_workspace_GiB=budget, minimum_large_GPU_MiB=minimum_gpu,
        maximum_memory_pilot_MiB=int(budget * 1024), host_required_bytes=256 * 2**30,
        SMD_memory_qualification_controls=[str(control_root / "controls" / str(i) / "status.json") for i in (1, 2)],
        inherited_energy_only=True, energy_only_execution=False, native_pipeline_enabled=True,
        purpose="Native thermochemistry from the exact checked electronic parent; no repeated optimization")
    # Validate the new execution controls before creating any output tree.
    controls(dict(api.read_plan(source), **changes))
    new = fresh_campaign(source, out, changes, {str(path): api.digest(path)})
    receipt = out / "energy/done" / (report["task"]["key"] + ".json")
    atomic_write_json(receipt, dict(key=report["task"]["key"], status_path=str(path),
        status_sha256=api.digest(path), reused_checked_energy=True, new_calculation_performed=False))
    atomic_write_json(out / "plan.json", new)
    records, invalid = api.admitted_records(out, "energy")
    if invalid or len(records) != 1:
        raise ValueError("Imported electronic parent failed normal campaign admission")
    return new


def annotate(path, directory, witnesses, plan):
    report = json.loads(path.read_text())
    report["workspace_execution"] = dict(adapter_sha256=api.digest(ADAPTER),
        workspace_GiB=plan["execution_workspace_GiB"], controls=witnesses,
        scientific_model_changed=False, precision_changed=False, grids_changed=False)
    report["output_sha256"].update({str(p): api.digest(p) for p in directory.rglob("*") if p.is_file()})
    report["output_sha256"].update(witnesses)
    atomic_write_json(path, report)


def execution():
    driver = api.load_module("workspace_standard_polymer_driver", ROOT / "scripts/pipeline/refine_polymer_conformers.py")
    driver.SCRIPT = SCRIPT
    campaign = Path(sys.argv[sys.argv.index("--campaign") + 1]).resolve()
    plan = api.read_plan(campaign)
    stage = sys.argv[sys.argv.index("--stage") + 1] if "--stage" in sys.argv else "energy"
    if sys.argv[1] == "launch":
        original_slots = driver.additional_host_slots
        required_gib = plan.get("host_required_bytes", 32 * 2**30) / 2**30 if stage == "energy" else 256
        driver.additional_host_slots = lambda available, promises, budget: original_slots(
            available, promises, max(budget, required_gib))
        atomic_write_json(campaign / stage / "workspace_host_budget.json", dict(
            required_host_GiB=required_gib, scope="Measured DF resource inventory or declared native-thermal budget"))
    if sys.argv[1] == "worker":
        from redox.workflows.resources import verify_gate_owner, host_capacity
        verify_gate_owner()
        host_capacity(2, minimum_ram_gib=plan.get("host_required_bytes", 32 * 2**30) / 2**30
                      if stage == "energy" else 256)
        witnesses = controls(plan)
        adapter = api.load_module("workspace_campaign_adapter", ADAPTER)
        old_energy, old_native = driver.calculate_energy, api.calculate_native_G

        def energy(folder, task, pilot=False):
            if pilot:
                raise ValueError("This branch reuses its already checked external energy resource point")
            directory = Path(plan["raw_root"]) / "workspace_execution" / (task["key"] + "_energy")
            with adapter.compact_smd_execution(directory, budget_GiB=plan["execution_workspace_GiB"]):
                result = old_energy(folder, task, pilot=False)
            annotate(folder / "energy/tasks" / task["key"] / "status.json", directory, witnesses, plan)
            return result

        def native(folder, record):
            directory = Path(plan["raw_root"]) / "workspace_execution" / (record["task"]["key"] + "_native")
            with adapter.compact_smd_execution(directory, budget_GiB=plan["execution_workspace_GiB"]):
                result = old_native(folder, record)
            annotate(folder / "thermal/tasks" / record["task"]["key"] / "status.json", directory, witnesses, plan)
            return result

        driver.calculate_energy = energy
        api.calculate_native_G = native
    driver.main()


def main():
    if len(sys.argv) > 1 and sys.argv[1] in ("worker", "launch", "collect"):
        execution()
        return
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("prepare-profile", "prepare-native"))
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--controls", type=Path)
    p.add_argument("--workspace-GiB", type=float, default=64)
    p.add_argument("--minimum-GPU-MiB", type=int, default=81920)
    args = p.parse_args()
    from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    if args.action == "prepare-profile":
        plan = prepare_profile(args.source.resolve(), args.out.resolve())
    else:
        if not args.controls:
            p.error("Native handoff requires its checked execution controls")
        plan = prepare_native(args.source.resolve(), args.out.resolve(), args.controls.resolve(),
                              args.workspace_GiB, args.minimum_GPU_MiB)
    print(json.dumps(dict(status="workspace_campaign_prepared", tasks=len(plan["tasks"]),
                         inherited_energy=plan.get("inherited_energy_only", False))))


if __name__ == "__main__":
    main()
