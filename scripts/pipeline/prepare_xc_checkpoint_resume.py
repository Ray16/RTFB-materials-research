#!/usr/bin/env python
"""Resume a stopped XC-workspace failure from its checked evaluated geometry.

The checkpoint remains unqualified as a minimum. The full optimizer, electronic
energy, structural/derivative checks and native thermal gates still run normally.
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


def prepare(parent, diagnostic, controls_root, out):
    check = json.loads((diagnostic / "gradient_qualification.json").read_text())
    point = api.read_plan(diagnostic)
    api.verify_sources(point["producer_sha256"])
    api.verify_sources(point["input_sha256"])
    if (check["status"] != "fixed_geometry_gradient_and_sampled_derivatives_passed"
            or check["free_energy_computed"] or check["native_Hessian_qualified"]
            or point["original_failed_campaign"] != str(parent)
            or not check["checks"]["passed"]):
        raise ValueError("Corrected checkpoint derivative evidence is incomplete or mismatched")
    api.verify_sources(check["output_sha256"])
    ownership = api.load_module("XC_resume_ownership", ROOT / "scripts/pipeline/prepare_polymer_postpilot_continuation.py")
    owners = ownership.verify_stopped(parent)
    factory = api.load_module("XC_resume_campaign_factory", ROOT / "scripts/pipeline/workspace_polymer_campaign.py")
    checker = api.load_module("XC_resume_control_pair", ROOT / "scripts/validation/qualify_xc_checkpoint.py")
    factory.SCRIPT = ROOT / "scripts/pipeline/xc_workspace_campaign.py"
    factory.ADAPTER = ROOT / "scripts/pipeline/xc_workspace_execution.py"
    factory.CONTROL_CHECKER = ROOT / "scripts/validation/qualify_xc_checkpoint.py"
    factory.controls = lambda plan: checker.checked_controls(controls_root, plan["execution_workspace_GiB"])
    old = api.read_plan(parent)
    controls = factory.controls(old)
    failure = parent / "energy/tasks" / old["tasks"][0]["key"] / "status.json"
    if json.loads(failure.read_text())["status"] != "failed_preserved":
        raise ValueError("Failed source must remain preserved")
    changes = dict(SMD_memory_qualification_controls=[str(controls_root / "controls" / str(i) / "status.json") for i in (1, 2)],
        purpose="Resume exact declared optimization from an evaluated checkpoint; all final checks remain required",
        restart_checkpoint=dict(diagnostic=str(diagnostic), evaluation=point["source_checkpoint"]["evaluation"],
            minimum_qualified=False, optimizer_state_restarted=True), native_pipeline_enabled=True)
    witnesses = {str(p): api.digest(p) for p in (diagnostic / "plan.json", diagnostic / "gradient_qualification.json", failure)}
    new = factory.fresh_campaign(parent, out, changes, witnesses)
    key = new["tasks"][0]["key"]
    ready_path = out / "live_inputs" / key / "task.json"
    ready = json.loads(ready_path.read_text())
    xyz = Path(ready["input_xyz"])
    shutil.copyfile(point["input_xyz"], xyz)
    ready["input_sha256"][str(xyz)] = api.digest(xyz)
    ready["source_restart_checkpoint"] = changes["restart_checkpoint"]
    atomic_write_json(ready_path, ready)
    new["input_sha256"].update(ready["input_sha256"])
    new["input_sha256"][str(ready_path)] = api.digest(ready_path)
    pilot = parent / "energy/pilot/status.json"
    original = json.loads(pilot.read_text())
    if original["status"] != "memory_pilot_passed":
        raise ValueError("Same-graph hybrid resource profile is not available")
    api.verify_sources(original["output_sha256"])
    copied = dict(original, resource_profile_reuse=dict(original_status=str(pilot), original_sha256=api.digest(pilot),
        corrected_gradient_status=str(diagnostic / "gradient_qualification.json"),
        corrected_gradient_sha256=api.digest(diagnostic / "gradient_qualification.json"),
        new_hybrid_calculation_performed=False,
        scope="Original same-graph hybrid footprint plus separately checked corrected checkpoint gradient; new full-state energies still required"))
    atomic_write_json(out / "energy/pilot/status.json", copied)
    shutil.copyfile(parent / "energy/pilot/energy_and_gradient.json", out / "energy/pilot/energy_and_gradient.json")
    for p in (out / "energy/pilot/status.json", out / "energy/pilot/energy_and_gradient.json", Path(point["input_xyz"])):
        new["input_sha256"][str(p)] = api.digest(p)
    for p in (Path(__file__).resolve(), ROOT / "scripts/pipeline/workspace_polymer_campaign.py"):
        new["producer_sha256"][str(p)] = api.digest(p)
    new["producer_sha256"].update(controls)
    atomic_write_json(out / "plan.json", new)
    atomic_write_json(out / "parent_ownership_verification.json", dict(checked=time.time(), hosts=owners))
    print(json.dumps(dict(status="XC_checkpoint_full_resume_prepared", evaluation=point["source_checkpoint"]["evaluation"],
                         final_geometry_qualified=False, native_G_computed=False)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--parent", type=Path, required=True)
    p.add_argument("--diagnostic", type=Path, required=True)
    p.add_argument("--controls", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    prepare(args.parent.resolve(), args.diagnostic.resolve(), args.controls.resolve(), args.out.resolve())


if __name__ == "__main__":
    main()
