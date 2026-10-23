#!/usr/bin/env python
"""Check corrected XC workspace handling at a preserved failed-run checkpoint.

Only geometry-level SCFs and sampled total-energy derivatives are evaluated.
No hybrid energy, optimized state, native Hessian or free energy is substituted.
"""
import argparse
import json
import math
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

ADAPTER = ROOT / "scripts/pipeline/xc_workspace_execution.py"


def checked_controls(folder, budget):
    paths = [str(folder / "controls" / str(i) / "status.json") for i in (1, 2)]
    checker = api.load_module("XC_control_identity", ROOT / "scripts/pipeline/run_smd_memory_recovery.py")
    checker.ADAPTER = ADAPTER
    sources = checker.checked_controls(dict(SMD_memory_qualification_controls=paths))
    for value in paths:
        p = Path(value).parent / "compact/execution/status.json"
        proof = json.loads(p.read_text())
        if (proof["budget_GiB"] != budget or proof.get("XC_equations_changed") is not False
                or proof.get("direct_LU_solver_retained") is not True
                or (json.loads(Path(value).read_text())["mult"] == 1 and proof.get("XC_density_batch_releases", 0) < 1)):
            raise ValueError("XC temporary-lifetime execution was not exercised by matched controls")
        api.verify_sources(proof["source_sha256"])
        sources[str(p)] = api.digest(p)
    return sources


def selected_checkpoint(parent):
    plan = api.read_plan(parent)
    task = plan["tasks"][0]
    folder = Path(plan["raw_root"]) / "energy" / task["key"] / "rest/geometry_checkpoints"
    trajectory, history = folder / "evaluated_trajectory.xyz", folder / "evaluations.jsonl"
    reader = api.load_module("checkpoint_plain_XYZ_reader", ROOT / "scripts/pipeline/sample_polymer_crest.py")
    frames = list(reader.xyz_frames(trajectory))
    rows = [json.loads(line) for line in history.read_text().splitlines()]
    if len(frames) != len(rows) or not rows:
        raise ValueError("Checkpoint geometry/evaluation inventory differs")
    for frame, row in zip(frames, rows):
        if (frame.splitlines()[1].strip() != f"Unqualified SCF evaluation {row['evaluation']}"
                or not row["SCF_converged"] or not math.isfinite(row["energy_Ha"])
                or not math.isfinite(row["gradient_max_Ha_Bohr"])):
            raise ValueError("Checkpoint frame, energy or SCF identity differs")
    index = min(range(len(rows)), key=lambda i: rows[i]["energy_Ha"])
    return frames[index], rows[index], {str(p): api.digest(p) for p in (trajectory, history)}


def prepare(args):
    parent = args.parent.resolve()
    old = api.read_plan(parent)
    api.verify_sources(old["producer_sha256"])
    api.verify_sources(old["input_sha256"])
    task = old["tasks"][0]
    failure = parent / "energy/tasks" / task["key"] / "status.json"
    data = json.loads(failure.read_text())
    if data["status"] != "failed_preserved" or "XC energy first derivative" not in data.get("error", ""):
        raise ValueError("An observed XC-gradient workspace failure is required")
    owner_check = api.load_module("XC_failed_owner_check", ROOT / "scripts/pipeline/prepare_polymer_postpilot_continuation.py")
    owners = owner_check.verify_stopped(parent)
    frame, row, witnesses = selected_checkpoint(parent)
    controls = checked_controls(args.controls.resolve(), old["execution_workspace_GiB"])
    out = args.campaign.resolve()
    out.mkdir(parents=True, exist_ok=False)
    xyz = out / "input.xyz"
    xyz.write_text(frame)
    from redox.validation.protonation import hydrogen_topology_status
    from redox.build.geometry import geometry_stereo_status
    if hydrogen_topology_status(task["smiles"], xyz) != "ok" or geometry_stereo_status(task["smiles"], xyz) not in ("ok", "not_applicable"):
        raise ValueError("Checkpoint has changed the declared graph/stereochemistry")
    reference = out / "checkpoint_reference.json"
    atomic_write_json(reference, dict(gradient_max_Ha_Bohr=row["gradient_max_Ha_Bohr"],
        evaluation=row["evaluation"], scope="Previously evaluated geometry-level gradient; not a hybrid energy or qualified minimum"))
    inputs = dict(old["input_sha256"], **witnesses)
    for p in (parent / "plan.json", failure, xyz, reference):
        inputs[str(p)] = api.digest(p)
    sources = dict(old["producer_sha256"], **controls)
    for p in (Path(__file__).resolve(), ADAPTER, ROOT / "scripts/validation/qualify_polymer_gradient_workspace.py",
              ROOT / "scripts/pipeline/prepare_polymer_postpilot_continuation.py", ROOT / "scripts/pipeline/sample_polymer_crest.py"):
        sources[str(p)] = api.digest(p)
    new = dict(old, input_sha256=inputs, producer_sha256=sources,
        raw_root=str(ROOT / "calcs/validation" / out.parent.name / out.name), input_xyz=str(xyz),
        previous_reference=str(reference), gradient_controls=str(args.controls.resolve()),
        minimum_large_GPU_MiB=args.minimum_GPU_MiB, maximum_memory_pilot_MiB=args.maximum_peak_GPU_MiB,
        purpose="Isolated checkpoint gradient resource/derivative qualification; no full-workflow or native-G admission",
        source_checkpoint=row, original_failed_campaign=str(parent), large_geometry_workflow_qualified=False)
    atomic_write_json(out / "plan.json", new)
    atomic_write_json(out / "parent_ownership_verification.json", dict(checked=time.time(), hosts=owners))
    print(json.dumps(dict(status="XC_checkpoint_gradient_check_prepared", selected_evaluation=row["evaluation"],
        minimum_device_MiB=args.minimum_GPU_MiB, full_workflow_qualified=False)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("prepare", "run"))
    p.add_argument("--campaign", type=Path, required=True)
    p.add_argument("--parent", type=Path)
    p.add_argument("--controls", type=Path)
    p.add_argument("--minimum-GPU-MiB", type=int, default=31744)
    p.add_argument("--maximum-peak-GPU-MiB", type=int, default=32768)
    args = p.parse_args()
    if args.action == "prepare":
        if not args.parent or not args.controls:
            p.error("Preparation requires failed parent and checked controls")
        from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
        with allocation_guard():
            host_capacity(1)
            os.sched_setaffinity(0, cpu_allocation(1))
        prepare(args)
    else:
        module = api.load_module("unchanged_checkpoint_gradient_checks", ROOT / "scripts/validation/qualify_polymer_gradient_workspace.py")
        module.ADAPTER = ADAPTER
        module.checked_controls = checked_controls
        module.run(args)


if __name__ == "__main__":
    main()
