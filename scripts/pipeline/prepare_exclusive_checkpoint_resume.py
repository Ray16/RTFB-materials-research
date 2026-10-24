#!/usr/bin/env python
"""Prepare a fresh full continuation after a verified device-sharing stop.

Reuse a complete evaluated checkpoint as an unqualified warm start. The prior
same-graph resource profile is not a new energy or exact-checkpoint benchmark.
All final optimization, energy, force and thermal gates remain mandatory.
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


def prepare(args):
    parent, out = args.parent.resolve(), args.out.resolve()
    old = api.read_plan(parent)
    api.verify_sources(old["producer_sha256"])
    api.verify_sources(old["input_sha256"])
    incident = json.loads(args.incident.read_text())
    if (incident["status"] != "owned_worker_stopped_checkpoints_preserved"
            or incident["worker_alive"] or incident["unrelated_processes_signaled"]
            or not incident["unrelated_processes_still_present"]):
        raise ValueError("Verified owned-worker stop and preserved peer processes required")
    if len(old["tasks"]) != 1 or (old["tasks"][0]["charge"], old["tasks"][0]["mult"]) != (0,1):
        raise ValueError("One existing neutral-singlet task is required")
    key = old["tasks"][0]["key"]
    source_worker = parent / "energy/workers" / f"{incident['host']}_gpu{incident['GPU']}_queue.json"
    worker = json.loads(source_worker.read_text())
    if (worker["pid"] != incident["pid"] or worker["process_start"] != incident["process_start"]
            or worker["status"] != "stopped_preserving_outputs"):
        raise ValueError("Resource incident does not identify the actual stopped source worker")
    owner_check = api.load_module("exclusive_resume_owner_check", ROOT / "scripts/pipeline/prepare_polymer_postpilot_continuation.py")
    ownership = owner_check.verify_stopped(parent)
    observer_check = api.load_module("exclusive_resume_monitor_identity", ROOT / "scripts/pipeline/run_polymer_spacing.py")
    monitor = json.loads(args.monitor.read_text())
    if not observer_check.owner_alive(monitor) or monitor["deadline"] <= time.time():
        raise ValueError("A live finite exclusivity observer is required")
    checkpoint = api.load_module("exclusive_resume_checkpoint_reader", ROOT / "scripts/validation/qualify_xc_checkpoint.py")
    frame, row, witnesses = checkpoint.selected_checkpoint(parent)
    factory = api.load_module("exclusive_resume_fresh_factory", ROOT / "scripts/pipeline/workspace_polymer_campaign.py")
    factory.SCRIPT = ROOT / "scripts/pipeline/xc_workspace_campaign.py"
    factory.ADAPTER = ROOT / "scripts/pipeline/xc_workspace_execution.py"
    factory.CONTROL_CHECKER = ROOT / "scripts/validation/qualify_xc_checkpoint.py"
    factory.controls = lambda plan: checkpoint.checked_controls(args.controls.resolve(), plan["execution_workspace_GiB"])
    witnesses.update({str(p): api.digest(p) for p in (args.incident.resolve(), source_worker)})
    changes = dict(recovery_scope="Fresh exclusive-device continuation from an evaluated checkpoint; original source and all final gates retained",
        interruption_kind="confirmed foreign GPU residency; not a molecular QC rejection",
        restart_checkpoint=dict(source_campaign=str(parent),evaluation=row["evaluation"],
            minimum_qualified=False,optimizer_state_restarted=True),
        # More conservative host allowance for possible screening changes at the
        # new geometry; the previous screened inventory is not relabeled exact.
        host_required_bytes=2*old["host_required_bytes"],
        checkpoint_resource_profile_exact_geometry=False)
    new = factory.fresh_campaign(parent,out,changes,witnesses)
    ready_path = out / "live_inputs" / key / "task.json"
    ready = json.loads(ready_path.read_text())
    xyz = Path(ready["input_xyz"])
    xyz.write_text(frame)
    from redox.validation.protonation import hydrogen_topology_status
    from redox.build.geometry import geometry_stereo_status
    topology = hydrogen_topology_status(ready["smiles"],xyz)
    stereo = geometry_stereo_status(ready["smiles"],xyz)
    if topology != "ok" or stereo not in ("ok","not_applicable"):
        raise ValueError("Checkpoint graph/stereochemistry does not match the declared model")
    ready["input_sha256"][str(xyz)] = api.digest(xyz)
    ready["source_restart_checkpoint"] = changes["restart_checkpoint"]
    atomic_write_json(ready_path,ready)
    new["input_sha256"].update(ready["input_sha256"])
    new["input_sha256"][str(ready_path)] = api.digest(ready_path)
    pilot = parent / "energy/pilot/status.json"
    resource = json.loads(pilot.read_text())
    if resource["status"] != "memory_pilot_passed":
        raise ValueError("Previous same-graph large-device resource qualification required")
    api.verify_sources(resource["output_sha256"])
    resource["resource_profile_reuse"] = dict(source=str(pilot),sha256=api.digest(pilot),
        new_hybrid_calculation_performed=False,exact_checkpoint_geometry_qualified=False,
        scope="Prior same-graph H100 hybrid footprint; unqualified checkpoint coordinates and doubled host allowance. Full-state energy and physical checks still required.")
    atomic_write_json(out / "energy/pilot/status.json",resource)
    shutil.copyfile(parent / "energy/pilot/energy_and_gradient.json",out / "energy/pilot/energy_and_gradient.json")
    atomic_write_json(out / "exclusivity_monitor_snapshot.json",monitor)
    atomic_write_json(out / "checkpoint_review.json",dict(checked=time.time(),source_evaluation=row,
        topology=topology,stereochemistry=stereo,unqualified_warm_start=True,
        parent_ownership=ownership,source_resource_incident=str(args.incident.resolve())))
    for p in (out / "energy/pilot/status.json",out / "energy/pilot/energy_and_gradient.json",
              out / "exclusivity_monitor_snapshot.json",out / "checkpoint_review.json"):
        new["input_sha256"][str(p)] = api.digest(p)
    for p in (Path(__file__).resolve(), ROOT / "scripts/fleet/watch_gpu_exclusivity.py"):
        new["producer_sha256"][str(p)] = api.digest(p)
    atomic_write_json(out / "plan.json",new)
    print(json.dumps(dict(status="exclusive_checkpoint_continuation_prepared",evaluation=row["evaluation"],
        host_required_GiB=new["host_required_bytes"]/2**30,minimum_qualified=False,G_computed=False)))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--parent",type=Path,required=True)
    p.add_argument("--incident",type=Path,required=True)
    p.add_argument("--monitor",type=Path,required=True)
    p.add_argument("--controls",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0,cpu_allocation(1))
    prepare(args)


if __name__=='__main__':
    main()
