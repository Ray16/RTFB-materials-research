#!/usr/bin/env python
"""Prepare/run a separate fine-step diagnostic with checked memory controls.

The original failed state remains unaccepted. This wrapper does not make an
energy receipt or change the production geometry/Hessian protocol.
"""
import argparse
from contextlib import redirect_stderr, redirect_stdout
import json
import os
from pathlib import Path
import shutil
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

ADAPTER = ROOT / "scripts/pipeline/xc_workspace_execution.py"
DIAGNOSTIC = ROOT / "scripts/validation/polymer_gradient_diagnostic.py"
CONTROLS = ROOT / "scripts/validation/qualify_xc_checkpoint.py"
THREADS = ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "MKL_NUM_THREADS",
           "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "XTB_THREADS")


def prepare(args):
    source = args.source.resolve()
    plan = api.read_plan(source)
    api.verify_sources(plan["producer_sha256"])
    api.verify_sources(plan["input_sha256"])
    workflow = source / "workspace_energy_status.json"
    failure = json.loads(workflow.read_text())
    if failure["status"] != "energy_attempt_failed_preserved" or "total force" not in failure["error"]:
        raise ValueError("An observed final total-force failure is required")
    if len(plan["tasks"]) != 1 or (plan["tasks"][0]["charge"], plan["tasks"][0]["mult"]) != (0, 1):
        raise ValueError("Only one existing neutral-singlet graph is in scope")
    key = plan["tasks"][0]["key"]
    state = Path(plan["raw_root"]) / "energy" / key / "rest"
    coarse = state / "total_energy_force_check/summary.json"
    old_check = json.loads(coarse.read_text())
    if old_check["passed"] or old_check["geometry_sha256"] != api.digest(state / "opt.xyz"):
        raise ValueError("Failed check or saved geometry identity differs")
    checker = api.load_module("fine_step_XC_controls", CONTROLS)
    controls = checker.checked_controls(args.controls.resolve(), 28)
    out = args.campaign.resolve()
    out.mkdir(parents=True, exist_ok=False)
    xyz = out / "input.xyz"
    shutil.copyfile(state / "opt.xyz", xyz)
    raw = ROOT / "calcs/validation" / out.name
    command = [str(DIAGNOSTIC), "--xyz", str(xyz), "--out", str(raw / "diagnostic"),
        "--grid-level", "3", "--steps", ".0001", ".00005", ".00001",
        "--baseline-repeats", "3", "--top-components", "2", "--omit-reference-axis", "--conservative-cds"]
    for atom, axis in sorted({(r["atom"], r["axis"]) for r in old_check["checks"]}):
        command.extend(["--probe", f"{atom}:{axis}"])
    inputs = {str(p): api.digest(p) for p in (source / "plan.json", workflow, coarse,
        state / "opt.xyz", state / "total_energy_force_check/baseline.json", xyz)}
    producers = dict(plan["producer_sha256"], **controls)
    producers.update({str(p): api.digest(p) for p in (Path(__file__).resolve(), ADAPTER, DIAGNOSTIC, CONTROLS)})
    new = dict(source_campaign=str(source), task=plan["tasks"][0], input_sha256=inputs,
        producer_sha256=producers, control_root=str(args.controls.resolve()), command=command,
        raw_root=str(raw), execution_workspace_GiB=28, minimum_GPU_MiB=31744,
        host_budget_GiB=max(32, plan["host_required_bytes"]/2**30),
        maximum_peak_GPU_MiB=32768, original_coarse_force_check_preserved=True,
        production_admission_changed=False, prepared=time.time(),
        acceptance=dict(maximum_total_derivative_error_Ha_Bohr=1e-6,
            finest_pair_derivative_difference_Ha_Bohr=1e-6, finest_steps_Bohr=[0.00005, 0.00001],
            maximum_baseline_energy_spread_Ha=1e-11,
            maximum_repeated_gradient_difference_Ha_Bohr=1e-7),
        scope="Separate fixed-geometry derivative-convergence diagnostic; no production state or G admission")
    atomic_write_json(out / "plan.json", new)
    print(json.dumps(dict(status="fine_step_diagnostic_prepared", campaign=str(out), source_state=str(state))))


def run(args):
    from redox.workflows.resources import verify_gate_owner, host_capacity, gpu_suitability
    verify_gate_owner()
    plan = json.loads((args.campaign / "plan.json").read_text())
    if len(os.sched_getaffinity(0)) != 2 or any(os.environ.get(k) != "2" for k in THREADS):
        raise ValueError("Two disjoint CPU cores and all native caps=2 required")
    host_capacity(2, minimum_ram_gib=plan["host_budget_GiB"])
    if gpu_suitability(int(os.environ["CUDA_VISIBLE_DEVICES"]))["memory_MiB"] < plan["minimum_GPU_MiB"]:
        raise ValueError("Device is below the separately controlled workspace profile")
    api.verify_sources(plan["producer_sha256"])
    api.verify_sources(plan["input_sha256"])
    checker = api.load_module("fine_step_runtime_controls", CONTROLS)
    checker.checked_controls(Path(plan["control_root"]), plan["execution_workspace_GiB"])
    raw = Path(plan["raw_root"])
    raw.mkdir(parents=True, exist_ok=False)
    status = dict(status="fine_step_diagnostic_running", host=os.uname().nodename, pid=os.getpid(),
        process_start=Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()[19],
        CPU_affinity=sorted(os.sched_getaffinity(0)), GPU=os.environ["CUDA_VISIBLE_DEVICES"],
        started=time.time(), production_admission_changed=False)
    target = args.campaign / "status.json"
    atomic_write_json(target, status)
    monitor = None
    try:
        adapter = api.load_module("fine_step_execution_adapter", ADAPTER)
        engine = api.load_module("fine_step_memory_observer", ROOT / "scripts/pipeline/run_pi_stacking.py")
        module = api.load_module("fine_step_original_diagnostic", DIAGNOSTIC)
        with engine.MemoryMonitor() as monitor, (raw / "calculation.log").open("x") as log, redirect_stdout(log), redirect_stderr(log):
            with adapter.compact_smd_execution(raw / "execution", budget_GiB=plan["execution_workspace_GiB"]):
                before = sys.argv
                try:
                    sys.argv = plan["command"]
                    module.main()
                finally:
                    sys.argv = before
        report = json.loads((raw / "diagnostic/status.json").read_text())
        if report["status"] != "diagnostic_complete_no_promotion":
            raise ValueError("Underlying diagnostic did not complete")
        api.verify_sources(plan["producer_sha256"])
        api.verify_sources(plan["input_sha256"])
        status.update(status="fine_step_diagnostic_completed_independent_review_pending",
            output_sha256={str(p): api.digest(p) for p in raw.rglob("*") if p.is_file()})
    except BaseException as error:
        status.update(status="fine_step_diagnostic_failed_preserved", error=repr(error))
        raise
    finally:
        status.update(finished=time.time(), observed_peak_GPU_MiB=monitor.peak if monitor else None)
        atomic_write_json(target, status)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("action", choices=("prepare", "run"))
    p.add_argument("--campaign", type=Path, required=True)
    p.add_argument("--source", type=Path)
    p.add_argument("--controls", type=Path)
    args = p.parse_args()
    if args.action == "prepare":
        if not args.source or not args.controls:
            p.error("Preparation requires the failed source and checked controls")
        from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
        with allocation_guard():
            host_capacity(1)
            os.sched_setaffinity(0, cpu_allocation(1))
        prepare(args)
    else:
        run(args)


if __name__ == "__main__":
    main()
