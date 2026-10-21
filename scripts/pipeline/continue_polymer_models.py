#!/usr/bin/env python
"""Continue checked polymer inputs without altering historical or live producers."""

import argparse
import json
import os
from pathlib import Path
import shlex
import shutil
import socket
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.polymer_continuation import execution_settings, review_memory_pilot

SCRIPT = Path(__file__).resolve()
OLD_SPACING = ROOT / "results/polymer_models/campaigns/2026-10-20_spacing_pilot"
SPACING = ROOT / "results/polymer_models/campaigns/2026-10-21_spacing_continuation"
OLD_RECOVERY = ROOT / "results/polymer_models/campaigns/2026-10-20_geometry_recovery"
RECOVERY = ROOT / "results/polymer_models/campaigns/2026-10-21_checkpoint_continuation"
KEY = "exp_phthalimide_control__n3__stereo1__open"
THREAD_KEYS = (
    "OMP_NUM_THREADS",
    "OMP_THREAD_LIMIT",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "XTB_THREADS",
)


def sources():
    paths = list((ROOT / "src/redox").rglob("*.py"))
    paths += list((ROOT / "config").glob("*.py")) + list(
        (ROOT / "config").glob("*.json")
    )
    paths += [
        SCRIPT,
        ROOT / "scripts/pipeline/run_polymer_spacing.py",
        ROOT / "scripts/pipeline/run_pi_stacking.py",
        ROOT / "scripts/pipeline/recover_pi_stacking.py",
        ROOT / "scripts/validation/audit_polymer_geometries.py",
    ]
    return {str(path): api.digest(path) for path in paths}


def prepare_spacing():
    if SPACING.exists():
        raise RuntimeError("Preserve the existing continuation")
    old = api.read_plan(OLD_SPACING)
    owner_check = api.load_module(
        "continuation_owner_check", ROOT / "scripts/pipeline/run_polymer_spacing.py"
    )
    owner = json.loads(
        (OLD_SPACING / "energy/workers/lambda13_gpu2_pilot.json").read_text()
    )
    if owner_check.owner_alive(owner):
        raise RuntimeError("Prior pilot owner is still live; preserve its outputs")
    api.verify_sources(old["producer_sha256"])
    api.verify_sources(old["input_sha256"])
    launch = json.loads(
        (OLD_SPACING / "energy/launches/lambda13_pilot.json").read_text()
    )
    source_status = OLD_SPACING / "energy/pilot/status.json"
    source_result = OLD_SPACING / "energy/pilot/energy_and_gradient.json"
    report, result = (
        json.loads(source_status.read_text()),
        json.loads(source_result.read_text()),
    )
    device = launch["workers"][0]["device"]["memory_MiB"]
    review = review_memory_pilot(report, result, device)
    SPACING.mkdir(parents=True)
    shutil.copyfile(OLD_SPACING / "cases.json", SPACING / "cases.json")
    tasks = []
    for task in old["tasks"]:
        directory = SPACING / "live_inputs" / task["key"]
        directory.mkdir(parents=True)
        previous = OLD_SPACING / "live_inputs" / task["key"] / "task.json"
        ready = json.loads(previous.read_text())
        api.verify_sources(ready["input_sha256"])
        shutil.copyfile(ready["input_xyz"], directory / "input.xyz")
        ready.update(
            input_xyz=str(directory / "input.xyz"),
            input_sha256={
                str(directory / "input.xyz"): api.digest(directory / "input.xyz"),
                str(previous): api.digest(previous),
            },
        )
        atomic_write_json(directory / "task.json", ready)
        tasks.append(task)
    inputs = [
        SPACING / "cases.json",
        OLD_SPACING / "plan.json",
        source_status,
        source_result,
        OLD_SPACING / "energy/launches/lambda13_pilot.json",
    ]
    inputs += list((SPACING / "live_inputs").glob("*/task.json"))
    prepared = dict(
        old,
        created=time.time(),
        tasks=tasks,
        raw_root=str(ROOT / "calcs/polymer_models" / SPACING.name),
        prepared_graphs=6,
        held_graphs=[],
        parent_campaign=str(OLD_SPACING),
        input_sha256={str(p): api.digest(p) for p in inputs},
        producer_sha256=sources(),
        maximum_memory_pilot_MiB=review["maximum_memory_pilot_MiB"],
        minimum_large_GPU_MiB=review["minimum_device_MiB"],
        execution_change=(
            "Release geometry SCF/scanner before active SP; "
            "same equations and acceptance gates"
        ),
        numerical_pilot_review=review,
    )
    atomic_write_json(SPACING / "preparation.json", prepared)
    api.freeze(SPACING)
    # The numerical pilot already completed. Preserve its failed old guard record
    # and reuse its actual energy/gradient evidence under the reviewed headroom.
    proof = dict(
        report,
        status="memory_pilot_passed",
        prior_guard_status=str(source_status),
        numerical_review=review,
        output_sha256={
            str(source_result): api.digest(source_result),
            str(source_status): api.digest(source_status),
        },
        scope="Numerical resource qualification only; geometry remains unoptimized",
    )
    atomic_write_json(SPACING / "energy/pilot/status.json", proof)
    plan = api.read_plan(SPACING)
    plan["producer_sha256"][str(SCRIPT)] = api.digest(SCRIPT)
    plan["input_sha256"][str(SPACING / "energy/pilot/status.json")] = api.digest(
        SPACING / "energy/pilot/status.json"
    )
    atomic_write_json(SPACING / "plan.json", plan)
    api.collect(SPACING)
    print(json.dumps(review))


def prepare_recovery():
    if RECOVERY.exists():
        raise RuntimeError("Preserve the existing checkpoint continuation")
    old = api.read_plan(OLD_RECOVERY)
    api.verify_sources(old["producer_sha256"])
    failure = OLD_RECOVERY / "failures" / f"{KEY}.json"
    record = json.loads(failure.read_text())
    owner_check = api.load_module(
        "continuation_owner_check", ROOT / "scripts/pipeline/run_polymer_spacing.py"
    )
    if owner_check.owner_alive(record["owner"]):
        raise RuntimeError("Original failed worker is still live")
    if "did not converge" not in record["error"]:
        raise ValueError(
            "This continuation addresses recorded optimizer exhaustion only"
        )
    if (OLD_RECOVERY / "claims" / KEY).exists():
        raise RuntimeError("The original task still has an owner")
    ready_path = OLD_RECOVERY / "live_inputs" / KEY / "task.json"
    ready = json.loads(ready_path.read_text())
    point = (
        Path(old["raw_root"]) / KEY / "rest/geometry_checkpoints/optimizer_returned.xyz"
    )
    audit = api.load_module(
        "continuation_audit", ROOT / "scripts/validation/audit_polymer_geometries.py"
    )
    sanity = audit.audit(ready, point)
    if sanity["issues"]:
        raise ValueError(f"Checkpoint geometry needs correction: {sanity['issues']}")
    for folder in (
        "claims",
        "tasks",
        "done",
        "failures",
        "workers",
        "launches",
        "pilot",
        "live_inputs",
    ):
        (RECOVERY / folder).mkdir(parents=True, exist_ok=True)
    directory = RECOVERY / "live_inputs" / KEY
    directory.mkdir()
    shutil.copyfile(point, directory / "input.xyz")
    ready.update(
        input_xyz=str(directory / "input.xyz"),
        checkpoint_continuation_of=str(point),
        input_sha256={
            str(point): api.digest(point),
            str(ready_path): api.digest(ready_path),
            str(directory / "input.xyz"): api.digest(directory / "input.xyz"),
        },
    )
    atomic_write_json(directory / "task.json", ready)
    atomic_write_json(RECOVERY / "checkpoint_geometry_review.json", sanity)
    shutil.copyfile(OLD_RECOVERY / "pilot/status.json", RECOVERY / "pilot/status.json")
    inputs = [
        failure,
        ready_path,
        point,
        directory / "task.json",
        RECOVERY / "pilot/status.json",
        OLD_RECOVERY / "plan.json",
        RECOVERY / "checkpoint_geometry_review.json",
    ]
    plan = dict(
        old,
        created=time.time(),
        tasks=[next(t for t in old["tasks"] if t["key"] == KEY)],
        raw_root=str(ROOT / "calcs/polymer_models" / RECOVERY.name),
        input_sha256={str(p): api.digest(p) for p in inputs},
        producer_sha256=sources(),
        recovery_settings=dict(
            max_opt_steps=200,
            geometry_SCF_conv_tol=1e-12,
            geometry_SCF_gradient_tol=1e-8,
            coordsys="cart",
            trust=0.01,
            tmin=1e-6,
        ),
        recovery_reason=(
            "Restart checked optimizer-returned coordinates with fresh Cartesian "
            "optimizer and tighter SCF; acceptance unchanged"
        ),
        execution_change=(
            "Separate geometry and active SP lifetimes; no scientific model change"
        ),
    )
    atomic_write_json(RECOVERY / "plan.json", plan)
    pointer = ROOT / "results/polymer_models/recovery_campaigns.json"
    registry = json.loads(pointer.read_text())
    registry["campaigns"].append(
        dict(
            path=str(RECOVERY.relative_to(ROOT)),
            plan_sha256=api.digest(RECOVERY / "plan.json"),
        )
    )
    atomic_write_json(pointer, registry)
    print(
        json.dumps(
            dict(checkpoint=point.as_posix(), geometry_QC=sanity["status"], cases=1)
        )
    )


def spacing_runner():
    runner = api.load_module(
        "continuation_spacing_runner", ROOT / "scripts/pipeline/run_polymer_spacing.py"
    )
    runner.SCRIPT = SCRIPT
    calculate = api.calculate_energy

    def evaluate(campaign, task, pilot=False):
        with execution_settings():
            return calculate(campaign, task, pilot=pilot)

    api.calculate_energy = evaluate
    return runner


def recovery_runner():
    runner = api.load_module(
        "continuation_recovery_runner", ROOT / "scripts/pipeline/recover_pi_stacking.py"
    )
    runner.BASE, runner.R.BASE, runner.R.SCRIPT = RECOVERY, RECOVERY, SCRIPT
    calculate = runner.R.calculate

    def evaluate(task, plan, pilot=False):
        with execution_settings(optimizer=True):
            return calculate(task, plan, pilot=pilot)

    runner.R.calculate = evaluate
    return runner.R


def coordinate_recovery(hours):
    from redox.workflows.resources import (
        allocation_guard,
        cpu_allocation,
        host_capacity,
        free_slots,
    )
    import fcntl

    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    with (RECOVERY / "coordinator.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        owner = dict(
            host=socket.gethostname(),
            pid=os.getpid(),
            started=time.time(),
            process_start=Path("/proc/self/stat")
            .read_text()
            .rsplit(")", 1)[1]
            .split()[19],
            deadline=time.time() + hours * 3600,
        )
        attempted = set()
        while time.time() < owner["deadline"]:
            plan = api.read_plan(RECOVERY)
            api.verify_sources(plan["producer_sha256"])
            done = len(list((RECOVERY / "done").glob("*.json")))
            failed = len(list((RECOVERY / "failures").glob("*.json")))
            if done + failed:
                atomic_write_json(
                    RECOVERY / "coordinator.json",
                    dict(
                        owner,
                        status="finite_attempt_terminal",
                        completed=done,
                        failures=failed,
                    ),
                )
                return
            available = free_slots()
            for host in sorted({host for host, _ in available}):
                request = (host, tuple(gpu for name, gpu in available if name == host))
                if request in attempted or (RECOVERY / "claims" / KEY).exists():
                    continue
                command = [
                    sys.executable,
                    str(
                        ROOT.parent
                        / ".codex-personal/skills/lambda-fleet/scripts"
                        / "preflight_node.py"
                    ),
                    host,
                    "--python",
                    sys.executable,
                    "--module",
                    "pyscf",
                    "--module",
                    "gpu4pyscf",
                    "--module",
                    "rdkit",
                    "--module",
                    "ase",
                    "--path",
                    str(RECOVERY / "plan.json"),
                    "--min-available-ram-mib",
                    "32768",
                    "--timeout",
                    "40",
                ]
                probe = subprocess.run(
                    command, capture_output=True, text=True, timeout=50
                )
                (RECOVERY / f"{host}_preflight.json").write_text(probe.stdout)
                if probe.returncode:
                    continue
                # Query hardware without creating a CUDA context or taking a claim.
                probe = subprocess.run(
                    [
                        "ssh",
                        "-n",
                        "-o",
                        "BatchMode=yes",
                        "-o",
                        "ConnectTimeout=8",
                        host,
                        "nvidia-smi --query-gpu=memory.total "
                        "--format=csv,noheader,nounits",
                    ],
                    capture_output=True,
                    text=True,
                    timeout=20,
                )
                if not any(
                    float(line) >= 80 * 1024 for line in probe.stdout.splitlines()
                ):
                    continue
                argv = [
                    "env",
                    "PYTHONPATH=" + str(ROOT / "src"),
                    *[key + "=1" for key in THREAD_KEYS],
                    sys.executable,
                    "-u",
                    str(SCRIPT),
                    "launch",
                    "--role",
                    "queue_large",
                ]
                remote = (
                    "nohup "
                    + shlex.join(argv)
                    + " > "
                    + shlex.quote(str(RECOVERY / f"launch_{host}.log"))
                    + " 2>&1 < /dev/null &"
                )
                subprocess.run(
                    [
                        "ssh",
                        "-n",
                        "-o",
                        "BatchMode=yes",
                        "-o",
                        "ConnectTimeout=8",
                        host,
                        remote,
                    ],
                    check=True,
                    timeout=20,
                )
                attempted.add(request)
            atomic_write_json(
                RECOVERY / "coordinator.json",
                dict(
                    owner,
                    status="active_or_waiting_for_free_suitable_GPU",
                    completed=done,
                    failures=failed,
                ),
            )
            time.sleep(30)
        atomic_write_json(
            RECOVERY / "coordinator.json",
            dict(owner, status="resource_window_ended_outputs_preserved"),
        )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "action",
        choices=(
            "prepare-spacing",
            "prepare-recovery",
            "launch",
            "worker",
            "coordinate",
            "coordinate-recovery",
            "collect",
        ),
    )
    parser.add_argument("--campaign", type=Path)
    parser.add_argument("--stage", choices=("energy", "thermal"), default="energy")
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--role", default="queue_large")
    parser.add_argument("--max-hours", type=float, default=48)
    args = parser.parse_args()
    if args.max_hours <= 0:
        parser.error("max-hours must be positive")
    if args.action == "prepare-spacing":
        prepare_spacing()
    elif args.action == "prepare-recovery":
        prepare_recovery()
    elif args.action == "coordinate-recovery":
        coordinate_recovery(args.max_hours)
    elif args.campaign:
        runner = spacing_runner()
        if args.action == "coordinate":
            runner.coordinate(args.campaign.resolve(), args.max_hours)
        elif args.action == "collect":
            api.collect(args.campaign.resolve())
            runner.render_report(args.campaign.resolve())
        else:
            getattr(runner, args.action)(
                args.campaign.resolve(), args.stage, args.pilot
            )
    else:
        runner = recovery_runner()
        if args.action == "launch":
            from redox.workflows import resources

            original = resources.free_slots
            # One remaining task can use one worker; the limit follows work count.
            resources.free_slots = lambda *a, **kw: original(*a, **kw)[:1]
            try:
                runner.launch(args.role)
            finally:
                resources.free_slots = original
        elif args.action == "worker":
            runner.worker(args.role)
        else:
            parser.error("Specify a spacing --campaign or a supported recovery action")


if __name__ == "__main__":
    main()
