#!/usr/bin/env python
"""Run a frozen polymer panel through energy and native-G resource gates.

Preparation is separate. Every GPU process enters through gpu_reserve and the
existing dependency probe. No result is admitted from a successful exit alone.
"""

import argparse
import fcntl
import json
import os
from pathlib import Path
import shlex
import signal
import socket
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as campaign_api
from redox.workflows.checkpoints import verify_state_checkpoint

SCRIPT = Path(__file__).resolve()
THREAD_KEYS = (
    "OMP_NUM_THREADS",
    "OMP_THREAD_LIMIT",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "XTB_THREADS",
)


def render_report(campaign):
    """Refresh presentation separately from numerical admission and producers."""
    subprocess.run(
        [
            str(ROOT / ".venv-figures/bin/python"),
            str(ROOT / "scripts/plotting/polymer_models/polymer_spacing_report.py"),
            "--campaign",
            str(campaign),
        ],
        check=True,
        timeout=60,
    )


def task_input(campaign, task):
    ready = json.loads(
        (campaign / "live_inputs" / task["key"] / "task.json").read_text()
    )
    if any(ready.get(key) != value for key, value in task.items()):
        raise ValueError("Prepared task differs from its frozen identity")
    campaign_api.verify_sources(ready["input_sha256"])
    return ready


def phase_tasks(campaign, plan, stage):
    if stage == "energy":
        return [task_input(campaign, task) for task in plan["tasks"]]
    if not campaign_api.energy_sampling_terminal(campaign, plan):
        return []
    selected, invalid = campaign_api.selected_graph_records(campaign)
    if invalid or len(selected) != len(plan["graph_cases"]):
        return []
    return [record["task"] for record in selected.values()]


def pending_tasks(campaign, plan, stage):
    directory = campaign / stage
    return [
        task
        for task in phase_tasks(campaign, plan, stage)
        if not (directory / "done" / f"{task['key']}.json").exists()
        and not (directory / "failures" / f"{task['key']}.json").exists()
        and not (directory / "claims" / task["key"]).exists()
    ]


def worker(campaign, stage, pilot):
    from redox.workflows.resources import verify_gate_owner
    from redox.workflows.processes import stop_owned_children
    from redox.qm.dft import _free_gpu_pool
    import pyscf.lib

    verify_gate_owner()
    pyscf.lib.num_threads(2)
    if len(os.sched_getaffinity(0)) != 2 or any(
        os.environ.get(key) != "2" for key in THREAD_KEYS
    ):
        raise RuntimeError(
            "Worker CPU allocation or native thread limits differ from the plan"
        )
    plan = campaign_api.read_plan(campaign)
    campaign_api.verify_sources(plan["producer_sha256"])
    campaign_api.verify_sources(plan["input_sha256"])
    owner = dict(
        host=socket.gethostname().split(".")[0],
        pid=os.getpid(),
        process_start=Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()[19],
        gpu=int(os.environ["CUDA_VISIBLE_DEVICES"]),
        cpus=sorted(os.sched_getaffinity(0)),
        stage=stage,
        pilot=pilot,
    )
    role = "pilot" if pilot else "queue"
    report_path = (
        campaign / stage / "workers" / f"{owner['host']}_gpu{owner['gpu']}_{role}.json"
    )
    report = dict(owner, status="starting", items=[], started=time.time())
    atomic_write_json(report_path, report)

    def interrupted(signum, frame):
        raise KeyboardInterrupt(f"signal {signum}")

    signal.signal(signal.SIGTERM, interrupted)
    consecutive_failures = 0
    try:
        if not pilot:
            qualification = json.loads(
                (campaign / stage / "pilot/status.json").read_text()
            )
            if qualification["status"] not in (
                "memory_pilot_passed",
                "thermal_pilot_passed",
            ):
                raise RuntimeError("The phase resource/quality pilot has not passed")
        tasks = phase_tasks(campaign, plan, stage)
        if pilot:
            key = (
                plan["pilot_key"]
                if stage == "energy"
                else max(tasks, key=lambda task: task["active_basis_functions"])["key"]
            )
            tasks = [next(task for task in tasks if task["key"] == key)]
        for task in tasks:
            key = task["key"]
            directory = campaign / stage
            receipt = directory / "done" / f"{key}.json"
            if receipt.exists():
                verify_state_checkpoint(receipt)
                continue
            if (directory / "failures" / f"{key}.json").exists():
                continue
            claim = directory / "claims" / key
            try:
                claim.mkdir()
            except FileExistsError:
                continue
            atomic_write_json(claim / "owner.json", owner)
            try:
                if receipt.exists():
                    verify_state_checkpoint(receipt)
                    continue
                if (directory / "tasks" / key).exists():
                    raise RuntimeError("A prior unreceipted attempt requires review")
                report.update(status="computing", active_key=key)
                atomic_write_json(report_path, report)
                print(
                    f"Actual polymer {stage} "
                    f"{'pilot' if pilot else 'calculation'} start: {key}",
                    flush=True,
                )
                try:
                    if stage == "energy":
                        campaign_api.calculate_energy(campaign, task, pilot=pilot)
                    else:
                        selected, invalid = campaign_api.selected_graph_records(
                            campaign
                        )
                        if invalid:
                            raise RuntimeError(
                                "Energy source records failed verification"
                            )
                        campaign_api.calculate_native_G(
                            campaign, selected[task["model_key"]]
                        )
                    if stage == "energy" and pilot:
                        report["items"].append(
                            dict(key=key, status="resource_pilot_complete")
                        )
                    else:
                        status_path = directory / "tasks" / key / "status.json"
                        atomic_write_json(
                            receipt,
                            dict(
                                key=key,
                                owner=owner,
                                status_path=str(status_path),
                                status_sha256=campaign_api.digest(status_path),
                            ),
                        )
                        verify_state_checkpoint(receipt)
                        report["items"].append(
                            dict(key=key, status="state_records_verified")
                        )
                        if pilot:
                            atomic_write_json(
                                directory / "pilot/status.json",
                                dict(
                                    status="thermal_pilot_passed",
                                    task=task,
                                    reused_result_receipt=str(receipt),
                                    receipt_sha256=campaign_api.digest(receipt),
                                    scope=(
                                        "Native thermochemistry admitted for the "
                                        "selected neutral cohort; conformer ensemble "
                                        "remains unqualified"
                                    ),
                                ),
                            )
                    consecutive_failures = 0
                except Exception as error:
                    status_path = directory / "tasks" / key / "status.json"
                    if status_path.exists():
                        state = json.loads(status_path.read_text())
                        state.update(status="failed_preserved", error=repr(error))
                        atomic_write_json(status_path, state)
                    if pilot:
                        pilot_path = directory / "pilot/status.json"
                        if stage == "thermal" or not pilot_path.exists():
                            atomic_write_json(
                                pilot_path,
                                dict(
                                    status="failed_preserved",
                                    key=key,
                                    error=repr(error),
                                    owner=owner,
                                ),
                            )
                    else:
                        atomic_write_json(
                            directory / "failures" / f"{key}.json",
                            dict(key=key, owner=owner, error=repr(error), attempts=1),
                        )
                    report["items"].append(
                        dict(key=key, status="failed_preserved", error=repr(error))
                    )
                    consecutive_failures += 1
                    if pilot:
                        raise
            finally:
                stop_owned_children()
                _free_gpu_pool()
                (claim / "owner.json").unlink()
                claim.rmdir()
                atomic_write_json(report_path, report)
            campaign_api.collect(campaign)
            if consecutive_failures >= 3:
                raise RuntimeError(
                    "Three consecutive actual failures; preserving remaining work"
                )
        report.update(status="no_unclaimed_unattempted_phase_work", active_key=None)
    except BaseException as error:
        report.update(status="stopped_preserving_outputs", error=repr(error))
        raise
    finally:
        stop_owned_children()
        report["finished"] = time.time()
        atomic_write_json(report_path, report)


def launch(campaign, stage, pilot):
    from redox.workflows.resources import (
        allocation_guard,
        free_slots,
        gpu_suitability,
        cpu_allocation,
        host_capacity,
        GATE,
    )
    import site

    plan = campaign_api.read_plan(campaign)
    campaign_api.verify_sources(plan["producer_sha256"])
    campaign_api.verify_sources(plan["input_sha256"])
    host = socket.gethostname().split(".")[0]
    preflight = json.loads((campaign / f"{host}_preflight.json").read_text())
    if not preflight["ok"] or preflight["python"] != sys.executable:
        raise RuntimeError("A successful node/interpreter preflight is required")
    pending = pending_tasks(campaign, plan, stage)
    if not pending:
        return
    phase = campaign / stage
    role = "pilot" if pilot else "queue"
    record_path = phase / "launches" / f"{host}_{role}.json"
    previous = json.loads(record_path.read_text()) if record_path.exists() else None
    if pilot and (previous or (phase / "pilot/status.json").exists()):
        return
    used = {row["gpu"] for row in previous["workers"]} if previous else set()
    with allocation_guard():
        suitable = []
        for _, gpu in free_slots([host]):
            if gpu in used:
                continue
            try:
                device = gpu_suitability(gpu)
            except RuntimeError:
                continue
            if device["memory_MiB"] >= plan["minimum_large_GPU_MiB"]:
                suitable.append((gpu, device))
        suitable = suitable[: 1 if pilot else len(pending)]
        if not suitable:
            return
        capacity = host_capacity(2 * len(suitable), minimum_ram_gib=32 * len(suitable))
        cpus = cpu_allocation(2 * len(suitable))
        environment = dict(os.environ, PYTHONPATH=str(ROOT / "src"))
        environment.update({key: "2" for key in THREAD_KEYS})
        environment["PATH"] = (
            str(Path(sys.executable).parent) + os.pathsep + environment.get("PATH", "")
        )
        libraries = [
            str(path)
            for base in site.getsitepackages()
            for path in Path(base).glob("nvidia/*/lib")
        ]
        environment["LD_LIBRARY_PATH"] = os.pathsep.join(
            libraries + [environment.get("LD_LIBRARY_PATH", "")]
        )
        record = previous or dict(
            host=host,
            stage=stage,
            role=role,
            workers=[],
            plan_sha256=campaign_api.digest(campaign / "plan.json"),
        )
        record["latest_resources"] = capacity
        for index, (gpu, device) in enumerate(suitable):
            assigned = cpus[2 * index : 2 * index + 2]
            destination = phase / f"reserved_{host}_gpu{gpu}_{role}"
            destination.mkdir()
            command = [
                "taskset",
                "--cpu-list",
                ",".join(map(str, assigned)),
                str(GATE),
                "run",
                str(gpu),
                "--label",
                "redox-polymer-framework",
                "--",
                sys.executable,
                "-m",
                "redox.workflows.gpu",
                "--out",
                str(destination),
                "--probe",
                "--",
                sys.executable,
                "-u",
                str(SCRIPT),
                "worker",
                "--campaign",
                str(campaign),
                "--stage",
                stage,
            ]
            if pilot:
                command.append("--pilot")
            with (phase / f"controller_{host}_gpu{gpu}_{role}.log").open("w") as stream:
                process = subprocess.Popen(
                    command,
                    cwd=ROOT,
                    env=environment,
                    stdin=subprocess.DEVNULL,
                    stdout=stream,
                    stderr=subprocess.STDOUT,
                    start_new_session=True,
                    preexec_fn=lambda assigned=assigned: os.sched_setaffinity(
                        0, assigned
                    ),
                )
            record["workers"].append(
                dict(
                    gpu=gpu,
                    device=device,
                    cpus=assigned,
                    gate_pid=process.pid,
                    argv=command,
                )
            )
            atomic_write_json(record_path, record)
            time.sleep(0.5)
    print(
        f"{host}: submitted {len(suitable)} reserved {stage} {role} workers", flush=True
    )


def owner_alive(owner):
    """Check PID and process start, refusing to guess when a host is unreachable."""
    if owner["host"] == socket.gethostname().split(".")[0]:
        path = Path(f"/proc/{owner['pid']}/stat")
        return path.exists() and (
            path.read_text().rsplit(")", 1)[1].split()[19] == owner["process_start"]
        )
    code = (
        "from pathlib import Path; "
        f"p=Path('/proc/{int(owner['pid'])}/stat'); "
        f"print(p.exists() and p.read_text().rsplit(')',1)[1].split()[19]"
        f"=={owner['process_start']!r})"
    )
    result = subprocess.run(
        [
            "ssh",
            "-n",
            "-o",
            "BatchMode=yes",
            "-o",
            "ConnectTimeout=8",
            owner["host"],
            shlex.join([sys.executable, "-c", code]),
        ],
        capture_output=True,
        text=True,
        timeout=20,
    )
    if result.returncode or result.stdout.strip() not in ("True", "False"):
        raise RuntimeError("Cannot establish previous coordinator ownership")
    return result.stdout.strip() == "True"


def coordinate(campaign, max_hours, resume=False):
    """Finite resource orchestration preserving failed attempts and live leases."""
    from redox.workflows.resources import (
        allocation_guard,
        cpu_allocation,
        host_capacity,
        free_slots,
    )

    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    plan = campaign_api.read_plan(campaign)
    path = campaign / "coordinator.json"
    lock = (campaign / "coordinator.lock").open("a+")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    if path.exists():
        previous = json.loads(path.read_text())
        if not resume or owner_alive(previous):
            raise RuntimeError("Existing coordinator requires review or is still live")
        history = campaign / "coordinator_history"
        history.mkdir(exist_ok=True)
        (history / f"{time.time_ns()}.json").write_bytes(path.read_bytes())
    owner = dict(
        host=socket.gethostname().split(".")[0],
        pid=os.getpid(),
        process_start=Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()[19],
        started=time.time(),
        deadline=time.time() + 3600 * max_hours,
    )
    requested = set()
    last_inventory = None
    current = None
    selected = {}
    invalid = []
    deferred_hosts = {}
    try:
        while time.time() < owner["deadline"]:
            campaign_api.verify_sources(plan["producer_sha256"])
            inventory = tuple(
                sorted(
                    (str(path), path.stat().st_mtime_ns)
                    for stage in ("energy", "thermal")
                    for folder in ("done", "failures")
                    for path in (campaign / stage / folder).glob("*.json")
                )
            )
            if inventory != last_inventory or current is None:
                current = campaign_api.collect(campaign)
                render_report(campaign)
                selected, invalid = campaign_api.selected_graph_records(campaign)
                last_inventory = inventory
            energy_terminal = campaign_api.energy_sampling_terminal(campaign, plan)
            stage = (
                "thermal"
                if energy_terminal and plan["thermochemistry"]["requested"]
                else "energy"
            )
            if energy_terminal and (
                invalid or len(selected) != len(plan["graph_cases"])
            ):
                atomic_write_json(
                    path,
                    dict(
                        owner,
                        status="terminal_missing_qualified_graphs",
                        detail=current["status"],
                    ),
                )
                return
            if energy_terminal and not plan["thermochemistry"]["requested"]:
                atomic_write_json(path, dict(owner, status="electronic_panel_terminal"))
                return
            tasks = (
                plan["tasks"]
                if stage == "energy"
                else [record["task"] for record in selected.values()]
            )
            terminal = tasks and all(
                (campaign / stage / "done" / f"{task['key']}.json").exists()
                or (campaign / stage / "failures" / f"{task['key']}.json").exists()
                for task in tasks
            )
            if stage == "thermal" and terminal:
                atomic_write_json(
                    path,
                    dict(
                        owner, status="finite_panel_terminal", detail=current["status"]
                    ),
                )
                return
            pilot_path = campaign / stage / "pilot/status.json"
            qualification = (
                json.loads(pilot_path.read_text())["status"]
                if pilot_path.exists()
                else None
            )
            if qualification == "failed_preserved":
                atomic_write_json(
                    path, dict(owner, status=f"{stage}_pilot_requires_review")
                )
                return
            pilot = qualification not in ("memory_pilot_passed", "thermal_pilot_passed")
            free = free_slots()
            by_host = {
                host: [gpu for name, gpu in free if name == host] for host, _ in free
            }
            for host, gpus in sorted(by_host.items()):
                if time.time() < deferred_hosts.get(host, 0):
                    continue
                request = (host, stage, pilot, tuple(gpus))
                if request in requested:
                    continue
                device_command = shlex.join(
                    [
                        "nvidia-smi",
                        "--query-gpu=index,memory.total",
                        "--format=csv,noheader,nounits",
                    ]
                )
                try:
                    probe = subprocess.run(
                        [
                            "ssh",
                            "-n",
                            "-o",
                            "BatchMode=yes",
                            "-o",
                            "ConnectTimeout=8",
                            host,
                            device_command,
                        ],
                        capture_output=True,
                        text=True,
                        timeout=25,
                    )
                except (OSError, subprocess.TimeoutExpired):
                    deferred_hosts[host] = time.time() + 300
                    continue
                eligible = [
                    int(line.split(",")[0])
                    for line in probe.stdout.splitlines()
                    if len(line.split(",")) == 2
                    and float(line.split(",")[1]) >= plan["minimum_large_GPU_MiB"]
                ]
                if not set(eligible) & set(gpus):
                    continue
                helper = (
                    ROOT.parent
                    / ".codex-personal/skills/lambda-fleet/scripts/preflight_node.py"
                )
                try:
                    preflight = subprocess.run(
                        [
                            sys.executable,
                            str(helper),
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
                            str(campaign / "plan.json"),
                            "--min-available-ram-mib",
                            "32768",
                            "--timeout",
                            "50",
                        ],
                        capture_output=True,
                        text=True,
                        timeout=55,
                    )
                except (OSError, subprocess.TimeoutExpired):
                    deferred_hosts[host] = time.time() + 300
                    continue
                (campaign / f"{host}_preflight.json").write_text(preflight.stdout)
                if preflight.returncode:
                    deferred_hosts[host] = time.time() + 300
                    continue
                command = [
                    "env",
                    "PYTHONPATH=" + str(ROOT / "src"),
                    *[key + "=1" for key in THREAD_KEYS],
                    sys.executable,
                    "-u",
                    str(SCRIPT),
                    "launch",
                    "--campaign",
                    str(campaign),
                    "--stage",
                    stage,
                ]
                if pilot:
                    command.append("--pilot")
                log = (
                    campaign
                    / stage
                    / f"launch_{host}_{'pilot' if pilot else 'queue'}.log"
                )
                remote = (
                    "nohup "
                    + shlex.join(command)
                    + " > "
                    + shlex.quote(str(log))
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
                requested.add(request)
                if pilot:
                    break
            phase_state = (
                "awaiting_resource_pilot"
                if pilot
                else "active_or_waiting_for_suitable_GPU"
            )
            atomic_write_json(
                path,
                dict(
                    owner,
                    status=f"{stage}_{phase_state}",
                    detail=current["status"],
                ),
            )
            time.sleep(20)
        atomic_write_json(
            path, dict(owner, status="resource_wait_window_ended_outputs_preserved")
        )
    except BaseException as error:
        atomic_write_json(
            path, dict(owner, status="stopped_preserving_outputs", error=repr(error))
        )
        raise
    finally:
        lock.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "action", choices=["freeze", "launch", "worker", "coordinate", "collect"]
    )
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--stage", choices=["energy", "thermal"], default="energy")
    parser.add_argument("--pilot", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--max-hours", type=float, default=48)
    args = parser.parse_args()
    if args.max_hours <= 0:
        parser.error("--max-hours must be positive")
    campaign = args.campaign.resolve()
    if args.action == "freeze":
        campaign_api.freeze(campaign)
        render_report(campaign)
    elif args.action == "collect":
        campaign_api.collect(campaign)
        render_report(campaign)
    elif args.action == "coordinate":
        coordinate(campaign, args.max_hours, args.resume)
    elif args.action == "worker":
        worker(campaign, args.stage, args.pilot)
    else:
        launch(campaign, args.stage, args.pilot)


if __name__ == "__main__":
    main()
