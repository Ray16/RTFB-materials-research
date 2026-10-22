#!/usr/bin/env python
"""Collect, independently audit and finish the selected-conformer thermal ladder.

The memory repair, SCF repair and missing SP run in their separate frozen
campaigns. This consumer only schedules an extra native G when the newly admitted
alternate conformer changes the selected electronic minimum. It never substitutes
one geometry's thermal correction onto another geometry's electronic energy.
"""

import argparse
import fcntl
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
from redox.workflows.checkpoints import verify_state_checkpoint

SCRIPT = Path(__file__).resolve()
THERMAL = ROOT / "results/polymer_models/campaigns/2026-10-21_spacing_bounded_pcm"
SP = ROOT / "results/polymer_models/campaigns/2026-10-21_spacing_singlepoint_repair"
EXTRA = ROOT / "results/polymer_models/campaigns/2026-10-21_spacing_extra_thermal"
SOURCE = ROOT / "scripts/pipeline/resume_bounded_polymer_thermal.py"
VALIDATOR = ROOT / "scripts/validation/validate_polymer_workflow.py"
KEY = "occ_none__pose1"
THREADS = (
    "OMP_NUM_THREADS",
    "OMP_THREAD_LIMIT",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "XTB_THREADS",
)


def extra_runner():
    module = api.load_module("extra_bounded_thermal", SOURCE)
    module.CAMPAIGN = EXTRA
    runner = module.configure().runner()
    runner.SCRIPT = SCRIPT
    return runner


def prepare_extra():
    if EXTRA.exists():
        return
    parent = api.read_plan(THERMAL)
    receipt = SP / "energy/done" / f"{KEY}.json"
    verify_state_checkpoint(receipt)
    EXTRA.mkdir()
    for stage in ("energy", "thermal"):
        for folder in (
            "done",
            "tasks",
            "claims",
            "failures",
            "workers",
            "pilot",
            "launches",
        ):
            (EXTRA / stage / folder).mkdir(parents=True)
    shutil.copyfile(receipt, EXTRA / "energy/done" / receipt.name)
    shutil.copyfile(
        THERMAL / "energy/pilot/status.json", EXTRA / "energy/pilot/status.json"
    )
    for name in ("cases.json", "reaction_terms.json", "replay_qualification.json"):
        shutil.copyfile(THERMAL / name, EXTRA / name)
    shutil.copytree(THERMAL / "controls", EXTRA / "controls")
    ready = EXTRA / "live_inputs" / KEY
    ready.mkdir(parents=True)
    shutil.copyfile(THERMAL / "live_inputs" / KEY / "task.json", ready / "task.json")
    graph = next(row for row in parent["graph_cases"] if row["key"] == "occ_none")
    task = next(row for row in parent["tasks"] if row["key"] == KEY)
    plan = dict(
        parent,
        tasks=[task],
        graph_cases=[graph],
        cycles=[],
        raw_root=str(ROOT / "calcs/polymer_models" / EXTRA.name),
        property_scope=(
            "Native G for a newly selected, replay-checked alternate conformer; "
            "no mixed geometries"
        ),
    )
    plan["producer_sha256"].update(
        {str(SCRIPT): api.digest(SCRIPT), str(SOURCE): api.digest(SOURCE)}
    )
    inputs = [receipt, THERMAL / "plan.json"] + list(EXTRA.rglob("*.json"))
    plan["input_sha256"] = {str(path): api.digest(path) for path in inputs}
    atomic_write_json(EXTRA / "plan.json", plan)
    api.collect(EXTRA)


def launch_extra():
    from redox.workflows.resources import free_slots

    for host in sorted({host for host, _ in free_slots()}):
        if (EXTRA / "thermal/launches" / f"{host}_pilot.json").exists():
            continue
        query = subprocess.run(
            [
                "ssh",
                "-n",
                "-o",
                "BatchMode=yes",
                "-o",
                "ConnectTimeout=8",
                host,
                "nvidia-smi --query-gpu=memory.total --format=csv,noheader,nounits",
            ],
            capture_output=True,
            text=True,
            timeout=20,
        )
        if not any(float(line) >= 95830 for line in query.stdout.splitlines()):
            continue
        helper = (
            ROOT.parent
            / ".codex-personal/skills/lambda-fleet/scripts/preflight_node.py"
        )
        probe = subprocess.run(
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
                str(EXTRA / "plan.json"),
                "--min-available-ram-mib",
                "131072",
                "--timeout",
                "40",
            ],
            capture_output=True,
            text=True,
            timeout=50,
        )
        (EXTRA / f"{host}_preflight.json").write_text(probe.stdout)
        if probe.returncode:
            continue
        argv = [
            "env",
            "PYTHONPATH=" + str(ROOT / "src"),
            *[key + "=1" for key in THREADS],
            sys.executable,
            "-u",
            str(SCRIPT),
            "launch",
            "--campaign",
            str(EXTRA),
            "--stage",
            "thermal",
            "--pilot",
        ]
        remote = (
            "nohup "
            + shlex.join(argv)
            + " > "
            + shlex.quote(str(EXTRA / f"launch_{host}.log"))
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
        return


def watch(hours):
    from redox.workflows.resources import (
        allocation_guard,
        cpu_allocation,
        host_capacity,
    )

    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    directory = ROOT / "results/polymer_models/completion"
    directory.mkdir(exist_ok=True)
    with (directory / "owner.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        owner = dict(
            host=socket.gethostname(),
            pid=os.getpid(),
            process_start=Path("/proc/self/stat")
            .read_text()
            .rsplit(")", 1)[1]
            .split()[19],
            deadline=time.time() + hours * 3600,
        )
        validator = api.load_module("independent_polymer_completion", VALIDATOR)
        previous = None
        try:
            while time.time() < owner["deadline"]:
                roots = [THERMAL / "thermal", SP / "energy", EXTRA / "thermal"]
                inventory = tuple(
                    sorted(
                        (str(path), path.stat().st_mtime_ns)
                        for root in roots
                        for kind in ("done", "failures")
                        for path in (root / kind).glob("*.json")
                    )
                )
                original_status = (
                    ROOT / "results/polymer_models/interaction_status.json"
                )
                inventory += (
                    (str(original_status), original_status.stat().st_mtime_ns),
                )
                watched = [
                    THERMAL / "thermal/pilot/status.json",
                    SP / "control.json",
                    SP / "target.json",
                ]
                inventory += tuple(
                    (str(path), path.stat().st_mtime_ns)
                    for path in watched
                    if path.exists()
                )
                if inventory != previous:
                    result = validator.validate(
                        THERMAL,
                        ROOT / "results/polymer_models/workflow_validation.json",
                    )
                    previous = inventory
                    if result["status"] == "finite_workflow_validated":
                        atomic_write_json(
                            directory / "status.json",
                            dict(owner, status="finite_workflow_validated"),
                        )
                        return
                else:
                    result = json.loads(
                        (
                            ROOT / "results/polymer_models/workflow_validation.json"
                        ).read_text()
                    )
                initial_G = len(list((THERMAL / "thermal/done").glob("*.json")))
                if initial_G == 6 and (SP / "energy/done" / f"{KEY}.json").exists():
                    selected = next(
                        row for row in result["cases"] if row["key"] == "occ_none"
                    )
                    if (
                        selected["selected_conformer"] == KEY
                        and selected["G_eV"] is None
                    ):
                        prepare_extra()
                        if (
                            not (EXTRA / "thermal/claims" / KEY).exists()
                            and not (EXTRA / "thermal/pilot/status.json").exists()
                        ):
                            launch_extra()
                failures = []
                for path in watched:
                    if path.exists():
                        state = json.loads(path.read_text())
                        if state.get("status") == "failed_preserved":
                            failures.append(
                                dict(path=str(path), error=state.get("error"))
                            )
                atomic_write_json(
                    directory / "status.json",
                    dict(
                        owner,
                        status="needs_review_with_other_work_preserved"
                        if failures
                        else "waiting_for_qualified_outputs",
                        current_failures=failures,
                        validation_status=result["status"],
                        original_checked=result["original_checked"],
                        spacing_conformers=result["qualified_spacing_conformers"],
                        native_G=result["qualified_native_G"],
                    ),
                )
                time.sleep(30)
            atomic_write_json(
                directory / "status.json",
                dict(owner, status="bounded_completion_window_ended_outputs_preserved"),
            )
        except BaseException as error:
            atomic_write_json(
                directory / "status.json",
                dict(
                    owner, status="needs_review_preserving_outputs", error=repr(error)
                ),
            )
            raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "action", choices=("watch", "launch", "worker"), nargs="?", default="watch"
    )
    parser.add_argument("--max-hours", type=float, default=48)
    parser.add_argument("--campaign", type=Path, default=EXTRA)
    parser.add_argument("--stage", default="thermal", choices=("thermal",))
    parser.add_argument("--pilot", action="store_true")
    args = parser.parse_args()
    if args.max_hours <= 0 or args.campaign.resolve() != EXTRA:
        parser.error("Use the frozen extra campaign and a positive time window")
    if args.action == "watch":
        watch(args.max_hours)
    else:
        getattr(extra_runner(), args.action)(EXTRA, args.stage, args.pilot)
