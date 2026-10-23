#!/usr/bin/env python
"""Finite, source-bound continuation of the all-pose native thermal campaign.

The initial reserved pilot must already have been submitted. Queue expansion
requires its checked result; existing task failures are not retried or removed.
Worker/GPU/RAM suitability remains enforced by the campaign launcher.
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
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--host", required=True)
    parser.add_argument("--max-hours", type=float, default=24)
    args = parser.parse_args()
    if not 0 < args.max_hours <= 48 or args.host.startswith("-"):
        parser.error("A valid host and finite duration up to 48 hours required")
    campaign = args.campaign.resolve()
    script = ROOT / "scripts/pipeline/polymer_conformer_thermals.py"
    source = {str(path): api.digest(path) for path in (Path(__file__), script, campaign / "plan.json")}
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    with (campaign / "coordinator.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        owner_path = campaign / "coordinator.json"
        if owner_path.exists():
            raise RuntimeError("Prior coordinator must be reviewed; no blind restart")
        owner = dict(
            host=socket.gethostname(), pid=os.getpid(),
            process_start=Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()[19],
            CPU_affinity=sorted(os.sched_getaffinity(0)), started=time.time(),
            deadline=time.time() + args.max_hours * 3600, source_sha256=source,
            status="waiting_for_submitted_native_pilot", requested_worker_host=args.host,
            scope="Finite two-pose control only; broader conformer convergence remains open",
        )
        atomic_write_json(owner_path, owner)
        def stop(signum, frame):
            raise KeyboardInterrupt(f"signal {signum}; quantum workers preserved")
        signal.signal(signal.SIGTERM, stop)
        previous_inventory = None
        try:
            while time.time() < owner["deadline"]:
                api.verify_sources(source)
                plan = api.read_plan(campaign)
                api.verify_sources(plan["producer_sha256"])
                api.verify_sources(plan["input_sha256"])
                pilot_path = campaign / "thermal/pilot/status.json"
                pilot = json.loads(pilot_path.read_text()) if pilot_path.exists() else {}
                if pilot.get("status") == "failed_preserved":
                    owner.update(status="pilot_failed_outputs_preserved", pilot=pilot)
                    return
                if pilot.get("status") == "thermal_pilot_passed":
                    argv = [
                        "env", "PYTHONPATH=" + str(ROOT / "src"),
                        *[name + "=1" for name in ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "XTB_THREADS")],
                        sys.executable, str(script), "launch", "--campaign", str(campaign),
                    ]
                    with (campaign / "coordinator_launch.log").open("a") as log:
                        # A timeout of this observation is not an input failure or
                        # permission to resubmit. Inspect the saved launch record.
                        proc = subprocess.run(
                            ["ssh", "-n", "-o", "BatchMode=yes", "-o", "ConnectTimeout=8", args.host, shlex.join(argv)],
                            stdout=log, stderr=subprocess.STDOUT, timeout=55,
                        )
                    if proc.returncode:
                        raise RuntimeError("Queue launch needs review; no calculation failure inferred")
                    owner["status"] = "native_queue_running_or_resource_pending"
                paths = sorted((campaign / "thermal/done").glob("*.json"))
                failures = sorted((campaign / "thermal/failures").glob("*.json"))
                inventory = tuple((str(path), api.digest(path)) for path in paths + failures)
                if inventory != previous_inventory:
                    collect = api.load_module("finite_ensemble_collector", script).collect
                    current = collect(campaign)
                    previous_inventory = inventory
                    owner.update(checked_native_G=current["checked_native_G"], expected_DFT_poses=current["expected_DFT_poses"], failures=len(failures))
                    if current["status"] == "finite_two_pose_partitions_complete":
                        owner["status"] = "finite_two_pose_partitions_complete"
                        return
                if len(paths) + len(failures) == len(plan["tasks"]):
                    owner["status"] = "finite_native_panel_terminal_with_explicit_failures"
                    return
                owner["updated"] = time.time()
                atomic_write_json(owner_path, owner)
                time.sleep(30)
            owner["status"] = "finite_window_ended_outputs_preserved"
        except BaseException as error:
            owner.update(status="observation_or_workflow_needs_review_outputs_preserved", error=repr(error))
            raise
        finally:
            owner["finished"] = time.time()
            atomic_write_json(owner_path, owner)


if __name__ == "__main__":
    main()
