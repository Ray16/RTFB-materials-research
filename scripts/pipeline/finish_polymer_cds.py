#!/usr/bin/env python
"""Coordinate the qualified CDS continuation and refresh independent polymer QC."""

import argparse
import fcntl
import json
import os
from pathlib import Path
import signal
import socket
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

CAMPAIGN = ROOT / "results/polymer_models/campaigns/2026-10-21_spacing_cds_recovery"
ORIGINAL = ROOT / "results/polymer_models/campaigns/2026-10-21_spacing_bounded_pcm"
DRIVER = ROOT / "scripts/pipeline/resume_polymer_cds.py"
VALIDATOR = ROOT / "scripts/validation/validate_polymer_workflow.py"


def inventory():
    paths = list(CAMPAIGN.glob("thermal/done/*.json")) + list(CAMPAIGN.glob("thermal/failures/*.json"))
    paths += [ROOT / "results/polymer_models/interaction_status.json"]
    return tuple((str(p), p.stat().st_mtime_ns) for p in sorted(paths) if p.exists())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--max-hours", type=float, default=48)
    args = parser.parse_args()
    if not 0 < args.max_hours <= 48:
        parser.error("Use a positive finite window up to 48 hours")
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity

    permitted = sorted(os.sched_getaffinity(0))
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    directory = ROOT / "results/polymer_models/completion"
    with (directory / "owner.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        owner = dict(host=socket.gethostname(), pid=os.getpid(),
                     process_start=Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()[19],
                     deadline=time.time() + 3600 * args.max_hours, controller_source_sha256=api.digest(__file__))
        child = None
        previous = None
        validator = api.load_module("CDS_independent_review", VALIDATOR)

        def stop(signum, frame):
            raise KeyboardInterrupt(f"signal {signum}")

        signal.signal(signal.SIGTERM, stop)
        try:
            while time.time() < owner["deadline"]:
                plan = api.read_plan(CAMPAIGN)
                api.verify_sources(plan["producer_sha256"])
                api.verify_sources(plan["input_sha256"])
                pilot_path = CAMPAIGN / "thermal/pilot/status.json"
                pilot = json.loads(pilot_path.read_text()) if pilot_path.exists() else {}
                if pilot.get("status") == "thermal_pilot_passed" and child is None:
                    with (directory / "CDS_thermal_coordinator.log").open("a") as log:
                        argv = [sys.executable, str(DRIVER), "coordinate", "--max-hours",
                                str(max(.001, (owner["deadline"] - time.time()) / 3600))]
                        child = subprocess.Popen(argv, stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
                                                 preexec_fn=lambda: os.sched_setaffinity(0, permitted))
                    atomic_write_json(directory / "CDS_coordinator_launch.json",
                                      dict(pid=child.pid, argv=argv, host=socket.gethostname(),
                                           process_start=Path(f"/proc/{child.pid}/stat").read_text().rsplit(")", 1)[1].split()[19]))
                current = inventory()
                if current != previous:
                    validator.validate(ORIGINAL, ROOT / "results/polymer_models/workflow_validation.json")
                    previous = current
                result = json.loads((ROOT / "results/polymer_models/workflow_validation.json").read_text())
                status = dict(owner, status="finite_workflow_validated" if result["status"] == "finite_workflow_validated"
                              else "calculations_and_physical_review_in_progress", validation_status=result["status"],
                              original_checked=result["original_checked"], spacing_conformers=result["qualified_spacing_conformers"],
                              native_G=result["qualified_native_G"], pilot_status=pilot.get("status", "pending"),
                              thermal_coordinator_returncode=child.poll() if child else None, updated=time.time())
                atomic_write_json(directory / "status.json", status)
                if status["status"] == "finite_workflow_validated":
                    return
                if pilot.get("status") == "failed_preserved":
                    raise RuntimeError("CDS pilot failed its unchanged physical QC; outputs preserved")
                if child is not None and child.poll() not in (None, 0):
                    raise RuntimeError("Thermal coordinator failed; inspect its log")
                time.sleep(30)
            atomic_write_json(directory / "status.json", dict(owner, status="finite_resource_window_ended_outputs_preserved"))
        except BaseException as error:
            atomic_write_json(directory / "status.json", dict(owner, status="needs_review_preserving_outputs", error=repr(error)))
            raise
        finally:
            if child is not None and child.poll() is None:
                child.terminate()
                child.wait(timeout=40)


if __name__ == "__main__":
    main()
