#!/usr/bin/env python
"""Dispatch one frozen oligomer continuation when a suitable GPU becomes free."""

import argparse
import fcntl
import json
import os
from pathlib import Path
import shlex
import socket
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--driver", type=Path, required=True)
    parser.add_argument("--max-hours", type=float, default=24)
    args = parser.parse_args()
    if not 0 < args.max_hours <= 48:
        parser.error("Finite resource window up to 48 hours required")
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity, free_slots

    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    campaign, driver = args.campaign.resolve(), args.driver.resolve()
    plan = api.read_plan(campaign)
    if len(plan["tasks"]) != 1 or plan["producer_sha256"].get(str(driver)) != api.digest(driver):
        raise ValueError("Expected one frozen continuation and its declared driver")
    with (campaign / "coordinator.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        owner = dict(host=socket.gethostname(), pid=os.getpid(),
                     process_start=Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()[19],
                     deadline=time.time() + args.max_hours * 3600)
        attempt = 0
        while time.time() < owner["deadline"]:
            api.verify_sources(plan["producer_sha256"])
            api.verify_sources(plan["input_sha256"])
            if list((campaign / "done").glob("*.json")) or list((campaign / "failures").glob("*.json")):
                atomic_write_json(campaign / "coordinator.json", dict(owner, status="finite_attempt_terminal"))
                return
            if not list((campaign / "claims").iterdir()):
                available = free_slots()
                launched = False
                for host in sorted({host for host, _ in available}):
                    probe = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=8", host,
                                            "nvidia-smi --query-gpu=index,memory.total --format=csv,noheader,nounits"],
                                           capture_output=True, text=True, timeout=20)
                    memory = {int(row.split(",")[0]): float(row.split(",")[1]) for row in probe.stdout.splitlines()} if probe.returncode == 0 else {}
                    eligible = [gpu for name, gpu in available if name == host and memory.get(gpu, 0) >= 81920]
                    if not eligible:
                        continue
                    preflight = [sys.executable, str(ROOT.parent / ".codex-personal/skills/lambda-fleet/scripts/preflight_node.py"),
                                 host, "--python", sys.executable, "--module", "pyscf", "--module", "gpu4pyscf",
                                 "--module", "rdkit", "--module", "ase", "--path", str(campaign / "plan.json"),
                                 "--min-available-ram-mib", "32768", "--timeout", "40"]
                    probe = subprocess.run(preflight, capture_output=True, text=True, timeout=50)
                    (campaign / f"{host}_preflight.json").write_text(probe.stdout)
                    if probe.returncode:
                        continue
                    attempt += 1
                    destination = campaign / "execution" / f"attempt_{attempt}_{time.time_ns()}"
                    caps = ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "XTB_THREADS")
                    argv = ["env", "PYTHONPATH=" + str(ROOT / "src"), *[key + "=1" for key in caps], sys.executable,
                            "-m", "redox.workflows.gpu", "--gpu", str(eligible[0]), "--out", str(destination), "--threads", "2",
                            "--", sys.executable, str(driver), "worker", "--campaign", str(campaign)]
                    log = campaign / f"launch_{attempt}.log"
                    command = "nohup " + shlex.join(argv) + " >" + shlex.quote(str(log)) + " 2>&1 </dev/null & echo $!"
                    launch = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=8", host, command],
                                            capture_output=True, text=True, check=True, timeout=20)
                    atomic_write_json(campaign / f"launch_{attempt}.json", dict(host=host, pid=launch.stdout.strip(),
                                      gpu=eligible[0], argv=argv, log=str(log), execution=str(destination), started=time.time()))
                    # A reservation race may refuse this wrapper. Recheck its outcome
                    # before another attempt; existing claims are never reclaimed here.
                    while time.time() < owner["deadline"]:
                        result = destination / "execution.json"
                        if result.exists() and json.loads(result.read_text()).get("status") in ("failed", "process_completed"):
                            break
                        time.sleep(30)
                    launched = True
                    break
                atomic_write_json(campaign / "coordinator.json", dict(owner, status="worker_terminal_recheck" if launched
                                  else "waiting_all_suitable_GPUs_occupied", attempts=attempt))
            time.sleep(30)
        atomic_write_json(campaign / "coordinator.json", dict(owner, status="finite_resource_window_ended_outputs_preserved"))


if __name__ == "__main__":
    main()
