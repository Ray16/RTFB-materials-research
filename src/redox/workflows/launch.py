"""Plan or run resumable DFT/UMA shards on all suitable free GPUs of this host.

    python -m redox.workflows.launch --module dft --out logs/runs/NAME
    python -m redox.workflows.launch --module dft --out logs/runs/NAME --execute

For multi-node static shards, supply the same --shard-count and disjoint
--shard-start ranges. Each node must pass its own reservation/resource preflight.
"""
from __future__ import annotations

import argparse
import hashlib
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
from redox.workflows.environment import worker_environment
from redox.workflows.resources import free_slots, gpu_suitability, cpu_allocation, host_capacity


def commands(module, gpus, destination, *, threads=2, shard_count=None, shard_start=0,
             extra=(), cpu_sets=None):
    count = shard_count if shard_count is not None else len(gpus)
    if not gpus or threads < 1 or count < 1 or shard_start < 0 or shard_start + len(gpus) > count:
        raise ValueError("Invalid or overlapping shard range")
    jobs = []
    for offset, gpu in enumerate(gpus):
        command = [sys.executable, "-m", f"redox.qm.{module}"]
        command += ["--all", "--backend", "gpu", "--nthreads", str(threads)] if module == "dft" else ["--device", "cuda"]
        command += ["--shard", f"{count}:{shard_start + offset}", *extra]
        assignment = ["--cpus", ",".join(map(str, cpu_sets[offset]))] if cpu_sets else []
        jobs.append([sys.executable, "-m", "redox.workflows.gpu", "--out",
                     str(destination / f"gpu_{gpu}"), "--gpu", str(gpu),
                     "--threads", str(threads), *assignment, "--", *command])
    return jobs


def main(argv=None, *, command_builder=None, report_metadata=None):
    """Run a launch policy without changing module globals."""
    command_builder = commands if command_builder is None else command_builder
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--module", choices=("dft", "uma"), required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--gpus", help="comma-separated local indices; default: all suitable free devices")
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--shard-count", type=int)
    parser.add_argument("--shard-start", type=int, default=0)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("extra", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    if args.threads < 1:
        parser.error("threads must be positive")
    host = socket.gethostname().split(".")[0]
    available = {gpu for name, gpu in free_slots([host]) if name == host}
    gpus = [int(gpu) for gpu in args.gpus.split(",")] if args.gpus else sorted(available)
    if len(set(gpus)) != len(gpus) or not set(gpus) <= available:
        parser.error("GPU selection is duplicated, occupied or reserved")
    suitable = []
    rejected = {}
    for gpu in gpus:
        try:
            suitable.append((gpu, gpu_suitability(gpu)))
        except (RuntimeError, ValueError, OSError, subprocess.SubprocessError) as exc:
            rejected[gpu] = str(exc)
    if args.gpus and rejected:
        parser.error(f"selected GPU failed suitability preflight: {rejected}")
    gpus = [gpu for gpu, _ in suitable]
    if not gpus:
        parser.error("no suitable free 32-GB-class GPUs on this host")
    extra = args.extra[1:] if args.extra[:1] == ["--"] else args.extra
    # A second shard/backend option would invalidate the reviewed assignments.
    owned = ("--shard", "--backend", "--device", "--nthreads")
    if any(value.split("=", 1)[0] in owned for value in extra):
        parser.error("shard/backend/thread assignment belongs to the launcher")
    destination = args.out.resolve()
    try:
        capacity = host_capacity(len(gpus) * args.threads, minimum_ram_gib=8 * len(gpus))
        cpus = cpu_allocation(len(gpus) * args.threads)
        cpu_sets = [cpus[offset:offset + args.threads] for offset in range(0, len(cpus), args.threads)]
        jobs = command_builder(args.module, gpus, destination, threads=args.threads,
                        shard_count=args.shard_count, shard_start=args.shard_start, extra=extra,
                        cpu_sets=cpu_sets)
    except (ValueError, RuntimeError) as exc:
        parser.error(str(exc))
    report = dict(status="planned", host=host, workers_submitted=0, quantum_jobs_started=None,
                  host_preflight=capacity, cpu_sets=cpu_sets,
                  scientifically_verified=False, suitable_gpus=[info for _, info in suitable],
                  rejected_gpus=rejected, argv=jobs, workers=[])
    if report_metadata is not None:
        report.update(report_metadata())
    if not args.execute:
        print(json.dumps(report, indent=2))
        return 0
    destination.mkdir(parents=True, exist_ok=False)
    source = ROOT / f"src/redox/qm/{args.module}.py"
    report["producer_sha256"] = hashlib.sha256(source.read_bytes()).hexdigest()
    status = destination / "launch.json"
    children = []

    def stop(signum, frame):
        raise KeyboardInterrupt(f"signal {signum}")

    signal.signal(signal.SIGTERM, stop)
    try:
        report["status"] = "running"
        atomic_write_json(status, report)
        env = worker_environment(ROOT)
        for gpu, command in zip(gpus, jobs):
            with (destination / f"gpu_{gpu}.log").open("w") as log:
                child = subprocess.Popen(command, cwd=ROOT, env=env, stdout=log,
                                         stderr=subprocess.STDOUT, start_new_session=True)
            children.append(child)
            report["workers"].append(dict(gpu=gpu, controller_pid=child.pid, returncode=None))
            report["workers_submitted"] += 1
            atomic_write_json(status, report)
        for child, worker in zip(children, report["workers"]):
            worker["returncode"] = child.wait()
            atomic_write_json(status, report)
        report["status"] = "processes_completed" if all(child.returncode == 0 for child in children) else "failed"
        return int(report["status"] == "failed")
    except BaseException as exc:
        report.update(status="failed", error=repr(exc))
        for child in children:
            if child.poll() is None:
                os.killpg(child.pid, signal.SIGTERM)
        # Allow the resource controllers to stop their reserved process trees.
        for child in children:
            try:
                child.wait(timeout=45)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
        raise
    finally:
        report["finished"] = time.time()
        atomic_write_json(status, report)


if __name__ == "__main__":
    raise SystemExit(main())
