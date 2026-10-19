"""Shared CPU allocation and authoritative GPU availability discovery.

Availability is a snapshot; only gpu_reserve's atomic run establishes ownership.
"""
from __future__ import annotations

import argparse
import csv
from contextlib import contextmanager
import fcntl
import os
from pathlib import Path
import socket
import subprocess

GATE = Path("/nfs/lambda_stor_01/homes/rzhu/bin/gpu_reserve")


@contextmanager
def allocation_guard():
    """Serialize CPU observation and process binding across project launch paths."""
    from redox.core.common import ROOT
    directory = ROOT / "logs/resource_locks"
    directory.mkdir(parents=True, exist_ok=True)
    host = socket.gethostname().split(".")[0]
    with (directory / f"{host}.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        yield


def verify_gate_owner():
    """Reject internal CUDA probes unless their parent owns the system GPU lease."""
    gpu = os.environ.get("CUDA_VISIBLE_DEVICES", "")
    if not gpu.isdigit():
        raise RuntimeError("GPU dependency probe requires a live gpu_reserve run")
    host = socket.gethostname().split(".")[0]
    owner = Path("/nfs/lambda_stor_01/homes/rzhu/.gpu_locks") / f"{host}_gpu{gpu}" / "owner"
    fields = owner.read_text().split()
    parent = os.getppid()
    start = Path(f"/proc/{parent}/stat").read_text().rsplit(")", 1)[1].split()[19]
    if len(fields) != 6 or fields[0] != str(parent) or fields[1] != host or fields[5] != start:
        raise RuntimeError("GPU dependency probe parent does not own the system reservation")


def gpu_suitability(gpu: int):
    """Inspect the required 32-GB-class device without creating a CUDA context."""
    result = subprocess.run(
        ["nvidia-smi", f"--id={gpu}", "--query-gpu=index,uuid,name,memory.total,compute_cap",
         "--format=csv,noheader,nounits"],
        check=True, text=True, capture_output=True, timeout=15,
    )
    rows = list(csv.reader(result.stdout.splitlines(), skipinitialspace=True))
    if len(rows) != 1 or len(rows[0]) != 5:
        raise RuntimeError("Cannot verify GPU identity and memory")
    index, uuid, name, memory, architecture = rows[0]
    if int(index) != gpu or float(memory) < 31 * 1024:
        raise RuntimeError("Full molecular DFT/Hessian requires a validated 32-GB-class GPU")
    return dict(index=gpu, uuid=uuid.strip(), name=name.strip(),
                memory_MiB=float(memory), compute_capability=architecture.strip())


def host_capacity(threads: int, *, minimum_ram_gib=8):
    import psutil
    online = os.cpu_count() or 1
    load = os.getloadavg()[0]
    available = psutil.virtual_memory().available
    if load + threads > online:
        raise RuntimeError("Insufficient shared-node CPU capacity")
    if available < minimum_ram_gib * 1024**3:
        raise RuntimeError("Insufficient available host RAM")
    return dict(online_cpus=online, load_1m=load, available_RAM_bytes=available,
                minimum_RAM_GiB=minimum_ram_gib)


def free_slots(hosts=None):
    argv = [str(GATE), "list"]
    if hosts is not None:
        argv.extend(["--hosts", ",".join(hosts)])
    result = subprocess.run(argv, check=True, text=True, capture_output=True, timeout=60)
    slots = []
    for line in result.stdout.splitlines():
        fields = line.split()
        if len(fields) != 2 or not fields[1].isdigit():
            raise RuntimeError("Unrecognized gpu_reserve availability output")
        slots.append((fields[0], int(fields[1])))
    return slots


def cpu_allocation(threads, requested=None):
    """Choose permitted CPUs outside existing narrow user-process bindings."""
    import psutil
    if threads < 1:
        raise ValueError("threads must be positive")
    permitted = set(os.sched_getaffinity(0))
    busy = set()
    for process in psutil.process_iter():
        if process.pid == os.getpid():
            continue
        try:
            if not process.cmdline():
                continue
            affinity = set(process.cpu_affinity())
        except psutil.NoSuchProcess:
            continue
        except psutil.AccessDenied as exc:
            raise RuntimeError("Cannot verify shared CPU allocations") from exc
        if 0 < len(affinity) <= 8:
            busy.update(affinity)
    selected = sorted(permitted - busy)[:threads] if requested is None else sorted(set(requested))
    if len(selected) != threads or not set(selected) <= permitted or set(selected) & busy:
        raise RuntimeError("Requested CPU allocation is unavailable or overlaps an existing worker")
    return selected


def main(argv=None):
    parser = argparse.ArgumentParser(description="Discover free GPUs through the system reservation policy")
    parser.add_argument("-n", type=int, help="optional user-requested limit")
    parser.add_argument("--hosts", help="comma-separated hosts")
    parser.add_argument("--all", action="store_true", help="query all hosts allowed by the system gate")
    parser.add_argument("--list", action="store_true", help="show free slots grouped by host")
    args = parser.parse_args(argv)
    if args.n is not None and args.n < 1:
        parser.error("-n must be positive")
    if args.all and args.hosts:
        parser.error("choose --all or --hosts")
    local = socket.gethostname().split(".")[0]
    hosts = args.hosts.split(",") if args.hosts else (None if args.all else [local])
    slots = free_slots(hosts)
    if args.n is not None:
        slots = slots[:args.n]
    if args.list:
        for host in sorted({host for host, _ in slots}):
            print(f"{host}: free {[gpu for name, gpu in slots if name == host]}")
    elif not args.all and not args.hosts:
        print(",".join(str(gpu) for _, gpu in slots))
    else:
        for host, gpu in slots:
            print(host, gpu)
    return 0 if slots else 1


if __name__ == "__main__":
    raise SystemExit(main())
