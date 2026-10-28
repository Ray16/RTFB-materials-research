"""Run one full molecular QM command with resource preflight and a GPU reservation.

    python -m redox.workflows.gpu --gpu 0 --out logs/runs/NAME -- \
        python -m redox.qm.sp --only ID --backend gpu

The output directory must be new. A completed process is execution evidence;
scientific outputs still need their domain's readiness and physical QC checks.
"""
from __future__ import annotations

import argparse
import importlib
import os
from pathlib import Path
import signal
import shutil
import socket
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows.environment import worker_environment
from redox.workflows.resources import (GATE, allocation_guard, cpu_allocation, free_slots,
                                       gpu_suitability, host_capacity, verify_gate_owner)


THREAD_KEYS = ("OMP_NUM_THREADS", "OMP_THREAD_LIMIT", "MKL_NUM_THREADS",
               "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "XTB_THREADS")


def _probe_and_exec(destination, command):
    """GPU library imports happen only inside the live system reservation."""
    verify_gate_owner()
    from redox.qm.native_thermal import installed_versions, NATIVE_THERMAL
    versions = installed_versions()
    if versions != {"pyscf": NATIVE_THERMAL["validated_pyscf_version"],
                    "gpu4pyscf": NATIVE_THERMAL["validated_gpu4pyscf_version"]}:
        raise RuntimeError("GPU QM libraries differ from the validated native implementation")
    for module in ("ase", "rdkit", "pyscf", "geometric", "gpu4pyscf", "cupy"):
        importlib.import_module(module)
    import cupy
    if cupy.cuda.runtime.getDeviceCount() != 1:
        raise RuntimeError("Expected exactly one GPU visible through the reservation")
    atomic_write_json(destination / "dependency_preflight.json", dict(status="passed", software=versions))
    os.execvpe(command[0], command, os.environ)


def run_gpu(destination, gpu, command, *, threads=2, requested_cpus=None):
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=False)
    host = socket.gethostname().split(".")[0]
    report = dict(status="preparing", host=host, gpu=gpu, pid=os.getpid(),
                  process_start=Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()[19],
                  interpreter=sys.executable, argv=command, started=time.time(),
                  scientific_outputs_verified=False)
    status = destination / "execution.json"
    atomic_write_json(status, report)
    child = None
    try:
        with allocation_guard():
            # Serialize observation/allocation/start on this host. Other projects
            # are also respected by inspecting existing narrow process bindings.
            cpus = cpu_allocation(threads, requested_cpus)
            report["host_preflight"] = host_capacity(threads)
            report["gpu_preflight"] = gpu_suitability(gpu)
            if (host, gpu) not in free_slots([host]):
                raise RuntimeError("GPU is occupied or reserved; no calculation started")
            env = dict(os.environ, PYTHONPATH=str(ROOT / "src"))
            env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + env.get("PATH", "")
            import site
            libraries = [str(path) for base in site.getsitepackages()
                         for path in Path(base).glob("nvidia/*/lib")]
            env["LD_LIBRARY_PATH"] = os.pathsep.join(libraries + [env.get("LD_LIBRARY_PATH", "")])
            for key in THREAD_KEYS:
                env[key] = str(threads)
            taskset = shutil.which("taskset")
            if not taskset:
                raise RuntimeError("taskset is required for the GPU worker tree")
            argv = [taskset, "--cpu-list", ",".join(map(str, cpus)),
                    str(GATE), "run", str(gpu), "--label", "redox", "--", sys.executable,
                    "-m", "redox.workflows.gpu", "--out", str(destination), "--probe", "--", *command]
            report.update(cpu_affinity=cpus, threads=threads, gate_argv=argv,
                          status="starting_reserved_command")
            atomic_write_json(status, report)
            with (destination / "command.log").open("w") as log:
                # This CLI is single-threaded. Set affinity before exec so the
                # next allocator observes the binding before our lock releases.
                child = subprocess.Popen(argv, cwd=ROOT, env=env, stdout=log,
                                         stderr=subprocess.STDOUT, start_new_session=True,
                                         preexec_fn=lambda: os.sched_setaffinity(0, cpus))
            report["gate_pid"] = child.pid
        returncode = child.wait()
        report.update(returncode=returncode, status="process_completed" if returncode == 0 else "failed")
        return returncode
    except BaseException as exc:
        report.update(status="failed", error=repr(exc))
        if child is not None and child.poll() is None:
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.wait(timeout=30)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
        raise
    finally:
        report["finished"] = time.time()
        atomic_write_json(status, report)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--gpu", type=int)
    parser.add_argument("--threads", type=int, default=2)
    parser.add_argument("--cpus", help="comma-separated disjoint CPU indices; default: allocate")
    parser.add_argument("--probe", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("command", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    command = args.command[1:] if args.command[:1] == ["--"] else args.command
    if not command or args.threads < 1:
        parser.error("a command and positive thread count are required")
    if args.probe:
        _probe_and_exec(args.out, command)
    if args.gpu is None or args.gpu < 0:
        parser.error("--gpu must name a local GPU index")
    cpus = [int(cpu) for cpu in args.cpus.split(",")] if args.cpus else None

    def stop(signum, frame):
        raise KeyboardInterrupt(f"signal {signum}")

    signal.signal(signal.SIGTERM, stop)
    return run_gpu(args.out, args.gpu, command, threads=args.threads, requested_cpus=cpus)


if __name__ == "__main__":
    raise SystemExit(main())
