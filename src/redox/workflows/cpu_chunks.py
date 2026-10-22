"""Launch a measured CPU chunk experiment with explicit per-host concurrency.

This is a CPU-only executor. Quantum stages retain the shared GPU reservation
workflow; a cheap descriptor completion cannot establish quantum readiness.
"""

import argparse
import json
import os
from pathlib import Path
import shlex
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import chunked


def launch(campaign, hosts, out):
    plan = chunked.read_plan(campaign, require_current=True)
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    helper = (
        ROOT.parent / ".codex-personal/skills/lambda-fleet/scripts/preflight_node.py"
    )
    report = dict(
        status="preflighting",
        campaign=str(Path(campaign).resolve()),
        manifest_sha256=plan["manifest_sha256"],
        executor_sha256=chunked.file_hash(__file__),
        requested_hosts=hosts,
        workers=[],
        started=time.time(),
        quantum_jobs_submitted=0,
        launcher_pid=os.getpid(),
    )
    atomic_write_json(out / "launch.json", report)
    for host, count in hosts.items():
        if count < 1 or not host.replace("-", "").replace(".", "").isalnum():
            raise ValueError(
                "Positive explicit CPU concurrency and valid host required"
            )
        probe = subprocess.run(
            [
                sys.executable,
                str(helper),
                host,
                "--python",
                sys.executable,
                "--module",
                "rdkit",
                "--module",
                "numpy",
                "--module",
                "scipy",
                "--path",
                str(Path(campaign).resolve() / "plan.json"),
                "--min-available-ram-mib",
                str(1024 * count),
                "--timeout",
                "40",
            ],
            text=True,
            capture_output=True,
            timeout=50,
        )
        (out / f"{host}_preflight.json").write_text(probe.stdout)
        if probe.returncode:
            raise RuntimeError(f"CPU stage preflight failed: {host}")
        # Verify aggregate CPU capacity; each worker then atomically chooses a
        # disjoint CPU under the existing allocation guard immediately at start.
        check = (
            "import os; from redox.workflows.resources import host_capacity; "
            f"host_capacity({count},minimum_ram_gib={count}); "
            "print('capacity verified')"
        )
        probe = subprocess.run(
            [
                "ssh",
                "-n",
                "-o",
                "BatchMode=yes",
                "-o",
                "ConnectTimeout=8",
                host,
                shlex.join(
                    [
                        "env",
                        "PYTHONPATH=" + str(ROOT / "src"),
                        "OPENBLAS_NUM_THREADS=1",
                        sys.executable,
                        "-c",
                        check,
                    ]
                ),
            ],
            capture_output=True,
            text=True,
            timeout=15,
        )
        if probe.returncode:
            raise RuntimeError(f"CPU stage aggregate capacity failed: {host}")
        for index in range(count):
            log = out / f"{host}_worker{index}.log"
            argv = [
                "env",
                "PYTHONPATH=" + str(ROOT / "src"),
                *[
                    key + "=1"
                    for key in (
                        "OMP_NUM_THREADS",
                        "OMP_THREAD_LIMIT",
                        "MKL_NUM_THREADS",
                        "OPENBLAS_NUM_THREADS",
                        "NUMEXPR_NUM_THREADS",
                        "XTB_THREADS",
                    )
                ],
                "timeout",
                "--signal=TERM",
                "--kill-after=30",
                "1h",
                sys.executable,
                "-u",
                "-m",
                "redox.screening.streaming_search",
                "worker",
                "--campaign",
                str(Path(campaign).resolve()),
                "--checkpoint-every",
                str(plan["chunk_size"]),
            ]
            remote = (
                "nohup "
                + shlex.join(argv)
                + " > "
                + shlex.quote(str(log))
                + " 2>&1 < /dev/null & echo $!"
            )
            result = subprocess.run(
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
                capture_output=True,
                text=True,
                timeout=15,
                check=True,
            )
            report["workers"].append(
                dict(
                    host=host,
                    remote_shell_pid=result.stdout.strip(),
                    argv=argv,
                    log=str(log),
                )
            )
            atomic_write_json(out / "launch.json", report)
    report.update(
        status="controllers_launched_actual_outputs_need_verification",
        finished=time.time(),
    )
    atomic_write_json(out / "launch.json", report)
    print(json.dumps(dict(workers=len(report["workers"]), out=str(out))))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--hosts",
        required=True,
        help="Explicit benchmark concurrency, e.g. lambda1:2,lambda2:2",
    )
    args = parser.parse_args()
    pairs = [value.split(":") for value in args.hosts.split(",")]
    if len({p[0] for p in pairs}) != len(pairs):
        parser.error("Duplicate host")
    try:
        launch(args.campaign, {host: int(count) for host, count in pairs}, args.out)
    except BaseException as error:
        path = args.out.resolve() / "launch.json"
        if path.exists():
            report = json.loads(path.read_text())
            if report.get("launcher_pid") == os.getpid():
                report.update(
                    status="launch_failed_preserving_existing_workers",
                    error=repr(error),
                )
                atomic_write_json(path, report)
        raise


if __name__ == "__main__":
    main()
