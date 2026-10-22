"""Resumable native-thermochemistry completion for frozen parent benchmarks.

Jobs are molecule-scoped so their native writer locks never compete internally.
Existing partial records are held for inspection. Results retain physical and
spin admission independently of process completion or experimental residuals.
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

from redox.core.common import ROOT, DFT, read_manifest
from redox.core.io import atomic_write_json
from redox.validation.benchmark_panel import digest, validate_panel, verify_hashes
from redox.workflows.checkpoints import verify_state_checkpoint


THREADS = (
    "OMP_NUM_THREADS",
    "OMP_THREAD_LIMIT",
    "MKL_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "XTB_THREADS",
)


def pending_states(panel, reviewed):
    """Select endpoint names from events, including reused TEMPO's actual rad state."""
    needed = {
        (row["id"], state)
        for row in panel["targets"]
        for state in row["event"].split("->")
    }
    return sorted(
        (gid, state)
        for gid, state in needed
        if reviewed.get(gid + ":" + state, {}).get("status") != "ready"
    )


def prepare(campaign, panel_path, review):
    from redox.core.protocol import load_record
    from redox.qm.dft import _opt_protocol_mismatch
    from redox.qm.native_thermal import cache_path

    panel = json.loads(Path(panel_path).read_text())
    validate_panel(panel)
    review = Path(review)
    reviewed = json.loads((review / "states.json").read_text())
    verify_hashes(json.loads((review / "input_inventory.json").read_text()))
    manifest = {(r["id"], r["state"]): r for r in read_manifest()}
    jobs, held = {}, []
    for gid, state in pending_states(panel, reviewed):
        row = manifest[gid, state]
        directory = DFT / gid / state
        raw = json.loads((directory / "result.json").read_text())
        q, mult = int(row["charge"]), int(raw["mult"])
        path, _ = cache_path(directory, q, mult)
        reason = None
        if raw["charge"] != q or _opt_protocol_mismatch(raw, q):
            reason = "source geometry protocol/charge needs review"
        elif load_record(directory, q, mult) is None:
            reason = "active single-point missing"
        elif path.exists() or path.with_suffix("").exists():
            reason = "existing native attempt preserved for inspection"
        item = dict(
            id=gid,
            state=state,
            manifest_row=row,
            input_sha256={
                str(directory / name): digest(directory / name)
                for name in ("opt.xyz", "result.json")
            },
        )
        if reason:
            held.append(dict(item, reason=reason))
        else:
            jobs.setdefault(gid, []).append(item)
    campaign = Path(campaign).resolve()
    campaign.mkdir(parents=True, exist_ok=False)
    for name in ("claims", "done", "tasks", "failures", "workers", "launches"):
        (campaign / name).mkdir()
    sources = [
        *sorted((ROOT / "src/redox/qm").glob("*.py")),
        *sorted((ROOT / "src/redox/core").glob("*.py")),
        *sorted((ROOT / "config").glob("*.py")),
        *sorted((ROOT / "config").glob("*.json")),
        Path(__file__),
        ROOT / "src/redox/validation/benchmark_panel.py",
        ROOT / "src/redox/validation/provenance.py",
        *[
            ROOT / "src/redox/workflows" / (name + ".py")
            for name in ("resources", "gpu", "processes", "checkpoints")
        ],
    ]
    plan = dict(
        schema=1,
        role="exact-parent benchmark method controls; no synthesis admission",
        panel_path=str(Path(panel_path).resolve()),
        panel_sha256=panel["panel_sha256"],
        producer_sha256={str(p): digest(p) for p in sources},
        jobs=jobs,
        held=held,
        created=time.time(),
        stop_new_work_at=time.time() + 48 * 3600,
        maximum_actual_attempts=1,
        physical_QC_unchanged=True,
        experimental_accuracy_claimed=False,
    )
    atomic_write_json(campaign / "plan.json", plan)
    return plan


def complete_state(item):
    from redox.qm.native_thermal import (
        ensure_native_thermal,
        free_energy_native,
        cache_path,
        load_native_thermal,
        select_policy,
    )
    from redox.validation.provenance import referenced_files

    gid, state = item["id"], item["state"]
    directory = DFT / gid / state
    current = [r for r in read_manifest() if r["id"] == gid and r["state"] == state]
    if current != [item["manifest_row"]]:
        raise ValueError("Declared benchmark state changed")
    verify_hashes(item["input_sha256"])
    select_policy(directory)
    with (ROOT / "logs/pipeline_locks" / f"{gid}_native_thermal.lock").open(
        "a+"
    ) as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        try:
            kind, record = ensure_native_thermal(gid, state)
        except RuntimeError as error:
            if not any(s in str(error) for s in ("not stationary", "not a minimum")):
                raise
            from redox.qm.native_thermal_repair import repair_minimum

            kind, record = repair_minimum(gid, state, str(error))
    total, quality = free_energy_native(gid, state)
    raw = json.loads((directory / "result.json").read_text())
    if (
        total is None
        or quality != "ok"
        or load_native_thermal(directory, raw["charge"], raw["mult"]) is None
    ):
        raise ValueError("Native surface/active-energy admission failed")
    path, _ = cache_path(directory, raw["charge"], raw["mult"])
    files = {
        directory / "opt.xyz",
        directory / "result.json",
        directory / "thermal_policy.json",
        path,
    }
    files |= referenced_files(json.loads(path.read_text()), ROOT)
    files |= set(path.with_suffix("").glob("*.npy"))
    files.add(Path(record["thermal_postprocessing_path"]))
    from redox.core.protocol import record_path

    files.add(record_path(directory, raw["charge"], raw["mult"]))
    return dict(
        id=gid,
        state=state,
        G_eV=total,
        kind=kind,
        status="native_surface_verified_spin_admission_separate",
        output_sha256={str(p): digest(p) for p in files},
    )


def worker(campaign):
    from redox.workflows.resources import verify_gate_owner
    from redox.workflows.processes import stop_owned_children

    verify_gate_owner()
    if len(os.sched_getaffinity(0)) != 2 or any(
        os.environ.get(k) != "2" for k in THREADS
    ):
        raise RuntimeError("Reserved two-CPU worker with capped native pools required")
    campaign = Path(campaign).resolve()
    plan = json.loads((campaign / "plan.json").read_text())
    host = socket.gethostname().split(".")[0]
    gpu = os.environ["CUDA_VISIBLE_DEVICES"]
    owner = dict(
        host=host,
        pid=os.getpid(),
        process_start=Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()[19],
        gpu=gpu,
        cpus=sorted(os.sched_getaffinity(0)),
        started=time.time(),
    )
    worker_path = campaign / "workers" / f"{host}_{gpu}_{os.getpid()}.json"
    state = dict(owner, status="running", completed=[])
    atomic_write_json(worker_path, state)

    def stop(sig, frame):
        raise KeyboardInterrupt(f"signal {sig}")

    signal.signal(signal.SIGTERM, stop)
    consecutive = 0
    try:
        for gid, items in plan["jobs"].items():
            if time.time() >= plan["stop_new_work_at"]:
                break
            verify_hashes(plan["producer_sha256"])
            done = campaign / "done" / f"{gid}.json"
            if done.exists():
                verify_state_checkpoint(done)
                continue
            if (campaign / "failures" / f"{gid}.json").exists():
                continue
            claim = campaign / "claims" / gid
            try:
                claim.mkdir()
            except FileExistsError:
                continue
            atomic_write_json(claim / "owner.json", owner)
            report = dict(
                status="computing", task=dict(key=gid), states=[], started=time.time()
            )
            report_path = campaign / "tasks" / f"{gid}.json"
            try:
                # Recheck under the claim before any scientific mutation.
                if done.exists():
                    verify_state_checkpoint(done)
                    continue
                state.update(active_id=gid)
                atomic_write_json(worker_path, state)
                for item in items:
                    atomic_write_json(report_path, report)
                    report["states"].append(complete_state(item))
                    stop_owned_children()
                report.update(
                    status="state_records_verified",
                    finished=time.time(),
                    output_sha256={
                        p: h
                        for r in report["states"]
                        for p, h in r["output_sha256"].items()
                    },
                )
                atomic_write_json(report_path, report)
                atomic_write_json(
                    done,
                    dict(
                        key=gid,
                        status_path=str(report_path),
                        status_sha256=digest(report_path),
                    ),
                )
                verify_state_checkpoint(done)
                state["completed"].append(gid)
                consecutive = 0
            except BaseException as error:
                report.update(
                    status="failed_preserved", error=repr(error), finished=time.time()
                )
                atomic_write_json(report_path, report)
                atomic_write_json(campaign / "failures" / f"{gid}.json", report)
                consecutive += 1
                if not isinstance(error, Exception):
                    raise
            finally:
                stop_owned_children()
                (claim / "owner.json").unlink()
                claim.rmdir()
                atomic_write_json(worker_path, state)
            if consecutive >= 3:
                raise RuntimeError(
                    "Three consecutive actual failures; worker circuit breaker"
                )
        state.update(status="worker_finished", finished=time.time())
    except BaseException as error:
        state.update(
            status="stopped_preserving_outputs", error=repr(error), finished=time.time()
        )
        raise
    finally:
        atomic_write_json(worker_path, state)


def launch(campaign):
    """Preflight and reserve all suitable free slots, up to remaining molecule jobs."""
    from redox.workflows.resources import free_slots

    campaign = Path(campaign).resolve()
    plan = json.loads((campaign / "plan.json").read_text())
    verify_hashes(plan["producer_sha256"])
    remaining = sum(
        not (campaign / "done" / f"{gid}.json").exists()
        and not (campaign / "failures" / f"{gid}.json").exists()
        and not (campaign / "claims" / gid).exists()
        for gid in plan["jobs"]
    )
    slots = free_slots()
    launch_id = str(time.time_ns())
    report = dict(
        status="preflighting",
        remaining_jobs=remaining,
        launched_controllers=[],
        rejected={},
        started=time.time(),
    )
    destination = campaign / "launches" / launch_id
    destination.mkdir()
    atomic_write_json(destination / "status.json", report)
    helper = (
        ROOT.parent / ".codex-personal/skills/lambda-fleet/scripts/preflight_node.py"
    )
    for host in sorted({host for host, _ in slots}):
        if len(report["launched_controllers"]) >= remaining:
            break
        try:
            query = subprocess.run(
                [
                    "ssh",
                    "-n",
                    "-o",
                    "BatchMode=yes",
                    "-o",
                    "ConnectTimeout=8",
                    host,
                    (
                        "nvidia-smi --query-gpu=index,memory.total "
                        "--format=csv,noheader,nounits"
                    ),
                ],
                capture_output=True,
                text=True,
                timeout=20,
                check=True,
            )
            memories = {
                int(line.split(",")[0]): float(line.split(",")[1])
                for line in query.stdout.splitlines()
            }
            suitable = [
                gpu
                for node, gpu in slots
                if node == host and memories.get(gpu, 0) >= 31 * 1024
            ]
            if not suitable:
                report["rejected"][host] = "No free validated 32-GB-class GPU"
                continue
            probe = subprocess.run(
                [
                    sys.executable,
                    str(helper),
                    host,
                    "--python",
                    sys.executable,
                    "--module",
                    "ase",
                    "--module",
                    "rdkit",
                    "--module",
                    "pyscf",
                    "--module",
                    "gpu4pyscf",
                    "--module",
                    "geometric",
                    "--path",
                    str(campaign / "plan.json"),
                    "--min-available-ram-mib",
                    str(8192 * len(suitable)),
                    "--timeout",
                    "40",
                ],
                capture_output=True,
                text=True,
                timeout=50,
            )
            (destination / f"{host}_preflight.json").write_text(probe.stdout)
            if probe.returncode:
                report["rejected"][host] = "Dependency/resource preflight failed"
                continue
            free_now = set(free_slots([host]))
            for gpu in suitable:
                if (host, gpu) not in free_now:
                    continue
                if len(report["launched_controllers"]) >= remaining:
                    break
                run = destination / f"{host}_gpu{gpu}"
                argv = [
                    "env",
                    "PYTHONPATH=" + str(ROOT / "src"),
                    *[name + "=1" for name in THREADS],
                    sys.executable,
                    "-u",
                    "-m",
                    "redox.workflows.gpu",
                    "--gpu",
                    str(gpu),
                    "--out",
                    str(run),
                    "--threads",
                    "2",
                    "--",
                    sys.executable,
                    "-u",
                    "-m",
                    "redox.workflows.benchmark_native",
                    "worker",
                    "--campaign",
                    str(campaign),
                ]
                remote = (
                    "cd "
                    + shlex.quote(str(ROOT))
                    + " && nohup "
                    + shlex.join(argv)
                    + " > "
                    + shlex.quote(str(run) + ".log")
                    + " 2>&1 < /dev/null & echo $!"
                )
                started = subprocess.run(
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
                    text=True,
                    capture_output=True,
                    timeout=20,
                    check=True,
                )
                report["launched_controllers"].append(
                    dict(
                        host=host,
                        gpu=gpu,
                        remote_shell_pid=started.stdout.strip(),
                        argv=argv,
                        execution=str(run / "execution.json"),
                        verified_scientific_starts=False,
                    )
                )
                atomic_write_json(destination / "status.json", report)
        except (OSError, ValueError, subprocess.SubprocessError) as error:
            report["rejected"][host] = repr(error)
        finally:
            atomic_write_json(destination / "status.json", report)
    report.update(
        status="controllers_launched_starts_need_verification", finished=time.time()
    )
    atomic_write_json(destination / "status.json", report)
    print(json.dumps(report, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "worker", "launch"))
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--panel", type=Path)
    parser.add_argument("--review", type=Path)
    args = parser.parse_args()
    if args.action == "prepare":
        if args.panel is None or args.review is None:
            parser.error("prepare requires --panel and --review")
        plan = prepare(args.campaign, args.panel, args.review)
        print(
            json.dumps(
                dict(
                    molecule_jobs=len(plan["jobs"]),
                    states=sum(map(len, plan["jobs"].values())),
                    held=len(plan["held"]),
                )
            )
        )
    elif args.action == "worker":
        worker(args.campaign)
    else:
        launch(args.campaign)


if __name__ == "__main__":
    main()
