"""Independently recompute selected historical SP adoptions on frozen geometries.

Preparation and collection are CPU-only. Workers run inside workflows.gpu and
claim a shared resumable queue. Collection can replace an agreeing historical
record only after preserving its exact original bytes and rechecking all inputs.
Changed energies remain explicit review items rather than silently replacing data.
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import gc
import hashlib
import json
import math
from importlib import metadata
import os
from pathlib import Path
import shutil
import signal
import socket
import time
import traceback

from redox.core.common import ROOT, DFT, HARTREE_EV, read_manifest
from redox.core.io import atomic_text_writer, atomic_write_json
from redox.core.protocol import ACTIVE_SP, geom_sha1, load_record, record_path, software_versions

AGREEMENT_EV = .001


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def read_json(path, default=None):
    try:
        return json.loads(Path(path).read_text())
    except FileNotFoundError:
        return default


def solution_diagnostics_ok(record):
    keys = ("smd_homo_eV", "smd_s_squared", "smd_spin_contam")
    return (all(type(record.get(key)) in (int, float) and math.isfinite(record[key]) for key in keys)
            and abs(record["smd_spin_contam"]) <= .1)


def prepare(output, rows=None, root=DFT, selection="adopted"):
    if selection not in ("adopted", "missing-smd-diagnostics"):
        raise ValueError("Unknown revalidation selection")
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=False)
    for name in ("inputs", "claims", "done", "failures", "workers"):
        (output / name).mkdir()
    tasks = []
    for row in read_manifest() if rows is None else rows:
        sd = Path(root) / row["id"] / row["state"]
        raw = read_json(sd / "result.json")
        if raw is None:
            continue
        old = load_record(sd, raw["charge"], raw["mult"])
        if old is None:
            continue
        if selection == "adopted" and not str(old.get("source", "")).startswith("adopted"):
            continue
        if selection == "missing-smd-diagnostics" and solution_diagnostics_ok(old):
            continue
        key = f"{row['id']}__{row['state']}"
        frozen = output / "inputs" / key
        frozen.mkdir()
        source_record = record_path(sd, raw["charge"], raw["mult"])
        inputs = [sd / "result.json", sd / "opt.xyz", source_record]
        before = {str(path.resolve()): sha(path) for path in inputs}
        for path, name in ((inputs[0], "result.json"), (inputs[1], "opt.xyz"),
                           (source_record, "previous_record.json")):
            shutil.copy2(path, frozen / name)
            if sha(frozen / name) != before[str(path.resolve())]:
                raise RuntimeError(f"Input changed while freezing {key}")
        if any(sha(path) != digest for path, digest in before.items()):
            raise RuntimeError(f"Input changed while freezing {key}")
        tasks.append(dict(key=key, id=row["id"], state=row["state"], charge=raw["charge"],
                          mult=raw["mult"], frozen=str(frozen), state_dir=str(sd.resolve()),
                          require_solution_diagnostics=selection == "missing-smd-diagnostics",
                          old_record=str(source_record.resolve()), inputs_sha256=before,
                          expected_protocol_hash=source_record.stem))
    source_paths = [Path(__file__), ROOT / "src/redox/qm/dft.py", ROOT / "src/redox/qm/sp.py",
                    ROOT / "src/redox/core/common.py", ROOT / "src/redox/core/protocol.py"]
    manifest = dict(schema=1, status="prepared_not_submitted", protocol=ACTIVE_SP, selection=selection,
                    agreement_eV=AGREEMENT_EV, max_attempts=2, tasks=tasks,
                    scientific_accuracy_validated=False, quantum_jobs_submitted=0,
                    code_sha256={str(path): sha(path) for path in source_paths}, created=time.time())
    atomic_write_json(output / "queue.json", manifest)
    return manifest


def verified_candidate(task):
    frozen = Path(task["frozen"])
    record = load_record(frozen, task["charge"], task["mult"])
    if record is None or record.get("source") not in ("computed:gpu:independent_revalidation",
                                                     "computed:cpu:independent_revalidation"):
        return None
    if record.get("protocol_hash") != task["expected_protocol_hash"]:
        return None
    if task.get("require_solution_diagnostics") and not solution_diagnostics_ok(record):
        return None
    for phase in ("smd", "gas"):
        try:
            hartree = record[f"e_{phase}_Ha"]
            if (type(hartree) not in (float, int) or not math.isfinite(hartree)
                    or abs(record[f"e_{phase}_eV"] - hartree * HARTREE_EV) > 1e-5):
                return None
        except (KeyError, TypeError):
            return None
    return record


def inputs_current(task):
    try:
        return all(sha(path) == digest for path, digest in task["inputs_sha256"].items())
    except OSError:
        return False


def done_valid(run, task):
    done = read_json(Path(run) / "done" / f"{task['key']}.json")
    candidate = record_path(Path(task["frozen"]), task["charge"], task["mult"])
    return bool(done and candidate is not None and candidate.is_file()
                and done.get("record_sha256") == sha(candidate)
                and verified_candidate(task) is not None)


def compute_task(task, backend="gpu", execution_context=None):
    from redox.qm import dft, sp
    frozen = Path(task["frozen"])
    existing = verified_candidate(task)
    if existing is not None:
        return existing
    if backend not in ("cpu", "gpu"):
        raise ValueError("Revalidation backend must be cpu or gpu")
    from redox.core import protocol
    old_cleanup, old_versions = dft._free_gpu_pool, protocol.software_versions
    cpu_versions = dict(pyscf=metadata.version("pyscf"), gpu4pyscf=None)
    try:
        if backend == "cpu":
            # GPU cleanup and version imports are irrelevant to CPU work and
            # must not create CUDA contexts on devices owned by other jobs.
            dft._free_gpu_pool = lambda: None
            protocol.software_versions = lambda: cpu_versions
        output = dft.dft_smd(frozen / "opt.xyz", task["charge"], task["mult"],
                         sp_xc=ACTIVE_SP["xc"], sp_basis=ACTIVE_SP["basis"],
                         sp_basis_anion=ACTIVE_SP["basis"], sp_disp=ACTIVE_SP["disp"],
                         sp_nlc=ACTIVE_SP["nlc"], solvent=ACTIVE_SP["solvent"],
                         do_opt=False, do_freq=False, do_gas=True, do_smd=True, backend=backend)
    finally:
        dft._free_gpu_pool, protocol.software_versions = old_cleanup, old_versions
    if (output.get("energy_protocol") != ACTIVE_SP
            or output.get("energy_geom_sha1") != geom_sha1(frozen / "opt.xyz")):
        raise RuntimeError("Producer protocol differs from the frozen revalidation protocol")
    path = record_path(frozen, task["charge"], task["mult"])
    record = dict(protocol=ACTIVE_SP, protocol_hash=path.stem,
                  geom_sha1=geom_sha1(frozen / "opt.xyz"), geom=str(frozen / "opt.xyz"),
                  charge=task["charge"], mult=task["mult"], id=task["id"], state=task["state"],
                  source=f"computed:{backend}:independent_revalidation",
                  software=cpu_versions if backend == "cpu" else software_versions(),
                  producer_source_sha256=sha(ROOT / "src/redox/qm/dft.py"),
                  **{key: output.get(key) for key in sp._ENERGY_KEYS})
    if execution_context:
        record["execution_context"] = dict(execution_context)
    atomic_write_json(path, record)
    if verified_candidate(task) is None:
        raise RuntimeError("Recomputed SP failed convergence, identity or energy-unit checks")
    return record


def worker(run, name, queue_file="queue.json"):
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner()
    run = Path(run).resolve()
    if Path(queue_file).name != queue_file:
        raise ValueError("Queue selection must be a file in this run")
    queue = read_json(run / queue_file)
    if queue["protocol"] != ACTIVE_SP or any(sha(path) != digest for path, digest in queue["code_sha256"].items()):
        raise RuntimeError("Prepared scientific code or protocol changed")
    execution = queue.get("execution_context", {})
    failure_dir = run / queue.get("failure_directory", "failures")
    failure_dir.mkdir(exist_ok=True)
    if execution.get("density_fit_tensors") == "host":
        from redox.workflows.resources import host_capacity
        host_capacity(2, minimum_ram_gib=64)
        certificate = read_json(execution["qualification_path"])
        if (sha(execution["qualification_path"]) != execution["qualification_sha256"]
                or certificate.get("status") != "verified_record"
                or max(abs(value) for value in certificate["control_differences_eV"].values()) > .001):
            raise RuntimeError("Host density-fitting storage has not passed its control")
        from gpu4pyscf.df import df
        from redox.qm import dft
        df.DF.use_gpu_memory = False
        old_free = dft._free_gpu_pool

        def free_cycles_and_pool():
            gc.collect()
            old_free()

        dft._free_gpu_pool = free_cycles_and_pool
    owner = dict(host=socket.gethostname().split(".")[0], pid=os.getpid(), worker=name,
                 process_start=Path("/proc/self/stat").read_text().rsplit(")", 1)[1].split()[19],
                 affinity=sorted(os.sched_getaffinity(0)), started=time.time())
    status_path = run / "workers" / f"{name}.json"
    status = dict(owner, status="running", completed=0, failures=0, active_task=None)
    claim = None
    consecutive = 0

    def stop(signum, frame):
        raise KeyboardInterrupt(f"signal {signum}")

    signal.signal(signal.SIGTERM, stop)
    try:
        while True:
            selected = None
            waiting = False
            for task in queue["tasks"]:
                if done_valid(run, task):
                    continue
                failure = read_json(failure_dir / f"{task['key']}.json", {})
                if failure.get("attempts", 0) >= queue["max_attempts"]:
                    continue
                if failure.get("retry_after", 0) > time.time():
                    waiting = True
                    continue
                proposed = run / "claims" / task["key"]
                try:
                    proposed.mkdir()
                except FileExistsError:
                    waiting = True
                    continue
                claim = proposed
                atomic_write_json(claim / "owner.json", owner)
                selected = task
                break
            if selected is None:
                if not waiting:
                    break
                atomic_write_json(status_path, dict(status, status="waiting_for_peers_or_retry"))
                time.sleep(5)
                continue
            task = selected
            status.update(status="computing", active_task=task["key"], task_started=time.time())
            atomic_write_json(status_path, status)
            print(f"[revalidate] START {task['key']} q={task['charge']} m={task['mult']}", flush=True)
            try:
                if not inputs_current(task):
                    raise RuntimeError("Canonical source changed; prepare a fresh task")
                compute_task(task, execution_context=execution)
                if not inputs_current(task):
                    raise RuntimeError("Canonical source changed during independent evaluation")
                path = record_path(Path(task["frozen"]), task["charge"], task["mult"])
                atomic_write_json(run / "done" / f"{task['key']}.json",
                                  dict(task=task["key"], record_sha256=sha(path), owner=owner, finished=time.time()))
                status["completed"] += 1
                consecutive = 0
                print(f"[revalidate] VERIFIED RECORD {task['key']}", flush=True)
            except Exception as exc:
                failure_path = failure_dir / f"{task['key']}.json"
                previous = read_json(failure_path, {})
                atomic_write_json(failure_path, dict(attempts=previous.get("attempts", 0) + 1,
                                  retry_after=time.time() + 60, error=repr(exc), owner=owner))
                traceback.print_exc()
                status["failures"] += 1
                consecutive += 1
            finally:
                shutil.rmtree(claim)
                claim = None
                status.update(active_task=None, updated=time.time())
                atomic_write_json(status_path, status)
            if consecutive >= 3:
                raise RuntimeError("Three consecutive actual failures; worker circuit breaker")
        status["status"] = "queue_terminal"
    except BaseException as exc:
        status.update(status="stopped_preserving_outputs", error=repr(exc))
        raise
    finally:
        if claim is not None:
            shutil.rmtree(claim)
        status["finished"] = time.time()
        atomic_write_json(status_path, status)


def collect(run, promote=False):
    run = Path(run).resolve()
    queue = read_json(run / "queue.json")
    rows = []
    for task in queue["tasks"]:
        row = dict(id=task["id"], state=task["state"], status="pending", promoted=False)
        if not done_valid(run, task):
            failure = read_json(run / "failures" / f"{task['key']}.json", {})
            row.update(attempts=failure.get("attempts", 0), error=failure.get("error", ""))
            rows.append(row)
            continue
        record = verified_candidate(task)
        previous = read_json(Path(task["frozen"]) / "previous_record.json")
        for phase in ("smd", "gas"):
            row[f"delta_{phase}_eV"] = record[f"e_{phase}_eV"] - previous[f"e_{phase}_eV"]
        agrees = max(abs(row[f"delta_{phase}_eV"]) for phase in ("smd", "gas")) <= queue["agreement_eV"]
        row["status"] = "independently_reproduced" if agrees else "changed_energy_requires_review"
        if promote and agrees:
            lockdir = ROOT / "logs/pipeline_locks"
            lockdir.mkdir(parents=True, exist_ok=True)
            with (lockdir / f"{task['id']}.lock").open("a+") as lock:
                fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                canonical = Path(task["old_record"])
                candidate = record_path(Path(task["frozen"]), task["charge"], task["mult"])
                evidence = dict(run=str(run), task=task["key"], record_sha256=sha(candidate),
                                original_sha256=task["inputs_sha256"][str(canonical)])
                current = read_json(canonical, {})
                geometry_current = all(sha(path) == digest for path, digest in task["inputs_sha256"].items()
                                       if path != str(canonical))
                if not geometry_current:
                    row["status"] = "source_changed_promotion_refused"
                elif current.get("independent_revalidation") == evidence:
                    row["promoted"] = True
                elif inputs_current(task):
                    archive = canonical.parent / "_superseded" / f"{canonical.stem}.{evidence['original_sha256']}.json"
                    archive.parent.mkdir(exist_ok=True)
                    if archive.exists() and sha(archive) != evidence["original_sha256"]:
                        raise RuntimeError("Historical SP archive collision")
                    if not archive.exists():
                        shutil.copy2(canonical, archive)
                    replacement = dict(record, independent_revalidation=evidence)
                    atomic_write_json(canonical, replacement)
                    if load_record(Path(task["state_dir"]), task["charge"], task["mult"]) != replacement:
                        raise RuntimeError("Promoted record failed current-address verification")
                    row["promoted"] = True
                else:
                    row["status"] = "source_changed_promotion_refused"
        rows.append(row)
    columns = ["id", "state", "status", "delta_smd_eV", "delta_gas_eV", "promoted", "attempts", "error"]
    with atomic_text_writer(run / "comparison.csv", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader(); writer.writerows(rows)
    summary = dict(tasks=len(rows), verified=sum(row["status"] == "independently_reproduced" for row in rows),
                   changed=sum(row["status"] == "changed_energy_requires_review" for row in rows),
                   promoted=sum(row["promoted"] for row in rows), pending=sum(row["status"] == "pending" for row in rows),
                   unresolved=sum(row["status"] != "independently_reproduced" for row in rows),
                   scientific_accuracy_validated=False)
    atomic_write_json(run / "summary.json", summary)
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="mode", required=True)
    plan = subparsers.add_parser("prepare"); plan.add_argument("--out", type=Path, required=True)
    plan.add_argument("--selection", choices=("adopted", "missing-smd-diagnostics"), default="adopted")
    plan.add_argument("--only", help="comma-separated molecule IDs or id:state selections")
    work = subparsers.add_parser("worker"); work.add_argument("--run", type=Path, required=True); work.add_argument("--name", required=True)
    work.add_argument("--queue-file", default="queue.json")
    gather = subparsers.add_parser("collect"); gather.add_argument("--run", type=Path, required=True); gather.add_argument("--promote", action="store_true")
    args = parser.parse_args(argv)
    if args.mode == "prepare":
        rows = read_manifest()
        if args.only:
            selected = set(args.only.split(","))
            available = {key for row in rows for key in (row["id"], f"{row['id']}:{row['state']}")}
            if selected - available:
                parser.error(f"unknown selections: {sorted(selected - available)}")
            rows = [row for row in rows if row["id"] in selected or f"{row['id']}:{row['state']}" in selected]
        queue = prepare(args.out, rows=rows, selection=args.selection)
        print(json.dumps(dict(tasks=len(queue["tasks"]), run=str(args.out), status=queue["status"])))
    elif args.mode == "worker":
        worker(args.run, args.name, args.queue_file)
    else:
        print(json.dumps(collect(args.run, args.promote)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
