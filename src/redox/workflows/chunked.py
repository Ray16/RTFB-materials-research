"""Bounded-memory task packing and atomic chunk ownership on shared storage.

Cheap in-process stages use chunk journals instead of per-molecule control files.
Quantum adapters must retain their separate reservations and physical admission;
an executor completion record is never a quantum-quality certificate.
"""

from contextlib import contextmanager
import hashlib
import json
import os
from pathlib import Path
import resource
import signal
import socket
import sqlite3
import tempfile
import time
import uuid

from redox.core.io import atomic_write_json

SOURCE_BYTES = Path(__file__).read_bytes()
SOURCE_SHA256 = hashlib.sha256(SOURCE_BYTES).hexdigest()


def canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def content_hash(value):
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def file_hash(path):
    result = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def local_scratch(path):
    """Keep the mutable SQLite dedup index off shared filesystems."""
    path = Path(path).resolve()
    mounts = []
    for line in Path("/proc/self/mountinfo").read_text().splitlines():
        fields = line.split()
        mount = fields[4].replace("\\040", " ")
        if str(path) == mount or str(path).startswith(mount.rstrip("/") + "/"):
            mounts.append((len(mount), fields[fields.index("-") + 1]))
    if not mounts or max(mounts)[1] in {"nfs", "nfs4", "cifs", "lustre", "gpfs"}:
        raise ValueError("Deduplication requires explicitly node-local scratch")
    path.mkdir(parents=True, exist_ok=True)
    return path


def task_record(record, stage):
    body = dict(
        identity=record["identity"],
        payload=record["payload"],
        stage=stage if isinstance(stage, str) else content_hash(stage),
    )
    return dict(task_id=content_hash(body), **body)


def prepare(
    records,
    destination,
    stage,
    *,
    chunk_size=1000,
    scratch="/tmp",
    max_chunk_bytes=8 * 1024**2,
):
    """Consume an iterator once; use local disk for exact dedup and stable sorting.

    Every recipe is kept in lineage.jsonl, including duplicate task inputs.
    Chunk ordering is independent of input order. Memory is bounded by SQLite's
    32-MiB page cache, one record and the O(number-of-chunks) manifest.
    """
    if type(chunk_size) is not int or chunk_size < 1:
        raise ValueError("Positive chunk size required")
    if type(max_chunk_bytes) is not int or max_chunk_bytes < 1:
        raise ValueError("Positive chunk byte limit required")
    if (
        not stage.get("name")
        or not stage.get("protocol")
        or not stage.get("quality_scope")
    ):
        raise ValueError("Named stage, protocol and quality scope required")
    destination = Path(destination).resolve()
    destination.mkdir(parents=True, exist_ok=False)
    for folder in (
        "inputs",
        "pending",
        "workers",
        "progress",
        "results",
        "done",
        "failures",
        "history",
    ):
        (destination / folder).mkdir()
    (destination / "producer_source.py").write_bytes(SOURCE_BYTES)
    started = time.time()
    atomic_write_json(
        destination / "preparation.json", dict(status="preparing", started=started)
    )
    counts = dict(input_records=0, unique_tasks=0, duplicate_recipes=0)
    chunks = []
    stage_reference = content_hash(stage)
    try:
        with tempfile.TemporaryDirectory(
            prefix="redox_dedup_", dir=local_scratch(scratch)
        ) as temporary:
            with sqlite3.connect(str(Path(temporary) / "dedup.sqlite")) as database:
                database.execute("PRAGMA cache_size=-32768")
                database.execute("PRAGMA temp_store=FILE")
                database.execute(
                    "CREATE TABLE tasks (id TEXT PRIMARY KEY, payload TEXT NOT NULL)"
                )
                with (destination / "lineage.jsonl").open("w") as lineage:
                    for record in records:
                        task = task_record(record, stage_reference)
                        encoded = canonical(task)
                        if len(encoded.encode()) + 1 > max_chunk_bytes:
                            raise ValueError("Single task exceeds the chunk byte limit")
                        inserted = database.execute(
                            "INSERT OR IGNORE INTO tasks VALUES (?, ?)",
                            (task["task_id"], encoded),
                        ).rowcount
                        counts["input_records"] += 1
                        if inserted:
                            counts["unique_tasks"] += 1
                        else:
                            previous = database.execute(
                                "SELECT payload FROM tasks WHERE id=?",
                                (task["task_id"],),
                            ).fetchone()[0]
                            if previous != encoded:
                                raise ValueError(
                                    "Task hash collision with different exact inputs"
                                )
                            counts["duplicate_recipes"] += 1
                        lineage.write(
                            canonical(
                                dict(
                                    task_id=task["task_id"],
                                    provenance=record.get("provenance", {}),
                                )
                            )
                            + "\n"
                        )
                        if counts["input_records"] % 4096 == 0:
                            database.commit()
                    lineage.flush()
                    os.fsync(lineage.fileno())
                database.commit()
                stream = None
                try:
                    chunk_bytes = 0
                    for _, encoded in database.execute(
                        "SELECT id,payload FROM tasks ORDER BY id"
                    ):
                        encoded_bytes = len(encoded.encode()) + 1
                        if (
                            stream is None
                            or chunks[-1]["count"] >= chunk_size
                            or chunk_bytes + encoded_bytes > max_chunk_bytes
                        ):
                            if stream is not None:
                                stream.flush()
                                os.fsync(stream.fileno())
                                stream.close()
                                chunks[-1]["sha256"] = file_hash(
                                    destination / chunks[-1]["path"]
                                )
                            key = f"chunk_{len(chunks):08d}"
                            relative = f"inputs/{key}.jsonl"
                            stream = (destination / relative).open("w")
                            chunks.append(dict(key=key, path=relative, count=0))
                            chunk_bytes = 0
                        stream.write(encoded + "\n")
                        chunk_bytes += encoded_bytes
                        chunks[-1]["count"] += 1
                    if stream is not None:
                        stream.flush()
                        os.fsync(stream.fileno())
                        stream.close()
                        stream = None
                        chunks[-1]["sha256"] = file_hash(
                            destination / chunks[-1]["path"]
                        )
                finally:
                    if stream is not None:
                        stream.close()
        plan = dict(
            schema=1,
            stage=stage,
            chunk_size=chunk_size,
            max_chunk_bytes=max_chunk_bytes,
            chunks=chunks,
            **counts,
            lineage_sha256=file_hash(destination / "lineage.jsonl"),
            executor_sha256=SOURCE_SHA256,
            prepared_seconds=time.time() - started,
        )
        plan["manifest_sha256"] = content_hash(plan)
        atomic_write_json(destination / "plan.json", plan)
        for chunk in chunks:
            atomic_write_json(destination / "pending" / (chunk["key"] + ".json"), chunk)
        atomic_write_json(
            destination / "preparation.json",
            dict(status="prepared", **counts, chunks=len(chunks)),
        )
        return plan
    except BaseException as error:
        atomic_write_json(
            destination / "preparation.json",
            dict(
                status="interrupted_preparation_preserved", error=repr(error), **counts
            ),
        )
        raise


def read_plan(campaign, *, require_current=False):
    preparation = json.loads((Path(campaign) / "preparation.json").read_text())
    if preparation["status"] != "prepared":
        raise ValueError("Task packing is not complete")
    plan = json.loads((Path(campaign) / "plan.json").read_text())
    body = {k: v for k, v in plan.items() if k != "manifest_sha256"}
    if content_hash(body) != plan["manifest_sha256"]:
        raise ValueError("Chunk manifest or executor source changed")
    current = plan["executor_sha256"] == SOURCE_SHA256 == file_hash(__file__)
    if require_current and not current:
        raise ValueError("Worker executor differs from frozen producer")
    if not current:
        snapshot = Path(campaign) / "producer_source.py"
        if not snapshot.exists() or file_hash(snapshot) != plan["executor_sha256"]:
            raise ValueError("Historical executor source snapshot missing or changed")
    return plan


def process_start(pid):
    return Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()[19]


@contextmanager
def owner(campaign):
    """Record ownership before any atomic pending->worker claim rename."""
    campaign = Path(campaign)
    host = socket.gethostname().split(".")[0]
    identity = dict(
        host=host,
        pid=os.getpid(),
        process_start=process_start(os.getpid()),
        nonce=uuid.uuid4().hex,
        cpus=sorted(os.sched_getaffinity(0)),
        executor_sha256=file_hash(__file__),
        thread_caps={
            k: os.environ.get(k)
            for k in (
                "OMP_NUM_THREADS",
                "OMP_THREAD_LIMIT",
                "MKL_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "NUMEXPR_NUM_THREADS",
                "XTB_THREADS",
            )
        },
    )
    key = f"{host}_{identity['pid']}_{identity['process_start']}_{identity['nonce']}"
    directory = campaign / "workers" / key
    (directory / "claims").mkdir(parents=True)
    atomic_write_json(directory / "owner.json", identity)
    try:
        yield directory
    finally:
        atomic_write_json(directory / "exit.json", dict(exited=time.time()))


def claim(campaign, worker, chunks, cursor=0):
    """Atomic NFS rename gives one claimant; no global mutable database or counter."""
    campaign, worker = Path(campaign), Path(worker)
    for offset in range(len(chunks)):
        index = (cursor + offset) % len(chunks)
        chunk = chunks[index]
        source = campaign / "pending" / (chunk["key"] + ".json")
        target = worker / "claims" / source.name
        try:
            source.rename(target)
        except FileNotFoundError:
            continue
        if json.loads(target.read_text()) != chunk:
            raise ValueError("Claim differs from the immutable chunk manifest")
        return chunk, (index + 1) % len(chunks)
    return None, cursor


def journal(campaign, chunk):
    path = Path(campaign) / "progress" / (chunk["key"] + ".json")
    if not path.exists():
        return {}
    saved = json.loads(path.read_text())
    if saved["chunk_sha256"] != chunk["sha256"] or saved["rows_sha256"] != content_hash(
        saved["rows"]
    ):
        raise ValueError("Checkpoint input or output integrity mismatch")
    if any(value.get("task_id") != key for key, value in saved["rows"].items()):
        raise ValueError("Checkpoint task identities do not match their keys")
    return saved["rows"]


def save_journal(campaign, chunk, rows):
    atomic_write_json(
        Path(campaign) / "progress" / (chunk["key"] + ".json"),
        dict(chunk_sha256=chunk["sha256"], rows_sha256=content_hash(rows), rows=rows),
    )
    with (Path(campaign) / "progress" / (chunk["key"] + ".json")).open("rb") as stream:
        os.fsync(stream.fileno())


def run_chunk(campaign, chunk, handler, *, checkpoint_every=100):
    if checkpoint_every < 1:
        raise ValueError("Positive checkpoint interval required")
    campaign = Path(campaign)
    path = campaign / chunk["path"]
    if file_hash(path) != chunk["sha256"]:
        raise ValueError("Task chunk changed")
    tasks = [json.loads(line) for line in path.open()]
    if len(tasks) != chunk["count"] or len({t["task_id"] for t in tasks}) != len(tasks):
        raise ValueError("Task count or uniqueness mismatch")
    rows = journal(campaign, chunk)
    if set(rows) - {t["task_id"] for t in tasks}:
        raise ValueError("Checkpoint contains an unrelated task")
    consecutive = 0
    since_save = 0
    try:
        for task in tasks:
            key = task.pop("task_id")
            if content_hash(task) != key:
                raise ValueError("Task identity differs from its exact inputs")
            if key in rows:
                continue
            try:
                value = handler(task)
                canonical(value)
                rows[key] = dict(task_id=key, status="completed", value=value)
                consecutive = 0
            except Exception as error:
                rows[key] = dict(task_id=key, status="failed", error=repr(error))
                consecutive += 1
            since_save += 1
            if since_save >= checkpoint_every:
                save_journal(campaign, chunk, rows)
                since_save = 0
            if consecutive >= 3:
                raise RuntimeError(
                    "Three consecutive actual task failures; circuit breaker"
                )
    finally:
        save_journal(campaign, chunk, rows)
    result = campaign / "results" / (chunk["key"] + ".jsonl")
    with result.open("w") as stream:
        for key in sorted(rows):
            stream.write(canonical(rows[key]) + "\n")
        stream.flush()
        os.fsync(stream.fileno())
    failed = sum(row["status"] == "failed" for row in rows.values())
    report = dict(
        status="executed_with_failures" if failed else "executed",
        chunk=chunk,
        results=str(result),
        results_sha256=file_hash(result),
        completed=len(rows) - failed,
        failed=failed,
        scientific_QC_inferred_from_execution=False,
    )
    atomic_write_json(campaign / "done" / (chunk["key"] + ".json"), report)
    return report


def verify_done(path, *, plan=None):
    report = json.loads(Path(path).read_text())
    campaign = Path(path).parent.parent
    plan = read_plan(campaign) if plan is None else plan
    expected = next(
        (c for c in plan["chunks"] if c["key"] == report["chunk"]["key"]), None
    )
    if expected != report["chunk"]:
        raise ValueError("Completion record differs from the frozen chunk")
    if file_hash(report["results"]) != report["results_sha256"]:
        raise ValueError("Chunk results changed after completion")
    if (
        file_hash(Path(path).parent.parent / report["chunk"]["path"])
        != report["chunk"]["sha256"]
    ):
        raise ValueError("Completed chunk inputs changed")
    if report["completed"] + report["failed"] != report["chunk"]["count"]:
        raise ValueError("Chunk completion count mismatch")
    with Path(report["results"]).open() as stream:
        rows = [json.loads(line) for line in stream]
    with (campaign / report["chunk"]["path"]).open() as stream:
        inputs = {json.loads(line)["task_id"] for line in stream}
    if (
        len(rows) != len(inputs)
        or {r["task_id"] for r in rows} != inputs
        or sum(r["status"] == "completed" for r in rows) != report["completed"]
        or sum(r["status"] == "failed" for r in rows) != report["failed"]
    ):
        raise ValueError("Actual completed task identities or statuses disagree")
    return report


def run_worker(campaign, handler, *, checkpoint_every=100):
    """Run synchronous CPU tasks with dynamic chunk allocation and one CPU."""
    from redox.workflows.resources import (
        allocation_guard,
        cpu_allocation,
        host_capacity,
    )

    campaign = Path(campaign).resolve()
    plan = read_plan(campaign, require_current=True)
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    for name in (
        "OMP_NUM_THREADS",
        "OMP_THREAD_LIMIT",
        "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS",
        "NUMEXPR_NUM_THREADS",
        "XTB_THREADS",
    ):
        if os.environ.get(name) != "1":
            raise ValueError("Cheap CPU stages require single-thread native pools")
    previous = signal.getsignal(signal.SIGTERM)

    def interrupted(sig, frame):
        raise KeyboardInterrupt(f"signal {sig}")

    signal.signal(signal.SIGTERM, interrupted)
    try:
        with owner(campaign) as directory:
            state = dict(status="running", started=time.time(), completed_chunks=[])
            cursor = int(content_hash(directory.name)[:8], 16) % max(
                1, len(plan["chunks"])
            )
            atomic_write_json(directory / "status.json", state)
            try:
                while True:
                    chunk, cursor = claim(campaign, directory, plan["chunks"], cursor)
                    if chunk is None:
                        break
                    state["active_chunk"] = chunk["key"]
                    atomic_write_json(directory / "status.json", state)
                    done = campaign / "done" / (chunk["key"] + ".json")
                    if done.exists():
                        verify_done(done)
                    else:
                        run_chunk(
                            campaign, chunk, handler, checkpoint_every=checkpoint_every
                        )
                    (directory / "claims" / (chunk["key"] + ".json")).unlink()
                    state["completed_chunks"].append(chunk["key"])
                    atomic_write_json(directory / "status.json", state)
                state.update(status="worker_finished", finished=time.time())
            except BaseException as error:
                state.update(
                    status="stopped_preserving_outputs",
                    error=repr(error),
                    finished=time.time(),
                )
                raise
            finally:
                usage = resource.getrusage(resource.RUSAGE_SELF)
                state.update(
                    peak_RSS_KiB=usage.ru_maxrss,
                    user_cpu_seconds=usage.ru_utime,
                    system_cpu_seconds=usage.ru_stime,
                )
                atomic_write_json(directory / "status.json", state)
    finally:
        signal.signal(signal.SIGTERM, previous)


def recover_local(campaign, worker):
    """Recover an in-process CPU worker only after exact local ownership is dead.

    Remote recovery requires a successful probe on that host; run this function
    there. Registered handlers must not detach child predictors. Quantum stages
    need their separate predictor/gate ownership audit before using this backend.
    """
    campaign, worker = Path(campaign), Path(worker)
    identity = json.loads((worker / "owner.json").read_text())
    if identity["host"] != socket.gethostname().split(".")[0]:
        raise ValueError("Ownership must be audited on its actual reachable host")
    try:
        live = process_start(identity["pid"]) == identity["process_start"]
    except FileNotFoundError:
        live = False
    if live:
        raise ValueError("Owner is still live; claim recovery refused")
    recovered = []
    for path in (worker / "claims").glob("*.json"):
        done = campaign / "done" / path.name
        if done.exists():
            verify_done(done)
            path.unlink()
        else:
            target = campaign / "pending" / path.name
            if target.exists():
                raise ValueError("Duplicate pending claim")
            path.rename(target)
            recovered.append(path.stem)
    atomic_write_json(
        worker / "recovery.json",
        dict(recovered=recovered, checked=time.time(), owner=identity),
    )
    return recovered
