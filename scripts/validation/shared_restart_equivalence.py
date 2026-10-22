#!/usr/bin/env python
"""Real descriptor execution, abrupt process death and cross-host NFS recovery."""

import argparse
from itertools import islice
import json
import os
from pathlib import Path

from redox.core.io import atomic_write_json
from redox.screening import streaming_search as search
from redox.workflows import chunked as queue


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "crash", "resume", "baseline", "recover", "audit"))
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity

    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    root = args.root.resolve()
    if args.action == "prepare":
        root.mkdir()
        records = list(islice(search.family_records(), 32))
        for name in ("baseline", "resumed"):
            queue.prepare(records, root / name, search.stage_spec("descriptors"), chunk_size=32)
        (root / "driver_source.py").write_bytes(Path(__file__).read_bytes())
        return
    name = "baseline" if args.action == "baseline" else "resumed"
    campaign = root / name
    plan = queue.read_plan(campaign, require_current=True)
    if args.action == "recover":
        owners = list((campaign / "workers").iterdir())
        if len(owners) != 1:
            raise ValueError("Expected the one abruptly terminated worker")
        before = queue.journal(campaign, plan["chunks"][0])
        if len(before) != 4:
            raise ValueError("Expected exactly four durable checked tasks")
        atomic_write_json(root / "before_recovery.json", before)
        recovered = queue.recover_local(campaign, owners[0])
        atomic_write_json(root / "recovery.json", dict(recovered=recovered, previous_owner=str(owners[0])))
        return
    if args.action == "audit":
        def outputs(folder):
            records = {}
            for done in (folder / "done").glob("*.json"):
                report = queue.verify_done(done)
                for line in Path(report["results"]).read_text().splitlines():
                    row = json.loads(line)
                    records[row["task_id"]] = row
            return records

        actual, baseline = outputs(campaign), outputs(root / "baseline")
        before = json.loads((root / "before_recovery.json").read_text())
        after = queue.journal(campaign, plan["chunks"][0])
        if len(actual) != plan["unique_tasks"] or actual != baseline:
            raise ValueError("Cross-host resumed outputs differ from uninterrupted execution")
        if not all(after[key] == value for key, value in before.items()):
            raise ValueError("Previously checked task changed during recovery")
        called = json.loads((root / "resume_calls.json").read_text())
        if len(called) != len(actual) - 4 or len(set(called)) != len(called) or set(called) & set(before):
            raise ValueError("Checked tasks repeated or missing tasks were skipped")
        # Independently replay serial descriptor definitions on the exact inputs.
        for chunk in plan["chunks"]:
            for line in (campaign / chunk["path"]).read_text().splitlines():
                task = json.loads(line)
                key = task.pop("task_id")
                if actual[key]["value"] != search.descriptor(task):
                    raise ValueError("Resumed descriptor differs from serial definition")
        atomic_write_json(root / "audit.json", dict(status="cross_host_shared_filesystem_restart_passed",
                          unique_tasks=len(actual), checked_before_crash=4, recomputed_after_restart=len(called),
                          repeated_checked_tasks=0, differences=0, physical_quantum_recovery_validated=False,
                          driver_sha256=queue.file_hash(__file__), resumed_manifest_sha256=plan["manifest_sha256"]))
        return
    calls = []

    def handler(task):
        calls.append(queue.content_hash(task))
        if args.action == "crash" and len(calls) == 5:
            os._exit(9)
        atomic_write_json(root / (args.action + "_calls.json"), calls)
        return search.descriptor(task)

    with queue.owner(campaign) as owner:
        chunk, _ = queue.claim(campaign, owner, plan["chunks"])
        if chunk is None:
            raise ValueError("Expected one uncompleted chunk")
        queue.run_chunk(campaign, chunk, handler, checkpoint_every=1)


if __name__ == "__main__":
    main()
