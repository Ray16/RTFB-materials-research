#!/usr/bin/env python
"""Independently check actual chunk results against exact inputs and known controls."""

import argparse
from functools import cache
import hashlib
import json
import os
from pathlib import Path
import time

from redox.core.io import atomic_write_json
from redox.workflows import chunked
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


@cache
def serial_tools():
    from redox.properties import capacity_and_proxies as capacity
    from redox.core.common import load_config

    return capacity, capacity._sascorer(), load_config("electrolyte").FARADAY


def serial_descriptor(source):
    from redox.screening.discovery import describe

    capacity, scorer, faraday = serial_tools()
    previous = capacity._sascorer
    capacity._sascorer = lambda: scorer
    try:
        payload = source["payload"]
        states = payload["states"]
        electrons = max(s[3] for s in states) - min(s[3] for s in states)
        return describe(payload["product_smiles"], electrons, faraday)
    finally:
        capacity._sascorer = previous


def audit(campaign, out, reference=None):
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    campaign = Path(campaign).resolve()
    plan = chunked.read_plan(campaign)
    expected = (
        json.loads(Path(reference).read_text())["expected"] if reference else None
    )
    checked = 0
    errors = []
    total_index = 0
    total_square = 0
    result_hashes = {}
    started = time.time()
    for part in plan["chunks"]:
        report = chunked.verify_done(
            campaign / "done" / (part["key"] + ".json"), plan=plan
        )
        if report["failed"]:
            raise ValueError("Benchmark contains failed tasks")
        with (campaign / part["path"]).open() as stream:
            inputs = {row["task_id"]: row for row in map(json.loads, stream)}
        with Path(report["results"]).open() as stream:
            for result in map(json.loads, stream):
                key = result["task_id"]
                source = inputs[key]
                exact = {k: v for k, v in source.items() if k != "task_id"}
                independently_hashed = hashlib.sha256(
                    json.dumps(
                        exact, sort_keys=True, separators=(",", ":"), allow_nan=False
                    ).encode()
                ).hexdigest()
                if independently_hashed != key:
                    raise ValueError("Input task content does not match its identity")
                value = result["value"]
                if plan["stage"]["name"] == "scheduler_control":
                    index = source["payload"]["index"]
                    if value["index"] != index or value["square"] != index**2:
                        errors.append(key)
                    total_index += index
                    total_square += value["square"]
                else:
                    if expected is not None and key not in expected:
                        raise ValueError("Exact serial molecular reference absent")
                    baseline = (
                        expected[key]
                        if expected is not None
                        else serial_descriptor(source)
                    )
                    if any(value.get(k) != v for k, v in baseline.items()):
                        errors.append(key)
                    if value.get("solution_properties_status") != "not computed":
                        raise ValueError(
                            "Cheap execution acquired an unsupported "
                            "solution-property claim"
                        )
                checked += 1
        result_hashes[report["results"]] = report["results_sha256"]
    if checked != plan["unique_tasks"] or (
        expected is not None and checked != len(expected)
    ):
        raise ValueError("Actual benchmark coverage mismatch")
    if plan["stage"]["name"] == "scheduler_control":
        n = checked
        if (
            total_index != n * (n - 1) // 2
            or total_square != n * (n - 1) * (2 * n - 1) // 6
        ):
            raise ValueError("Independent aggregate arithmetic failed")
    workers = []
    for directory in (campaign / "workers").iterdir():
        owner = json.loads((directory / "owner.json").read_text())
        status = json.loads((directory / "status.json").read_text())
        workers.append(dict(owner=owner, status=status))
    files = [p for p in campaign.rglob("*") if p.is_file()]
    report = dict(
        status="verified_execution_control" if not errors else "failed_control",
        checked_tasks=checked,
        discrepancies=len(errors),
        error_examples=errors[:20],
        stage=plan["stage"],
        worker_count=len(workers),
        workers=workers,
        chunks=len(plan["chunks"]),
        input_records=plan["input_records"],
        duplicate_recipes=plan["duplicate_recipes"],
        raw_root=str(campaign),
        raw_file_count=len(files),
        raw_bytes=sum(p.stat().st_size for p in files),
        plan_sha256=chunked.file_hash(campaign / "plan.json"),
        result_sha256=result_hashes,
        reference_sha256=chunked.file_hash(reference) if reference else None,
        serial_reference_mode="frozen_fixture"
        if reference
        else "existing_serial_definition_on_exact_inputs",
        validator_sha256=chunked.file_hash(__file__),
        checked_at=time.time(),
        audit_seconds=time.time() - started,
        scientific_prediction_accuracy_validated=False,
        scope=(
            "Synthetic scheduler/storage and arithmetic control; "
            "no chemical throughput claim"
            if plan["stage"]["name"] == "scheduler_control"
            else "Exact serial/parallel bookkeeping descriptor equivalence; "
            "no redox or route validation"
        ),
    )
    atomic_write_json(out, report)
    print(
        json.dumps(
            {
                k: report[k]
                for k in (
                    "status",
                    "checked_tasks",
                    "discrepancies",
                    "chunks",
                    "raw_file_count",
                    "audit_seconds",
                )
            }
        )
    )
    if errors:
        raise ValueError("Benchmark values disagree with the independent control")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--reference", type=Path)
    args = parser.parse_args()
    audit(args.campaign, args.out, args.reference)
