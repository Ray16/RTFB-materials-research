"""Stream checked descriptor outputs into the existing exact-route gate.

This stage prepares a successor queue; it registers no candidates and submits no
quantum jobs. Each decision binds an exact precursor/product and attachment bond.
"""

import argparse
from contextlib import contextmanager
from functools import cache
import importlib.util
import json
import os
from pathlib import Path
import sqlite3
import tempfile

from importlib.metadata import version

from redox.build.attachment_identity import exact_attachment
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import chunked as queue


@cache
def gate():
    path = ROOT / "scripts/pipeline/screen_grafting_routes.py"
    spec = importlib.util.spec_from_file_location("streamed_prospective_route_gate", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.screen


@contextmanager
def catalog(path, kind, scratch="/tmp"):
    """Use node-local disk for large NDJSON catalogs; bound legacy JSON size."""
    path = Path(path)
    queue.local_scratch(Path(scratch))
    with tempfile.TemporaryDirectory(prefix="redox_route_catalog_", dir=scratch) as folder:
        database = sqlite3.connect(str(Path(folder) / "catalog.sqlite"))
        database.execute("PRAGMA cache_size=-8192")
        database.execute("CREATE TABLE entries (id TEXT PRIMARY KEY, value TEXT NOT NULL)")
        try:
            if path.suffix in (".jsonl", ".ndjson"):
                def rows():
                    with path.open() as stream:
                        for line in stream:
                            row = json.loads(line)
                            yield row["id"], row[kind]
            else:
                if path.stat().st_size > 16 * 1024**2:
                    raise ValueError("Use an NDJSON route catalog above 16 MiB")
                data = json.loads(path.read_text())["candidates"]
                def rows():
                    yield from data.items()
            for identifier, value in rows():
                if not isinstance(identifier, str) or not isinstance(value, dict):
                    raise ValueError("Catalog requires a string id and object value")
                database.execute("INSERT INTO entries VALUES (?,?)", (identifier, queue.canonical(value)))
            database.commit()

            def lookup(identifier):
                row = database.execute("SELECT value FROM entries WHERE id=?", (identifier,)).fetchone()
                return json.loads(row[0]) if row else {}

            yield lookup
        finally:
            database.close()


def source_paths():
    return [Path(__file__), ROOT / "src/redox/build/attachment_identity.py",
            ROOT / "src/redox/build/structures.py", ROOT / "src/redox/build/contracts.py",
            ROOT / "src/redox/build/grafting.py", ROOT / "scripts/pipeline/screen_grafting_routes.py"]


def stage(parent, evidence, assessment):
    parent = Path(parent).resolve()
    plan = queue.read_plan(parent)
    if plan["stage"]["name"] != "descriptors":
        raise ValueError("Route stage requires a molecular descriptor campaign")
    inputs = [Path(evidence).resolve(), Path(assessment).resolve(), parent / "plan.json"]
    for chunk in plan["chunks"]:
        receipt = parent / "done" / (chunk["key"] + ".json")
        if not receipt.exists():
            raise ValueError("Complete the descriptor campaign before freezing its successor")
        verified = queue.verify_done(receipt, plan=plan)
        if verified["failed"]:
            raise ValueError("Descriptor failures require review before a successor is frozen")
        inputs += [receipt, Path(verified["results"]), parent / chunk["path"]]
    return dict(name="exact_routes", quality_scope="prospective_route_review_only_no_quantum_submission",
                protocol=dict(schema=1, source_sha256={str(p): queue.file_hash(p) for p in source_paths()},
                              input_sha256={str(p): queue.file_hash(p) for p in inputs}, rdkit=version("rdkit"),
                              upstream_manifest_sha256=plan["manifest_sha256"]))


def records(parent, evidence, assessment, scratch="/tmp"):
    plan = queue.read_plan(parent)
    with catalog(evidence, "evidence", scratch) as review, catalog(assessment, "assessment", scratch) as selection:
        for chunk in plan["chunks"]:
            receipt = Path(parent) / "done" / (chunk["key"] + ".json")
            verified = queue.verify_done(receipt, plan=plan)
            values = {r["task_id"]: r for r in map(json.loads, Path(verified["results"]).read_text().splitlines())}
            for line in (Path(parent) / chunk["path"]).open():
                task = json.loads(line)
                descriptor_id = task.pop("task_id")
                if queue.content_hash(task) != descriptor_id:
                    raise ValueError("Upstream task differs from its exact identity")
                result = values[descriptor_id]
                if result["status"] != "completed":
                    raise ValueError("Unqualified descriptor task cannot advance")
                attachment = exact_attachment(task)
                identifier = "graft_" + attachment["exact_design_key"]
                yield dict(identity=dict(exact_design_key=attachment["exact_design_key"], parent_task_id=descriptor_id),
                           payload=dict(id=identifier, fragment=task["payload"]["attachment_fragment"],
                                        attachment=attachment, descriptor=result["value"],
                                        reviewed=review(identifier), assessment=selection(identifier)),
                           provenance=dict(parent_campaign=str(Path(parent).resolve()), parent_chunk=chunk["key"],
                                           parent_receipt_sha256=queue.file_hash(receipt)))


def screen(task):
    payload = task["payload"]
    identifier = payload["id"]
    reviewed = payload["reviewed"]
    decision = gate()(dict(id=identifier, frag=payload["fragment"]), {identifier: reviewed},
                      dict(candidates={identifier: payload["assessment"]}))
    expected = payload["attachment"]["route_binding"]
    exact = reviewed.get("route_binding") == expected
    if not exact:
        if decision["route_status"] == "proceed":
            decision["route_status"] = "hold"
        decision["proceed_to_full_redox"] = False
        decision["reasons"].append("Route evidence must bind this exact precursor, product and mapped new bond")
    return dict(decision, exact_design_key=payload["attachment"]["exact_design_key"],
                exact_route_binding_verified=exact, next_stage="eligible_for_declared_conformer_selection"
                if decision["proceed_to_full_redox"] else "route_review_or_exclusion",
                quantum_jobs_submitted=0, experimental_grafting_validated=False)


def verify(campaign):
    campaign = Path(campaign)
    plan = queue.read_plan(campaign)
    if plan["stage"]["name"] != "exact_routes":
        raise ValueError("Expected an exact-route campaign")
    counts = {}
    total = 0
    missing = []
    for chunk in plan["chunks"]:
        receipt = campaign / "done" / (chunk["key"] + ".json")
        if not receipt.exists():
            missing.append(chunk["key"])
            continue
        report = queue.verify_done(receipt, plan=plan)
        if report["failed"]:
            raise ValueError("Route-stage failures remain unqualified")
        for row in map(json.loads, Path(report["results"]).read_text().splitlines()):
            decision = row["value"]
            key = decision["route_status"]
            counts[key] = counts.get(key, 0) + 1
            total += 1
    result = dict(status="route_stage_complete" if not missing else "route_stage_pending", completed=total,
                  expected=plan["unique_tasks"], counts=counts, missing_chunks=missing, quantum_jobs_submitted=0)
    atomic_write_json(campaign / "verification.json", result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "worker", "verify"))
    parser.add_argument("--campaign", type=Path, required=True)
    parser.add_argument("--descriptors", type=Path)
    parser.add_argument("--evidence", type=Path, default=ROOT / "results/grafting/route_evidence.json")
    parser.add_argument("--assessment", type=Path, default=ROOT / "results/grafting/selection_assessment.json")
    parser.add_argument("--chunk-size", type=int, default=64)
    parser.add_argument("--checkpoint-every", type=int, default=100,
                        help="Completed tasks per durable checkpoint; compatible with the shared CPU launcher")
    args = parser.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    if args.action == "prepare":
        if args.descriptors is None:
            parser.error("--descriptors is required")
        spec = stage(args.descriptors, args.evidence, args.assessment)
        queue.prepare(records(args.descriptors, args.evidence, args.assessment), args.campaign, spec, chunk_size=args.chunk_size)
    elif args.action == "worker":
        plan = queue.read_plan(args.campaign, require_current=True)
        if plan["stage"]["name"] != "exact_routes" or plan["stage"]["protocol"]["rdkit"] != version("rdkit"):
            raise ValueError("Wrong route stage or chemical toolkit version")
        for group in ("source_sha256", "input_sha256"):
            if any(queue.file_hash(p) != h for p, h in plan["stage"]["protocol"][group].items()):
                raise ValueError("Frozen route source or input changed")
        queue.run_worker(args.campaign, screen, checkpoint_every=args.checkpoint_every)
    else:
        print(json.dumps(verify(args.campaign)))


if __name__ == "__main__":
    main()
