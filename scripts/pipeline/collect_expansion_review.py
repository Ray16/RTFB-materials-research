#!/usr/bin/env python
"""Publish the systematic design/evaluation loop from existing checked outputs.

This read-only collector neither changes registries nor submits calculations.
Every publication snapshots its exact inputs and refuses concurrent input changes.
"""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json, atomic_text_writer
from redox.screening.expansion_review import build_review

INPUTS = (
    "pool",
    "feed",
    "route_evidence",
    "selection_assessment",
    "grafting",
    "base_assisted",
    "polymer",
    "spacing",
)


def digest(content):
    return hashlib.sha256(content).hexdigest()


def write_csv(path, rows):
    if not rows:
        return
    with atomic_text_writer(path, newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    key: json.dumps(value, sort_keys=True)
                    if isinstance(value, (list, dict))
                    else value
                    for key, value in row.items()
                }
            )


def collect(plan_path, out):
    import fcntl
    from redox.workflows.resources import (
        allocation_guard,
        host_capacity,
        cpu_allocation,
    )

    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    out.mkdir(parents=True, exist_ok=True)
    with (out / "collector.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        content = plan_path.read_bytes()
        plan = json.loads(content)
        paths = {key: (ROOT / plan[key]).resolve() for key in INPUTS}
        blobs = {key: path.read_bytes() for key, path in paths.items()}
        values = {key: json.loads(blob) for key, blob in blobs.items()}
        route_path = ROOT / "scripts/pipeline/screen_grafting_routes.py"
        producers = [
            Path(__file__).resolve(),
            route_path,
            ROOT / "src/redox/screening/expansion_review.py",
            ROOT / "src/redox/screening/comparisons.py",
            ROOT / "src/redox/build/grafting.py",
            ROOT / "scripts/plotting/candidate_expansion_report.py",
        ]
        producer_hashes = {str(path): digest(path.read_bytes()) for path in producers}
        spec = importlib.util.spec_from_file_location(
            "expansion_route_gate", route_path
        )
        gate = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(gate)
        evidence = values["route_evidence"].get("candidates", {})
        assessment = values["selection_assessment"]
        review = build_review(
            values["pool"],
            values["feed"],
            values["grafting"],
            values["base_assisted"],
            assessment,
            values["polymer"],
            plan,
            lambda row: gate.screen(row, evidence, assessment),
        )
        # A live polymer stage is context only. Its numerical fields are never
        # copied into capped-site properties or assigned to uncomputed motifs.
        review["spacing_campaign_status"] = values["spacing"]["status"]
        now = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S_%f")
        snapshot = out / "snapshots" / now
        source_dir = snapshot / "sources"
        source_dir.mkdir(parents=True)
        for key, blob in blobs.items():
            (source_dir / f"{key}.json").write_bytes(blob)
        (source_dir / "plan.json").write_bytes(content)
        if plan_path.read_bytes() != content or any(
            path.read_bytes() != blobs[key] for key, path in paths.items()
        ):
            raise RuntimeError(
                "Inputs changed during collection; archived inputs retained, "
                "publication refused"
            )
        if any(
            digest(Path(path).read_bytes()) != expected
            for path, expected in producer_hashes.items()
        ):
            raise RuntimeError(
                "Producer changed during collection; publication refused"
            )
        review["provenance"] = dict(
            created_UTC=now,
            snapshot=str(snapshot),
            source_candidate_snapshot=values["feed"].get("snapshot_path"),
            input_sha256={str(paths[key]): digest(blob) for key, blob in blobs.items()},
            frozen_input_sha256={
                str(source_dir / f"{key}.json"): digest(blob)
                for key, blob in blobs.items()
            },
            plan_sha256=digest(content),
            producer_sha256=producer_hashes,
            raw_roots=[
                str(ROOT / "calcs/dft"),
                str(ROOT / "calcs/grafting"),
                str(ROOT / "calcs/polymer_models"),
            ],
            input_stability="checked_before_publication",
            calculations_submitted=0,
        )
        atomic_write_json(snapshot / "review.json", review)
        for filename, rows in (
            ("designs.csv", review["records"]),
            ("comparisons.csv", review["comparisons"]),
            ("family_coverage.csv", review["family_coverage"]),
            ("design_cells.csv", review["design_cells"]),
        ):
            write_csv(snapshot / filename, rows)
        atomic_write_json(out / "current.json", review)
        atomic_write_json(
            out / "status.json",
            dict(
                status="published_existing_evaluation_snapshot",
                counts=review["counts"],
                source_candidate_as_of=review["source_candidate_as_of"],
                snapshot=str(snapshot),
                snapshot_sha256=digest((snapshot / "review.json").read_bytes()),
                calculations_submitted=0,
            ),
        )
        for filename in (
            "designs.csv",
            "comparisons.csv",
            "family_coverage.csv",
            "design_cells.csv",
        ):
            shutil.copyfile(snapshot / filename, out / filename)
        subprocess.run(
            [
                str(ROOT / ".venv-figures/bin/python"),
                str(ROOT / "scripts/plotting/candidate_expansion_report.py"),
                "--out",
                str(out),
            ],
            check=True,
        )
        print(json.dumps(review["counts"], indent=2))


def watch(plan_path, out, interval, hours):
    """Bounded downstream refresh only when one of the evaluator outputs changes."""
    import fcntl
    import socket

    out.mkdir(parents=True, exist_ok=True)
    with (out / "watch.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        owner = dict(
            host=socket.gethostname(),
            pid=os.getpid(),
            process_start=Path("/proc/self/stat")
            .read_text()
            .rsplit(")", 1)[1]
            .split()[19],
            deadline=time.time() + hours * 3600,
            publications=0,
        )
        previous = None
        frozen_producers = {
            str(path): digest(path.read_bytes())
            for path in (
                Path(__file__).resolve(),
                ROOT / "src/redox/screening/expansion_review.py",
                ROOT / "src/redox/screening/comparisons.py",
                ROOT / "src/redox/build/grafting.py",
                ROOT / "scripts/pipeline/screen_grafting_routes.py",
                ROOT / "scripts/plotting/candidate_expansion_report.py",
            )
        }
        try:
            while time.time() < owner["deadline"]:
                if any(
                    digest(Path(path).read_bytes()) != expected
                    for path, expected in frozen_producers.items()
                ):
                    raise RuntimeError(
                        "Watcher producer changed; restart with reviewed code"
                    )
                content = plan_path.read_bytes()
                plan = json.loads(content)
                fingerprint = (
                    digest(content),
                    tuple(digest((ROOT / plan[key]).read_bytes()) for key in INPUTS),
                )
                if fingerprint != previous:
                    collect(plan_path, out)
                    previous = fingerprint
                    owner["publications"] += 1
                owner["cpus"] = sorted(os.sched_getaffinity(0))
                owner["native_threads"] = 1
                atomic_write_json(
                    out / "watch_status.json",
                    dict(owner, status="watching_evaluation_outputs"),
                )
                time.sleep(min(interval, max(0, owner["deadline"] - time.time())))
            atomic_write_json(
                out / "watch_status.json", dict(owner, status="bounded_watch_complete")
            )
        except BaseException as error:
            atomic_write_json(
                out / "watch_status.json",
                dict(owner, status="stopped_preserving_outputs", error=repr(error)),
            )
            raise


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--plan", type=Path, default=ROOT / "config/candidate_expansion_plan.json"
    )
    parser.add_argument(
        "--out", type=Path, default=ROOT / "results/discovery/expansion"
    )
    parser.add_argument("--watch", action="store_true")
    parser.add_argument("--interval-seconds", type=float, default=120)
    parser.add_argument("--max-hours", type=float, default=24)
    args = parser.parse_args()
    if args.interval_seconds < 30 or args.max_hours <= 0:
        parser.error("Use an interval of at least 30 seconds and positive max-hours")
    if args.watch:
        watch(
            args.plan.resolve(),
            args.out.resolve(),
            args.interval_seconds,
            args.max_hours,
        )
    else:
        collect(args.plan.resolve(), args.out.resolve())
