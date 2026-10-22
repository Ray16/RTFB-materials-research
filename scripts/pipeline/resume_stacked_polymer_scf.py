#!/usr/bin/env python
"""Continue the failed stacked PMDI trimer from its last checked geometry."""

import argparse
import json
from pathlib import Path
import shutil

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

SCRIPT = Path(__file__).resolve()
MODULE = ROOT / "src/redox/qm/stable_geometry_scf.py"
PARENT = ROOT / "results/polymer_models/campaigns/2026-10-21_checkpoint_pmdi_stacked"
CAMPAIGN = (
    ROOT / "results/polymer_models/campaigns/2026-10-21_checkpoint_pmdi_stacked_scf"
)
KEY = "design_pmdi_ethyl__n3__stereo0__stacked"


def prepare():
    old = api.read_plan(PARENT)
    api.verify_sources(old["producer_sha256"])
    failure = PARENT / "failures" / f"{KEY}.json"
    issue = json.loads(failure.read_text())
    ownership = api.load_module(
        "scf_owner", ROOT / "scripts/pipeline/run_polymer_spacing.py"
    )
    if ownership.owner_alive(issue["owner"]) or (PARENT / "claims" / KEY).exists():
        raise RuntimeError("Prior worker still owns the task")
    if "Nuclear gradients" not in issue["error"]:
        raise ValueError("This continuation requires the recorded SCF failure")
    ready_path = PARENT / "live_inputs" / KEY / "task.json"
    task = json.loads(ready_path.read_text())
    source = (
        Path(old["raw_root"]) / KEY / "rest/geometry_checkpoints/last_evaluated.xyz"
    )
    review = api.load_module(
        "scf_geometry_audit", ROOT / "scripts/validation/audit_polymer_geometries.py"
    ).audit(task, source)
    if review["issues"]:
        raise ValueError("Checkpoint fails geometry checks")
    CAMPAIGN.mkdir()
    for name in (
        "claims",
        "tasks",
        "done",
        "failures",
        "workers",
        "launches",
        "pilot",
        "live_inputs",
    ):
        (CAMPAIGN / name).mkdir()
    directory = CAMPAIGN / "live_inputs" / KEY
    directory.mkdir()
    shutil.copyfile(source, directory / "input.xyz")
    task.update(
        input_xyz=str(directory / "input.xyz"),
        input_sha256={
            str(source): api.digest(source),
            str(ready_path): api.digest(ready_path),
            str(directory / "input.xyz"): api.digest(directory / "input.xyz"),
        },
    )
    atomic_write_json(directory / "task.json", task)
    atomic_write_json(CAMPAIGN / "checkpoint_geometry_review.json", review)
    shutil.copyfile(PARENT / "pilot/status.json", CAMPAIGN / "pilot/status.json")
    plan = dict(
        old,
        raw_root=str(ROOT / "calcs/polymer_models" / CAMPAIGN.name),
        recovery_reason=(
            "Scanner-safe SCF stabilization at unchanged tolerances; "
            "preserve failed parent"
        ),
    )
    plan["producer_sha256"].update(
        {str(SCRIPT): api.digest(SCRIPT), str(MODULE): api.digest(MODULE)}
    )
    inputs = [
        failure,
        ready_path,
        source,
        directory / "task.json",
        CAMPAIGN / "checkpoint_geometry_review.json",
        CAMPAIGN / "pilot/status.json",
        PARENT / "plan.json",
    ]
    plan["input_sha256"] = {str(path): api.digest(path) for path in inputs}
    atomic_write_json(CAMPAIGN / "plan.json", plan)
    pointer = ROOT / "results/polymer_models/recovery_campaigns.json"
    registry = json.loads(pointer.read_text())
    registry["campaigns"].append(
        dict(
            path=str(CAMPAIGN.relative_to(ROOT)),
            plan_sha256=api.digest(CAMPAIGN / "plan.json"),
        )
    )
    atomic_write_json(pointer, registry)


def runner():
    from redox.qm.stable_geometry_scf import stabilize
    import redox.qm.dft as dft

    continuation = api.load_module(
        "stacked_scf_continuation", ROOT / "scripts/pipeline/continue_polymer_models.py"
    )
    continuation.RECOVERY, continuation.SCRIPT = CAMPAIGN, SCRIPT
    base = continuation.recovery_runner()
    original = dft._build_smd_mf

    def build(*args, **kwargs):
        mean_field = original(*args, **kwargs)
        if (args[1] if len(args) > 1 else kwargs.get("xc")) == "r2scan":
            return stabilize(
                mean_field,
                Path(api.read_plan(CAMPAIGN)["raw_root"])
                / KEY
                / "rest/scf_recovery.jsonl",
            )
        return mean_field

    dft._build_smd_mf = build
    return base


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "launch", "worker"))
    parser.add_argument("--role", default="queue_large", choices=("queue_large",))
    args = parser.parse_args()
    if args.action == "prepare":
        prepare()
    elif args.action == "launch":
        from redox.workflows import resources

        original = resources.free_slots
        resources.free_slots = lambda *a, **kw: original(*a, **kw)[:1]
        runner().launch(args.role)
    else:
        runner().worker(args.role)
