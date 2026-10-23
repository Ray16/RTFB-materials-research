#!/usr/bin/env python
"""Preserve a failed resource pilot and freeze a control-qualified execution branch."""
import argparse
import json
from pathlib import Path
import shutil
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--controls", type=Path, required=True)
    args = parser.parse_args()
    parent, out = args.parent.resolve(), args.out.resolve()
    old = api.read_plan(parent)
    failure_path = parent / "energy/pilot/status.json"
    failure = json.loads(failure_path.read_text())
    if failure["status"] != "failed_preserved" or "out of memory" not in failure.get("error", ""):
        raise ValueError("This branch requires the reviewed actual resource failure")
    controls = [args.controls.resolve() / "controls" / str(mult) / "status.json" for mult in (1, 2)]
    recovery = api.load_module("qualified_SMD_memory_recovery", ROOT / "scripts/pipeline/run_smd_memory_recovery.py")
    recovery.checked_controls(dict(SMD_memory_qualification_controls=list(map(str, controls))))
    out.mkdir(parents=True, exist_ok=False)
    shutil.copytree(parent / "live_inputs", out / "live_inputs")
    shutil.copytree(parent / "controls", out / "controls")
    for stage in ("energy", "thermal"):
        for name in ("tasks", "done", "claims", "failures", "workers", "pilot", "launches"):
            (out / stage / name).mkdir(parents=True)
    inputs = {str(parent / "plan.json"): api.digest(parent / "plan.json"), str(failure_path): api.digest(failure_path)}
    for task in old["tasks"]:
        path = out / "live_inputs" / task["key"] / "task.json"
        ready = json.loads(path.read_text())
        ready["input_xyz"] = str(path.parent / "input.xyz")
        ready["input_sha256"][str(path.parent / "input.xyz")] = api.digest(path.parent / "input.xyz")
        atomic_write_json(path, ready)
        inputs[str(path)] = api.digest(path)
        inputs.update(ready["input_sha256"])
    plan = dict(old, raw_root=str(ROOT / "calcs/polymer_models" / out.parent.name / out.name),
        source_failed_resource_pilot=str(failure_path), SMD_memory_qualification_controls=list(map(str, controls)),
        execution_driver=str(ROOT / "scripts/pipeline/run_smd_memory_recovery.py"),
        input_sha256=inputs, producer_sha256=dict(old["producer_sha256"]), recovery_created=time.time(),
        recovery_scope="Same chemical/numerical model and starting geometry; temporary-buffer lifetime reduction only")
    for path in [Path(__file__).resolve(), ROOT / "scripts/pipeline/run_smd_memory_recovery.py", ROOT / "scripts/pipeline/smd_memory_execution.py", *controls]:
        plan["producer_sha256"][str(path)] = api.digest(path)
    api.verify_sources(plan["producer_sha256"])
    api.verify_sources(inputs)
    atomic_write_json(out / "plan.json", plan)
    print(json.dumps(dict(status="qualified_execution_recovery_prepared", tasks=len(plan["tasks"]), unchanged_memory_safety_budget_MiB=plan["maximum_memory_pilot_MiB"])))


if __name__ == "__main__":
    main()
