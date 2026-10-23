#!/usr/bin/env python
"""Continue an unchanged checked resource profile in fresh geometry-force storage.

The pilot and old FileExists failure remain immutable. Only the full-calculation
storage changes; a copied pilot receipt explicitly identifies its original run.
"""
import argparse
import json
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def verify_stopped(parent):
    groups = {}
    for p in parent.glob("energy/workers/*.json"):
        d = json.loads(p.read_text())
        groups.setdefault(d["host"], []).append(dict(pid=d["pid"], start=d["process_start"]))
    if not groups:
        raise ValueError("An actual parent execution and its ownership records are required")
    code = '''import json,os
from pathlib import Path
c=json.load(__import__('sys').stdin);rows=[]
for r in c['owners']:
 p=Path('/proc')/str(r['pid'])/'stat'
 live=p.exists() and p.read_text().rsplit(')',1)[1].split()[19]==r['start']
 rows.append(dict(r,same_process_live=live))
matches=[]
for p in Path('/proc').glob('[0-9]*/cmdline'):
 try:command=p.read_bytes().replace(bytes([0]),b' ').decode()
 except (OSError,UnicodeError):continue
 if c['parent'] in command:matches.append(dict(pid=int(p.parent.name),command=command))
print(json.dumps(dict(host=os.uname().nodename,owners=rows,matching_processes=matches)))'''
    records = []
    for host, owners in groups.items():
        r = subprocess.run(["ssh", "-o", "BatchMode=yes", "-o", "ConnectTimeout=8", host,
            shlex.join([sys.executable, "-c", code])], input=json.dumps(dict(parent=str(parent), owners=owners)),
            text=True, capture_output=True, timeout=20, check=True)
        d = json.loads(r.stdout)
        if any(row["same_process_live"] for row in d["owners"]) or d["matching_processes"]:
            raise RuntimeError("Parent owner or predictor is still live; preserve its workspace")
        records.append(d)
    return records


def prepare(parent, out):
    old = api.read_plan(parent)
    api.verify_sources(old["producer_sha256"])
    api.verify_sources(old["input_sha256"])
    if len(old["tasks"]) != 1:
        raise ValueError("One explicitly declared baseline task is required")
    key = old["tasks"][0]["key"]
    failure_path = parent / "energy/failures" / (key + ".json")
    failure = json.loads(failure_path.read_text())
    pilot_path = parent / "energy/pilot/status.json"
    pilot = json.loads(pilot_path.read_text())
    forces = Path(old["raw_root"]) / "energy" / (key + "_CDS.evaluations")
    if (not failure["error"].startswith("FileExistsError(") or not forces.is_dir()
            or (parent / "energy/tasks" / key).exists()
            or pilot["status"] != "memory_pilot_passed"
            or not 0 < pilot["observed_peak_GPU_MiB"] <= old["maximum_memory_pilot_MiB"]):
        raise ValueError("The checked-pilot force-storage collision is not established")
    api.verify_sources(pilot["output_sha256"])
    ownership = verify_stopped(parent)
    out.mkdir(parents=True, exist_ok=False)
    shutil.copytree(parent / "live_inputs", out / "live_inputs")
    shutil.copytree(parent / "controls", out / "controls")
    for stage in ("energy", "thermal"):
        for name in ("tasks", "done", "claims", "failures", "workers", "pilot", "launches"):
            (out / stage / name).mkdir(parents=True)
    inputs = dict(old["input_sha256"])
    for p in (parent / "plan.json", pilot_path, failure_path):
        inputs[str(p)] = api.digest(p)
    for task in old["tasks"]:
        p = out / "live_inputs" / task["key"] / "task.json"
        ready = json.loads(p.read_text())
        ready["input_xyz"] = str(p.parent / "input.xyz")
        ready["input_sha256"][ready["input_xyz"]] = api.digest(ready["input_xyz"])
        atomic_write_json(p, ready)
        inputs[str(p)] = api.digest(p)
        inputs.update(ready["input_sha256"])
    copied = dict(pilot, resource_profile_reuse=dict(original_status=str(pilot_path),
        original_sha256=api.digest(pilot_path), new_calculation_performed=False,
        scope="Identical task, numerical model and execution adapter; fresh full-geometry storage"))
    atomic_write_json(out / "energy/pilot/status.json", copied)
    energy = parent / "energy/pilot/energy_and_gradient.json"
    shutil.copyfile(energy, out / "energy/pilot/energy_and_gradient.json")
    inputs[str(out / "energy/pilot/status.json")] = api.digest(out / "energy/pilot/status.json")
    inputs[str(out / "energy/pilot/energy_and_gradient.json")] = api.digest(energy)
    sources = dict(old["producer_sha256"])
    sources[str(Path(__file__).resolve())] = api.digest(__file__)
    new = dict(old, producer_sha256=sources, input_sha256=inputs,
        raw_root=str(ROOT / "calcs/polymer_models" / out.parent.name / out.name),
        parent_campaign=str(parent), recovery_scope="Fresh full-calculation geometry-force storage after a checked resource pilot; model/tasks unchanged",
        created=time.time())
    atomic_write_json(out / "plan.json", new)
    atomic_write_json(out / "parent_ownership_verification.json", dict(checked=time.time(), hosts=ownership))
    print(json.dumps(dict(status="postpilot_continuation_prepared", task=key,
        original_resource_pilot_reused=True, prior_failure_preserved=True)))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--parent", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    args = p.parse_args()
    from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    prepare(args.parent.resolve(), args.out.resolve())


if __name__ == "__main__":
    main()
