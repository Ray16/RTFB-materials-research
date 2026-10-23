#!/usr/bin/env python
"""Freeze a new branch after an observed OOM or unchanged-budget pilot failure.

Reuse the existing recovery preparation body with guarded changes to execution
paths and the accepted resource-failure category. Numerical task/model fields
remain byte-for-byte equal to the declared parent plan.
"""
import argparse
import inspect
from pathlib import Path
import textwrap

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--controls",type=Path,required=True)
    args = parser.parse_args()
    source = ROOT / "scripts/pipeline/prepare_smd_memory_recovery.py"
    old = api.load_module("preserved_resource_recovery_preparation",source)
    body = inspect.getsource(old.main)
    guard = 'if failure["status"] != "failed_preserved" or "out of memory" not in failure.get("error", ""):'
    if body.count(guard) != 1:
        raise ValueError("Recovery preparation failure guard changed")
    body = body.replace(guard,'if failure["status"] != "failed_preserved" or not any(token in failure.get("error", "") for token in ("out of memory", "memory pilot exceeds admitted safety budget")):')
    for before,after in (("scripts/pipeline/run_smd_memory_recovery.py","scripts/pipeline/run_bounded_smd_recovery.py"),
                         ("scripts/pipeline/smd_memory_execution.py","scripts/pipeline/smd_bounded_execution.py")):
        if before not in body:
            raise ValueError("Recovery execution-path template changed")
        body = body.replace(before,after)
    body = body.replace("temporary-buffer lifetime reduction only","bounded workspace and equivalent solvent matrix-vector contraction")
    namespace = dict(vars(old),__file__=str(Path(__file__).resolve()))
    exec(compile(textwrap.dedent(body),str(Path(__file__).resolve())+":adapted_preparation","exec"),namespace)
    namespace["main"]()
    out = args.out.resolve()
    snapshot = out / "preparation_execution_body.py"
    snapshot.write_text(body)
    plan = api.read_plan(out)
    for path in (source, ROOT / "scripts/pipeline/run_smd_memory_recovery.py",ROOT / "scripts/pipeline/smd_memory_execution.py",snapshot):
        plan["producer_sha256"][str(path)] = api.digest(path)
    runner = api.load_module("bounded_recovery_control_witnesses",ROOT / "scripts/pipeline/run_bounded_smd_recovery.py")
    plan["producer_sha256"].update(runner.checked_controls(plan))
    atomic_write_json(out / "plan.json",plan)


if __name__ == "__main__":
    main()
