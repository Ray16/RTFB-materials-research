#!/usr/bin/env python
"""Run the existing full spin control with an explicit execution-memory budget."""
import argparse
import hashlib
import importlib.util
import inspect
import json
import math
from pathlib import Path
import sys
import textwrap


def adapted_body(body):
    changes = (
        ('source = ROOT / "scripts/pipeline/smd_memory_execution.py"',
         'source = ROOT / "scripts/pipeline/smd_workspace_execution.py"'),
        ('parser.add_argument("--mult", type=int, choices=(1, 2), required=True)',
         'parser.add_argument("--mult", type=int, choices=(1, 2), required=True)\n    parser.add_argument("--workspace-GiB", type=float, required=True)'),
        ('memory.compact_smd_execution(arm / "execution")',
         'memory.compact_smd_execution(arm / "execution", budget_GiB=args.workspace_GiB)'),
    )
    for before,after in changes:
        if body.count(before) != 1:
            raise ValueError("Full spin-control template changed")
        body = body.replace(before,after)
    return body


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mult",type=int,choices=(1,2),required=True)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--workspace-GiB",type=float,required=True)
    args=parser.parse_args()
    if not math.isfinite(args.workspace_GiB) or args.workspace_GiB<=0:
        raise ValueError("Positive finite workspace budget required")
    root=Path(__file__).resolve().parents[2]
    original=root/"scripts/validation/qualify_smd_memory_execution.py"
    adapter=root/"scripts/pipeline/smd_workspace_execution.py"
    spec=importlib.util.spec_from_file_location("unchanged_full_workspace_spin_control",original)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    body=adapted_body(inspect.getsource(module.main))
    snapshot=args.out.parent/f"workspace_qualifier_{args.workspace_GiB:g}GiB_mult{args.mult}.py"
    snapshot.parent.mkdir(parents=True,exist_ok=True)
    with snapshot.open("x") as stream:stream.write(body)
    namespace=dict(vars(module),__file__=str(Path(__file__).resolve()))
    exec(compile(textwrap.dedent(body),str(snapshot),"exec"),namespace)
    try:
        namespace["main"]()
    finally:
        path=args.out/"status.json"
        if path.exists():
            value=json.loads(path.read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
            value.setdefault("source_sha256",{}).update({str(p):sha(p) for p in (Path(__file__).resolve(),original,adapter,snapshot)})
            value.update(execution_workspace_GiB=args.workspace_GiB,
                control_body_change="Execution adapter path and workspace argument only; quantum models and acceptance thresholds unchanged")
            path.write_text(json.dumps(value,indent=2)+"\n")


if __name__=="__main__":main()
