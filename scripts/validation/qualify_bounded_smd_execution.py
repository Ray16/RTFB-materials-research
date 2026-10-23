#!/usr/bin/env python
"""Reuse the established two-spin energy/geometry/gradient/Hessian control.

The only change to the control body selects the new run-local execution adapter.
All generated and original producer bytes are bound into the final control proof.
"""
import argparse
import hashlib
import importlib.util
import inspect
import json
from pathlib import Path
import textwrap


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mult",type=int,choices=(1,2),required=True)
    parser.add_argument("--out",type=Path,required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[2]
    original = root / "scripts/validation/qualify_smd_memory_execution.py"
    adapter = root / "scripts/pipeline/smd_bounded_execution.py"
    spec = importlib.util.spec_from_file_location("unchanged_full_spin_control",original)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    before = inspect.getsource(module.main)
    old = 'source = ROOT / "scripts/pipeline/smd_memory_execution.py"'
    if before.count(old) != 1:
        raise ValueError("Original control body differs from the qualified source substitution")
    after = before.replace(old,'source = ROOT / "scripts/pipeline/smd_bounded_execution.py"')
    snapshot = args.out.parent / ("bounded_qualifier_mult"+str(args.mult)+".py")
    snapshot.parent.mkdir(parents=True,exist_ok=True)
    with snapshot.open("x") as stream:
        stream.write(after)
    namespace = dict(vars(module),__file__=str(Path(__file__).resolve()))
    exec(compile(textwrap.dedent(after),str(snapshot),"exec"),namespace)
    try:
        namespace["main"]()
    finally:
        status = args.out / "status.json"
        if status.exists():
            value = json.loads(status.read_text())
            sha = lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
            value.setdefault("source_sha256",{}).update({str(p):sha(p) for p in (Path(__file__).resolve(),original,adapter,snapshot)})
            value["control_body_change"] = "Adapter path only; quantum models and acceptance thresholds unchanged"
            status.write_text(json.dumps(value,indent=2)+"\n")


if __name__ == "__main__":
    main()
