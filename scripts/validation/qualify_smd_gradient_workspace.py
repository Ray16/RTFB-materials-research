#!/usr/bin/env python
"""Reuse the full workspace spin control with the bounded gradient adapter."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mult',type=int,choices=(1,2),required=True)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--workspace-GiB',type=float,required=True)
    args=parser.parse_args();root=Path(__file__).resolve().parents[2]
    previous=root/'scripts/validation/qualify_smd_workspace.py'
    adapter=root/'scripts/pipeline/smd_gradient_workspace.py'
    spec=importlib.util.spec_from_file_location('preserved_full_workspace_control',previous)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    old=module.adapted_body
    def adapted(body):
        value=old(body);before='scripts/pipeline/smd_workspace_execution.py'
        if value.count(before)!=1:raise ValueError('Workspace adapter template changed')
        return value.replace(before,'scripts/pipeline/smd_gradient_workspace.py')
    module.adapted_body=adapted
    try:module.main()
    finally:
        path=args.out/'status.json'
        if path.exists():
            result=json.loads(path.read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
            result.setdefault('source_sha256',{}).update({str(p):sha(p) for p in (Path(__file__).resolve(),previous,adapter)})
            result['gradient_execution_variant']='Installed bounded contractions; energy direct-LU flag unchanged'
            path.write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
