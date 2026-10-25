#!/usr/bin/env python
"""Measure GPU-library import residency only inside an actual reservation."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

from redox.workflows.resources import verify_gate_owner


def residents():
    result=subprocess.run(['nvidia-smi','--query-compute-apps=pid,gpu_uuid,used_memory',
        '--format=csv,noheader,nounits'],capture_output=True,text=True,check=True)
    return [line for line in result.stdout.splitlines() if line.split(',')[0].strip()==str(os.getpid())]


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args();verify_gate_owner()
    if len(os.sched_getaffinity(0))!=1:raise ValueError('Exactly one allocated CPU is required')
    before=residents()
    import gpu4pyscf
    after=residents()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    with args.out.open('x') as stream:
        json.dump(dict(status='reserved_GPU_library_import_residency_measured',checked=time.time(),
            host=os.uname().nodename,pid=os.getpid(),CPU_affinity=sorted(os.sched_getaffinity(0)),
            CUDA_VISIBLE_DEVICES=os.environ['CUDA_VISIBLE_DEVICES'],GPU4PySCF_version=gpu4pyscf.__version__,
            own_residency_before_import=before,own_residency_after_import=after,
            explicit_quantum_calculation_run=False),stream,indent=2);stream.write('\n')
    print(json.dumps(dict(before=len(before),after=len(after),pid=os.getpid())))


if __name__=='__main__':main()
