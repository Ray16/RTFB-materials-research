#!/usr/bin/env python3
"""Collect only exact, checked and model-matched stationary E cycles."""
import argparse
import json
import os
from pathlib import Path

from redox.core.io import atomic_write_json
from redox.workflows.ionic_stationary_cycles import analyze


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--preparation',type=Path,required=True);p.add_argument('--review',type=Path,action='append',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out.exists():raise FileExistsError('Preserve earlier electronic-cycle analyses')
    result=analyze(a.preparation.resolve(),[path.resolve() for path in a.review]);atomic_write_json(a.out,result)
    print(json.dumps({k:result[k] for k in ('status','computed_cycles','native_G_qualified')}))
