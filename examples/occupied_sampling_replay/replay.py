#!/usr/bin/env python3
"""Replay the frozen sampling diagnostic from shipped inputs in this checkout."""
import json
import os
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','BLIS_NUM_THREADS'):
    os.environ[key]='1'
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
with allocation_guard():
    host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
from redox.validation.sampling_publication import audit
print(json.dumps(audit(Path(__file__).parent/'data')))
