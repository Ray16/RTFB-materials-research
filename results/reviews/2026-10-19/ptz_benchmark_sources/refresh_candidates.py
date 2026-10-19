"""CPU-only snapshot refresh with a disjoint allocation and capped native pools."""
import os
from pathlib import Path
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity

for key in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS',
            'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    assert os.environ.get(key) == '2'
out = ROOT / 'results/reviews/2026-10-19/post_parent_registration_candidates'
assert not out.exists(), 'Keep existing snapshots intact'
with allocation_guard():
    capacity = host_capacity(2)
    cpus = cpu_allocation(2)
    os.sched_setaffinity(0, cpus)
    atomic_write_json(Path(__file__).with_name('candidate_refresh_allocation.json'),
                      dict(pid=os.getpid(), cpu_affinity=cpus, host_preflight=capacity,
                           native_threads=2, started=time.time()))
os.execv(sys.executable, [sys.executable, str(ROOT / 'scripts/mining/collect_proposal_coverage.py'),
                          '--out', str(out)])
