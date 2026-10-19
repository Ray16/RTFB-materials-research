"""Check only ready parent states, sequentially under their molecule spin lock."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import socket
import time
import traceback

from redox.core.io import atomic_write_json
from redox.validation import spins
from redox.workflows.resources import verify_gate_owner

ap=argparse.ArgumentParser();ap.add_argument('--parent',type=int,required=True);args=ap.parse_args()
base=Path(__file__).resolve().parent
plan=json.loads((base/'plan.json').read_text());task=plan['parents'][args.parent]
out=base/'parents'/task['id'];out.mkdir(parents=True,exist_ok=False)
report=dict(status='starting',id=task['id'],host=socket.gethostname(),pid=os.getpid(),
            process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
            cpu_affinity=sorted(os.sched_getaffinity(0)),started=time.time(),cases=[],
            experimental_accuracy_validated=False)
def save():atomic_write_json(out/'status.json',report)
save()
try:
    verify_gate_owner();assert len(os.sched_getaffinity(0))==2
    for name in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        assert os.environ[name]=='2'
    for path,digest in plan['producer_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==digest,path
    failures=0
    for case in task['cases']:
        assert spins.snapshot(task['id'],case['state'])==case['chosen_identity'],'Selected source changed'
        report.update(status='computing',state=case['state']);save()
        result=spins.check(task['id'],case['state']);report['cases'].append(result);save()
        if result['status'] in ('failed','source_stale'):
            failures+=1
        elif result['status'] in ('spin_gap_verified','lower_spin_promoted'):
            current=spins.read_comparison(task['id'],case['state'])
            assert current and current['status']==result['status']
            failures=0
        elif result['status'] in ('near_gap_waiting_chosen_native','lower_spin_waiting_motif_release'):
            failures=0
        else:
            raise RuntimeError('No completed or explicitly gated spin result')
        if failures>=3:raise RuntimeError('Three consecutive real spin failures; circuit breaker')
    report.update(status='completed_with_gates' if any(c['status'] not in
        ('spin_gap_verified','lower_spin_promoted') for c in report['cases']) else 'comparisons_verified')
except BaseException as exc:
    report.update(status='failed_preserved',error=repr(exc),traceback=traceback.format_exc());raise
finally:
    report['finished']=time.time();save()
