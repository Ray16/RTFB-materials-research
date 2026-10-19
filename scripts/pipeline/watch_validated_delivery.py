#!/usr/bin/env python
"""Refresh audited tables and existing plots/slides in fresh, CPU-only processes."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import signal
import socket
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True)
    ap.add_argument('--max-hours',type=float,default=12);ap.add_argument('--interval',type=int,default=300)
    ap.add_argument('--audit-script',type=Path,help='independent full-value audit before publishing changed results')
    a=ap.parse_args();run=a.run.resolve();out=run/'delivery_watch';out.mkdir(exist_ok=True)
    guard=(out/'watch.lock').open('a+');fcntl.flock(guard,fcntl.LOCK_EX|fcntl.LOCK_NB)
    assert sorted(os.sched_getaffinity(0))==[76]
    status=dict(status='watching',host=socket.gethostname(),pid=os.getpid(),
                process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
                started=time.time(),deadline=time.time()+a.max_hours*3600,interval_seconds=a.interval,
                quantum_jobs_submitted=0,iterations=0)
    path=out/'status.json';child=None;last_published=None;last_signature=None
    def save():
        tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(status,indent=2)+'\n');tmp.replace(path)
    def stop(sig,frame):raise KeyboardInterrupt(f'signal{sig}')
    signal.signal(signal.SIGTERM,stop)
    env=dict(os.environ,PYTHONPATH=str(ROOT/'src'),MPLBACKEND='Agg')
    for name in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):env[name]='1'
    def execute(argv,log,timeout):
        nonlocal child
        with log.open('w') as output:
            child=subprocess.Popen(argv,cwd=ROOT,env=env,stdout=output,stderr=subprocess.STDOUT,start_new_session=True)
            rc=child.wait(timeout=timeout);child=None
        return rc
    save()
    try:
        while time.time()<status['deadline']:
            stamp=str(time.time_ns());status['iterations']+=1;save()
            collect=['taskset','-c','76',sys.executable,str(ROOT/'scripts/pipeline/collect_verified_candidates.py'),
                     '--out',str(run/'verified_collection'),'--queue',str(run)]
            rc=execute(collect,out/('collect_'+stamp+'.log'),180)
            status['last_collector_returncode']=rc
            if rc==0:
                feed=run/'verified_collection/latest.json';payload=json.loads(feed.read_text());s=payload['summary']
                status['last_snapshot']=payload['snapshot_path']
                status['last_admitted_candidates']=s['original39_native_potential_molecules']
                status['last_arithmetic_discrepancies']=s['arithmetic_discrepancies']
                signature=hashlib.sha256(json.dumps(dict(original39=payload['original39'],additional_candidates=payload.get('additional_candidates',payload.get('extra_NDI',[])),
                    inputs=s.get('admitted_input_sha256',{}),sources=s.get('source_sha256',{})),sort_keys=True).encode()).hexdigest()
                # A source-changing snapshot cannot replace the checked delivered figures.
                if (s['arithmetic_discrepancies']==0 and not s['inputs_changed_during_snapshot']
                        and s['original39_native_potential_molecules']>0 and signature!=last_signature):
                    if a.audit_script:
                        arc=execute(['taskset','-c','72,73',sys.executable,str(a.audit_script.resolve())],
                                    out/('audit_'+stamp+'.log'),600)
                        status['last_full_value_audit_returncode']=arc
                        if arc:
                            status.update(updated=time.time(),publication_held_for_full_value_audit=True);save()
                            time.sleep(min(a.interval,max(0,status['deadline']-time.time())))
                            continue
                        status['publication_held_for_full_value_audit']=False
                    publish=['taskset','-c','77',sys.executable,str(ROOT/'scripts/pipeline/update_verified_slides.py'),
                             '--out',str(run/'slides_publisher'),'--feed',str(feed),'--update-plots']
                    prc=execute(publish,out/('publish_'+stamp+'.log'),300)
                    status['last_publisher_returncode']=prc
                    sp=run/'slides_publisher/status.json'
                    if prc==0 and sp.exists() and json.loads(sp.read_text()).get('status')=='published_checked_pdf':
                        last_published=payload['snapshot_path'];last_signature=signature;status['last_published_snapshot']=last_published
            status['updated']=time.time();save()
            time.sleep(min(a.interval,max(0,status['deadline']-time.time())))
        status.update(status='deadline_reached',finished=time.time());save()
    except BaseException as exc:
        status.update(status='stopped_preserving_checked_outputs',error=repr(exc),finished=time.time());save();raise
    finally:
        if child is not None and child.poll() is None:
            os.killpg(child.pid,signal.SIGTERM)
            try:child.wait(timeout=30)
            except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
        guard.close()


if __name__=='__main__':main()
