#!/usr/bin/env python
"""Reserved resumable validation queue; larger budgets wait without occupying GPUs."""
import argparse
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import socket
import subprocess
import sys
import time

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'src'))
from redox.validation.spins import atomic_json
GATE='/nfs/lambda_stor_01/homes/rzhu/bin/gpu_reserve'


def save(path,data):atomic_json(path,data)


def archive_partial(run,job):
    from redox.core.common import DFT
    from redox.qm.native_thermal import cache_path,load_native_thermal
    archived=[]
    with (ROOT/'logs/pipeline_locks'/f"{job['gid']}.lock").open('a+') as ml, \
         (ROOT/'logs/pipeline_locks'/f"{job['gid']}_native_thermal.lock").open('a+') as tl:
        fcntl.flock(ml,fcntl.LOCK_EX|fcntl.LOCK_NB);fcntl.flock(tl,fcntl.LOCK_EX|fcntl.LOCK_NB)
        for st in job['states']:
            sd=DFT/job['gid']/st;raw=json.loads((sd/'result.json').read_text())
            if load_native_thermal(sd,raw['charge'],raw['mult']):continue
            path,_=cache_path(sd,raw['charge'],raw['mult'])
            if path.exists():raise RuntimeError('Invalid complete native record needs separate review')
            partial=path.with_suffix('')
            if not partial.exists():continue
            target=run/'archives'/job['name']/st/partial.name;target.parent.mkdir(parents=True,exist_ok=True)
            sums={str(p.relative_to(partial)):hashlib.sha256(p.read_bytes()).hexdigest() for p in partial.rglob('*') if p.is_file()}
            if target.exists():raise RuntimeError('Archive collision; prior evidence retained')
            shutil.move(partial,target)
            assert sums=={str(p.relative_to(target)):hashlib.sha256(p.read_bytes()).hexdigest() for p in target.rglob('*') if p.is_file()}
            archived.append(dict(source=str(partial),archive=str(target),sha256=sums,
                                 reason=job['reason'],owner_audit=job['owner_audit']))
    save(run/'archives'/job['name']/'inventory.json',archived)


def verify(job,run):
    if job['kind']=='control':
        r=json.loads((run/'memory_controls'/(job['state']+'_control.json')).read_text())
        if r['status']!='control_verified':raise RuntimeError('Control did not qualify')
        with (run/'qualification.lock').open('a+') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            if all((run/'memory_controls'/(s+'_control.json')).exists() for s in ('neu','red1')):
                from redox.qm.hessian_memory import CERTIFICATE
                if not CERTIFICATE.exists():
                    rc=subprocess.call([sys.executable,str(ROOT/'scripts/validation/qualify_hessian_memory.py'),
                                        '--out',str(run/'memory_controls'),'--seal'],cwd=ROOT)
                    if rc:raise RuntimeError('Memory qualification certificate failed')
        return r
    if job['kind']=='spin':
        from redox.validation.spins import read_comparison
        r=read_comparison(job['gid'],job['states'][0])
        if not r or r['status'] not in ('spin_gap_verified','lower_spin_promoted','near_gap_waiting_chosen_native','lower_spin_waiting_motif_release'):
            raise RuntimeError('No current source-addressed completed spin comparison')
        return {k:v for k,v in r.items() if k in ['status','electronic_gap_eV','thermochemical_gap_eV','ensemble_primary_integration_pending']}
    from redox.core.common import DFT
    from redox.core.protocol import load_record
    from redox.qm.native_thermal import load_native_thermal
    results=[]
    for st in job['states']:
        sd=DFT/job['gid']/st;r=json.loads((sd/'result.json').read_text())
        th=load_native_thermal(sd,r['charge'],r['mult'])
        if not th or not load_record(sd,r['charge'],r['mult']):raise RuntimeError('Native/SP coverage failed after zero exit')
        results.append(dict(state=st,thermal=th['g_thermal_total_eV'],execution_memory=th.get('execution_memory')))
    return results


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);ap.add_argument('--worker',required=True);a=ap.parse_args()
    run=a.run.resolve();q=json.loads((run/'queue.json').read_text());allocation=next(w for w in q['workers'] if w['name']==a.worker)
    assert socket.gethostname().split('.')[0]==allocation['host'] and sorted(os.sched_getaffinity(0))==[allocation['controller_cpu']]
    for path,expected in q['source_sha256'].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest()==expected,'Prepared source changed'
    directory=run/'workers'/a.worker;directory.mkdir(exist_ok=True)
    guard=(directory/'worker.lock').open('a+');fcntl.flock(guard,fcntl.LOCK_EX|fcntl.LOCK_NB)
    status=dict(status='running',pid=os.getpid(),host=socket.gethostname(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
                allocation=allocation,started=time.time(),jobs=[])
    def update():save(directory/'status.json',status)
    def stop(sig,frame):raise KeyboardInterrupt(f'signal {sig}')
    signal.signal(signal.SIGTERM,stop)
    child=None;claim=None;failures=0;update()
    env=dict(os.environ,PYTHONPATH=str(ROOT/'src'))
    env['PATH']=str(Path(sys.executable).parent)+':'+env.get('PATH','')
    env['LD_LIBRARY_PATH']=(run/'library_path.txt').read_text().strip()+':'+env.get('LD_LIBRARY_PATH','')
    for name in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):env[name]='2'
    try:
        while time.time()<q['deadline']:
            selected=None;pending=False
            for job in q['jobs']:
                jp=run/'job_status'/(job['name']+'.json');prior=json.loads(jp.read_text()) if jp.exists() else {}
                if prior.get('status') in ('verified','failed_preserved','dependency_failed'):continue
                pending=True
                if job.get('memory_cap_GiB')==24:
                    from redox.qm.hessian_memory import CERTIFICATE
                    if not CERTIFICATE.exists():
                        control_status=[json.loads(p.read_text()).get('status') for p in (run/'job_status').glob('control_*.json')]
                        if 'failed_preserved' in control_status:
                            save(jp,dict(status='dependency_failed',reason='Execution budget controls failed; no GPU claimed'))
                        continue
                proposed=run/'claims'/job['name']
                try:proposed.mkdir()
                except FileExistsError:continue
                claim=proposed;selected=job;save(claim/'owner.json',dict(host=status['host'],pid=os.getpid(),process_start=status['process_start'],worker=a.worker));break
            if selected is None:
                if pending:time.sleep(15);continue
                break
            job=selected;jp=run/'job_status'/(job['name']+'.json');entry=dict(job=job,status='preparing',worker=a.worker,started=time.time());status['jobs'].append(entry);save(jp,entry);update()
            try:
                if job['kind']=='native':archive_partial(run,job)
                taskenv=dict(env)
                taskenv['REDOX_HESSIAN_MEMORY_GIB']=str(job.get('memory_cap_GiB',6))
                if job['kind']=='control':
                    argv=[sys.executable,'-u',str(ROOT/'scripts/validation/qualify_hessian_memory.py'),'--out',str(run/'memory_controls'),'--state',job['state']]
                elif job['kind']=='spin':
                    argv=[sys.executable,'-u',str(ROOT/'scripts/validation/queued_candidate_spin_checks.py'),'--worker','--only',job['gid']+':'+job['states'][0],'--cpus',','.join(map(str,allocation['compute_cpus']))]
                else:
                    argv=[sys.executable,'-u','-m','redox.qm.native_thermal','--only',','.join(job['gid']+':'+s for s in job['states']),
                          '--threads','2','--select-policy','--repair-minimum']
                argv=[GATE,'run',str(allocation['gpu']),'--label','redox-resume-'+job['name'],'--','taskset','-c',','.join(map(str,allocation['compute_cpus'])),
                      'timeout','--signal=TERM','--kill-after=30','6h',*argv]
                entry.update(status='computing',argv=argv);save(jp,entry);update();contention=0
                while True:
                    log=run/'logs'/(job['name']+f'_{contention:02d}.log');entry['log']=str(log);save(jp,entry)
                    with log.open('w') as output:
                        child=subprocess.Popen(argv,cwd=ROOT,env=taskenv,stdout=output,stderr=subprocess.STDOUT,start_new_session=True)
                        entry['gate_pid']=child.pid;save(jp,entry);rc=child.wait();child=None
                    if rc and 'BUSY/reserved -> not running' in log.read_text():
                        if time.time()>=q['deadline']:raise RuntimeError('Contention deadline; no scientific attempt consumed')
                        entry['resource_contentions']=contention+1;save(jp,entry);contention+=1;time.sleep(15);continue
                    break
                if rc:raise RuntimeError(f'Validation failed rc={rc}; all outputs preserved')
                evidence=verify(job,run)
                entry.update(status='verified',evidence=evidence,finished=time.time(),experimental_accuracy_validated=False);failures=0
            except Exception as exc:
                entry.update(status='failed_preserved',error=repr(exc),finished=time.time());failures+=1
            finally:
                save(jp,entry);shutil.rmtree(claim);claim=None;update()
            if failures>=3:raise RuntimeError('Three consecutive real validation failures; worker circuit breaker')
        status.update(status='finished',finished=time.time());update()
    except BaseException as exc:
        status.update(status='stopped_preserving_outputs',error=repr(exc),finished=time.time());update();raise
    finally:
        if child is not None and child.poll() is None:
            os.killpg(child.pid,signal.SIGTERM)
            try:child.wait(timeout=30)
            except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
        if claim is not None:shutil.rmtree(claim)
        guard.close()


if __name__=='__main__':main()
