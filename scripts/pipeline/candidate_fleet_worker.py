#!/usr/bin/env python
"""Resumable candidate queue; compute stages use fresh atomic GPU reservations."""
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
sys.path.insert(0,str(ROOT/'src'))
from redox.core.common import read_manifest


def save(path,data):
    tmp=path.with_name(path.name+f'.{os.getpid()}.tmp')
    tmp.write_text(json.dumps(data,indent=2)+'\n');tmp.replace(path)


def process_start():
    return Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]


def run(argv,log,env,timeout):
    with log.open('w') as stream:
        child=subprocess.Popen(argv,cwd=ROOT,env=env,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
        try:return child.wait(timeout=timeout)
        except BaseException:
            os.killpg(child.pid,signal.SIGTERM)
            try:child.wait(timeout=30)
            except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
            raise


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--run',type=Path,required=True);ap.add_argument('--worker',required=True)
    a=ap.parse_args();out=a.run.resolve();queue=json.loads((out/'queue.json').read_text())
    allocation=next(w for w in queue['workers'] if w['name']==a.worker)
    host=socket.gethostname().split('.')[0]
    if host!=allocation['host'] or set(os.sched_getaffinity(0))!={allocation['controller_cpu']}:
        raise RuntimeError('Worker host/CPU identity mismatch')
    for path,digest in queue['source_sha256'].items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest()!=digest:
            raise RuntimeError('Prepared input/code changed before worker startup')
    worker_dir=out/'workers'/a.worker;worker_dir.mkdir(parents=True,exist_ok=True)
    ownership=(worker_dir/'worker.lock').open('a+')
    fcntl.flock(ownership,fcntl.LOCK_EX|fcntl.LOCK_NB)
    status=dict(status='running',worker=a.worker,host=host,pid=os.getpid(),process_start=process_start(),
                allocation=allocation,started=time.time(),items=[])
    sp=worker_dir/'status.json';save(sp,status)
    env=dict(os.environ,PYTHONPATH=str(ROOT/'src'))
    env['PATH']=str(Path(sys.executable).parent)+os.pathsep+env.get('PATH','')
    env['LD_LIBRARY_PATH']=(out/'library_path.txt').read_text().strip()+':'+env.get('LD_LIBRARY_PATH','')
    for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        env[key]='1'
    def interrupted(signum,frame):raise KeyboardInterrupt(f'signal {signum}')
    signal.signal(signal.SIGTERM,interrupted)
    consecutive=0;claim=None
    try:
        while time.time()<queue['stop_new_work_at_unix']:
            selected=None;cooldown=False;peers=False
            for gid in queue['ids']:
                done=out/'done'/f'{gid}.json'
                if done.exists():continue
                failure=out/'failures'/f'{gid}.json'
                prior=json.loads(failure.read_text()) if failure.exists() else {}
                if prior.get('attempts',0)>=2:continue
                if prior.get('retry_after',0)>time.time():cooldown=True;continue
                proposed=out/'claims'/gid
                try:proposed.mkdir()
                except FileExistsError:peers=True;continue
                claim=proposed;selected=gid
                save(claim/'owner.json',dict(host=host,pid=os.getpid(),process_start=process_start(),worker=a.worker))
                if done.exists():
                    import shutil;shutil.rmtree(claim);claim=None;selected=None;continue
                break
            if selected is None:
                if cooldown or peers:time.sleep(15);continue
                break
            gid=selected
            rows=[r for r in read_manifest() if r['id']==gid]
            if len(rows)!=queue['expected_states'][gid]:raise RuntimeError('Incomplete prepared manifest')
            item=out/'items'/gid;item.mkdir(exist_ok=True)
            attempt=len(list(item.glob('attempt_*')))+1;ad=item/f'attempt_{attempt:02d}';ad.mkdir()
            entry=dict(id=gid,status='preparing',attempt=attempt,started=time.time(),log_dir=str(ad))
            status['items'].append(entry);save(sp,status)
            try:
                archive=[sys.executable,str(ROOT/'scripts/pipeline/stage_uma_recalculation.py'),
                         '--only',gid,'--archive-root',str(ad/'superseded_uma'),'--execute']
                if run(archive,ad/'archive.log',env,180):raise RuntimeError('UMA archive preparation failed')
                pipeline=item/'pipeline'
                argv=['taskset','-c',','.join(map(str,allocation['compute_cpus'])),sys.executable,'-u',
                      str(ROOT/'scripts/pipeline/run_motifs.py'),'--config',str(ROOT/'config/custom_candidates.py'),
                      '--only',gid,'--out',str(pipeline),'--execute','--gpu',str(allocation['gpu']),
                      '--threads','2','--cpus',','.join(map(str,allocation['compute_cpus'])),
                      '--nconf','3','--thermal-definition','native','--skip-review']
                if pipeline.exists():argv.append('--resume')
                entry.update(status='computing',argv=argv);save(sp,status)
                contention=0
                while True:
                    rc=run(argv,ad/f'pipeline_{contention:02d}.log',env,12*3600)
                    manifest=json.loads((pipeline/'run_manifest.json').read_text()) if (pipeline/'run_manifest.json').exists() else {}
                    if rc==0:break
                    steps=manifest.get('steps',[])
                    last=Path(steps[-1]['log']) if steps else None
                    if last is None or not last.exists() or 'BUSY/reserved -> not running' not in last.read_text():
                        raise RuntimeError(f'Candidate stage failed, rc={rc}; outputs preserved')
                    if time.time()>=queue['stop_new_work_at_unix']:raise RuntimeError('Reservation contention deadline')
                    entry['status']='waiting_for_reservation';save(sp,status)
                    if '--resume' not in argv:argv.append('--resume')
                    contention+=1;time.sleep(15)
                from redox.validation.readiness import inspect_states
                qc=inspect_states(rows)
                if (manifest.get('status')!='selected_outputs_verified' or len(qc)!=len(rows)
                        or any(r['status']!='ready' for r in qc)):
                    raise RuntimeError('Zero exit did not produce the required verified state records')
                entry.update(status='selected_records_verified',finished=time.time(),readiness=qc,
                             pipeline_manifest=str(pipeline/'run_manifest.json'),scientific_accuracy_validated=False)
                save(out/'done'/f'{gid}.json',entry);consecutive=0
            except Exception as exc:
                failure=out/'failures'/f'{gid}.json'
                prior=json.loads(failure.read_text()) if failure.exists() else {}
                attempts=prior.get('attempts',0)+1
                entry.update(status='failed_preserved',error=repr(exc),finished=time.time())
                save(failure,dict(id=gid,attempts=attempts,retry_after=time.time()+300,last=entry))
                consecutive+=1
            finally:
                import shutil;shutil.rmtree(claim);claim=None;save(sp,status)
            if consecutive>=3:raise RuntimeError('Three consecutive scientific failures; worker circuit breaker')
        status.update(status='worker_finished',finished=time.time());save(sp,status)
    except (Exception,KeyboardInterrupt) as exc:
        status.update(status='stopped_preserving_outputs',error=repr(exc),finished=time.time());save(sp,status);raise
    finally:
        if claim is not None:
            import shutil;shutil.rmtree(claim)
        ownership.close()


if __name__=='__main__':main()
