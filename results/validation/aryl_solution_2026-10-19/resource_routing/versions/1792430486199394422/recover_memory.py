"""Bounded H100 recovery; original scientific producers and failed output are retained."""
import argparse
import json
import os
from pathlib import Path
import runpy
import shlex
import signal
import socket
import subprocess
import sys
import time

from redox.core.common import ROOT, DFT
from redox.core.io import atomic_write_json
from redox.qm.uma import file_sha256
from redox.workflows.checkpoints import verify_state_checkpoint
from redox.workflows.resources import (allocation_guard, cpu_allocation, free_slots,
    gpu_suitability, host_capacity, verify_gate_owner)

BASE = Path(__file__).resolve().parent
RECOVERY = BASE / 'resource_routing'


def identity():
    return dict(host=socket.gethostname().split('.')[0], pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19],
        cpus=sorted(os.sched_getaffinity(0)))


def old_owner_absent(owner, key):
    """Fail closed on unreachable hosts, PID reuse, or an extant state predictor."""
    code = '''import json,pathlib,sys,psutil
o=json.loads(sys.argv[1]); key=sys.argv[2]; p=pathlib.Path('/proc')/str(o['pid'])/'stat'
if p.exists():
 f=p.read_text().rsplit(')',1)[1].split()
 if f[19]==o['process_start'] and f[0]!='Z':raise RuntimeError('Original owner still live')
selection=key.replace('__',':')
for proc in psutil.process_iter(['cmdline']):
 argv=proc.info['cmdline'] or []
 if selection in argv:raise RuntimeError('A state predictor is still live')
print('original_owner_and_state_predictor_absent')'''
    argv=[sys.executable, '-c', code, json.dumps(owner), key]
    if owner['host'] != socket.gethostname().split('.')[0]:
        argv=['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=8', owner['host'], shlex.join(argv)]
    result=subprocess.run(argv, check=True, text=True, capture_output=True, timeout=20)
    if result.stdout.strip() != 'original_owner_and_state_predictor_absent':
        raise RuntimeError('Cannot establish prior ownership release')


def archive_attempt(task, old, destination):
    """Only incomplete, owner-released records may leave their canonical paths."""
    key=task['key']; sd=DFT/task['id']/task['state']
    if (BASE/'done'/f'{key}.json').exists():
        raise RuntimeError('Completed state must be verified and retained')
    if old.get('status') != 'failed_preserved':
        raise RuntimeError('An existing attempt is not a preserved terminal failure')
    old_owner_absent(old, key)
    # All four resource-stopped attempts failed before canonical state publication.
    # A result requires a separate, scientifically reviewed resume rather than replay.
    if (sd/'result.json').exists():
        raise RuntimeError('Published state result requires separate resume review')
    paths=[BASE/'tasks'/key, BASE/'failures'/f'{key}.json', sd]
    source_hashes={str(p):file_sha256(p) for path in paths
                   for p in (path.rglob('*') if path.is_dir() else [path]) if p.is_file()}
    destination.mkdir(parents=True, exist_ok=False)
    atomic_write_json(destination/'archive.json', dict(key=key, source_sha256=source_hashes,
        reason='Memory-class migration; replay identical three seeds on H100',
        ownership_verified_absent=True, physical_thresholds_changed=False, archived=time.time()))
    for path,label in zip(paths, ('task', 'failure.json', 'calculation')):
        if path.exists():path.rename(destination/label)
    for source,h in source_hashes.items():
        original=Path(source)
        for path,label in zip(paths, ('task','failure.json','calculation')):
            if original == path or path in original.parents:
                target=destination/label/ original.relative_to(path) if original!=path else destination/label
                if file_sha256(target)!=h:raise RuntimeError('Archive bytes changed')
                break


def worker(index):
    verify_gate_owner()
    device=gpu_suitability(int(os.environ['CUDA_VISIBLE_DEVICES']))
    if 'H100' not in device['name'] or device['memory_MiB']<90*1024:
        raise RuntimeError('Memory recovery requires the observed H100 NVL memory class')
    plan=json.loads((BASE/'plan.json').read_text()); task=plan['tasks'][index]; key=task['key']
    assert all(file_sha256(p)==h for p,h in plan['producer_sha256'].items())
    from redox.qm.native_thermal import implementation_controls
    implementation_controls()
    claim=BASE/'claims'/key
    try:claim.mkdir()
    except FileExistsError:return
    owner=identity(); owner['gpu']=int(os.environ['CUDA_VISIBLE_DEVICES'])
    atomic_write_json(claim/'owner.json',owner)
    status_path=RECOVERY/f'recovery_{key}.json'
    report=dict(owner,key=key,status='preparing',started=time.time(),attempt_limit=1)
    atomic_write_json(status_path,report)
    def stop(signum,frame):raise KeyboardInterrupt(f'signal {signum}')
    signal.signal(signal.SIGTERM,stop)
    try:
        complete=BASE/'done'/f'{key}.json'
        if complete.exists():
            verify_state_checkpoint(complete); report['status']='existing_completion_verified'; return
        state=BASE/'tasks'/key/'status.json'
        if state.exists():
            old=json.loads(state.read_text())
            archive_attempt(task,old,RECOVERY/'archives'/key/str(time.time_ns()))
        elif (DFT/task['id']/task['state']).exists():
            raise RuntimeError('Calculation without task evidence requires manual ownership review')
        report.update(status='computing',device=device)
        atomic_write_json(status_path,report)
        sys.argv=[str(BASE/'worker.py'),'--task',str(index)]
        runpy.run_path(str(BASE/'worker.py'),run_name='__main__')
        result=json.loads(state.read_text());assert result['status']=='state_records_verified'
        atomic_write_json(complete,dict(key=key,owner=owner,status_path=str(state),
            status_sha256=file_sha256(state),finished=time.time(),recovery=str(status_path)))
        verify_state_checkpoint(complete)
        report['status']='state_records_verified'
    except BaseException as exc:
        report.update(status='failed_preserved',error=repr(exc))
        atomic_write_json(BASE/'failures'/f'{key}.json',dict(key=key,owner=owner,error=repr(exc),
            state_status=str(BASE/'tasks'/key/'status.json'),attempts=1,recovery=str(status_path),finished=time.time()))
        raise
    finally:
        report['finished']=time.time();atomic_write_json(status_path,report)
        # worker.py joins its children before returning or raising.
        (claim/'owner.json').unlink();claim.rmdir()


def dispatch(hours):
    import fcntl
    lock=(RECOVERY/'dispatcher.lock').open('a+')
    fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    plan=json.loads((RECOVERY/'recovery_plan.json').read_text())
    assert file_sha256(__file__)==plan['recovery_source_sha256']
    assert file_sha256(BASE/'plan.json')==plan['original_plan_sha256']
    preflight=json.loads((RECOVERY/'lambda13_preflight.json').read_text())
    assert preflight['ok'] and preflight['host']==socket.gethostname().split('.')[0]
    assert preflight['python']==sys.executable
    report=dict(identity(),status='checking_resources',workers=[],started=time.time(),
        deadline=time.time()+hours*3600,plan=str(RECOVERY/'recovery_plan.json'))
    atomic_write_json(RECOVERY/'dispatcher.json',report)
    children=[]
    def stop(signum,frame):raise KeyboardInterrupt(f'signal {signum}')
    signal.signal(signal.SIGTERM,stop)
    try:
        while time.time()<report['deadline']:
            inflight=set();setup_failures=0
            for submitted in report['workers']:
                stat=Path(f"/proc/{submitted['controller_pid']}/stat")
                if stat.exists():
                    fields=stat.read_text().rsplit(')',1)[1].split()
                    if fields[19]==submitted['process_start'] and fields[0]!='Z':
                        inflight.add(submitted['key']);continue
                evidence=RECOVERY/f"recovery_{submitted['key']}.json"
                if not evidence.exists():
                    execution=Path(submitted['reservation'])/'execution.json'
                    record=json.loads(execution.read_text()) if execution.exists() else {}
                    # Contention has no molecular attempt and remains eligible.
                    contention='occupied or reserved' in record.get('error','')
                    log=Path(submitted['reservation'])/'command.log'
                    if log.exists() and 'BUSY/reserved -> not running' in log.read_text():contention=True
                    submitted['preexecution_outcome']='contention_no_attempt' if contention else 'setup_failure_no_attempt'
                    if not contention:setup_failures+=1
            if setup_failures>=3:
                raise RuntimeError('Three pre-execution failures; stop dispatcher without consuming molecular inputs')
            pending=[]
            for target in plan['tasks']:
                key=target['key']
                done=BASE/'done'/f'{key}.json'
                if done.exists():verify_state_checkpoint(done);continue
                if (RECOVERY/f'recovery_{key}.json').exists():continue
                if key in inflight:continue
                if (BASE/'claims'/key).exists():continue
                pending.append(target)
            if not pending:
                report['status']='all_recovery_targets_started_or_owned_by_original_peers';break
            with allocation_guard():
                slots=free_slots([socket.gethostname().split('.')[0]])
                suitable=[]
                for host,gpu in slots:
                    try:device=gpu_suitability(gpu)
                    except RuntimeError:continue
                    if 'H100' in device['name'] and device['memory_MiB']>=90*1024:suitable.append((gpu,device))
                count=min(len(pending),len(suitable))
                if count:
                    resources=host_capacity(2*count,minimum_ram_gib=32*count)
                    cpus=cpu_allocation(2*count)
                    for i,(target,(gpu,device)) in enumerate(zip(pending,suitable)):
                        assigned=cpus[2*i:2*i+2]
                        reserved=RECOVERY/f"reserved_{target['key']}_{time.time_ns()}"
                        argv=[sys.executable,'-u','-m','redox.workflows.gpu','--gpu',str(gpu),
                            '--cpus',','.join(map(str,assigned)),'--threads','2','--out',str(reserved),
                            '--',sys.executable,'-u',str(Path(__file__).resolve()),'worker','--task',str(target['index'])]
                        env=dict(os.environ,PYTHONPATH=str(ROOT/'src'),REDOX_HESSIAN_MEMORY_GIB='24')
                        for name in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[name]='2'
                        with reserved.with_suffix('.log').open('w') as stream:
                            proc=subprocess.Popen(argv,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,
                                stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
                        children.append(proc)
                        report['workers'].append(dict(key=target['key'],gpu=gpu,device=device,
                            cpu_affinity=assigned,controller_pid=proc.pid,argv=argv,reservation=str(reserved),
                            process_start=Path(f'/proc/{proc.pid}/stat').read_text().rsplit(')',1)[1].split()[19],
                            resources=resources,submitted=time.time()))
                        atomic_write_json(RECOVERY/'dispatcher.json',report)
            report.update(status='waiting_for_free_H100' if not count else 'workers_submitted_waiting_for_remaining_H100',
                pending_keys=[t['key'] for t in pending],checked=time.time())
            atomic_write_json(RECOVERY/'dispatcher.json',report)
            time.sleep(30)
        else:report['status']='deadline_reached_remaining_targets_preserved'
    except BaseException as exc:
        report.update(status='dispatcher_stopped_workers_retained',error=repr(exc));raise
    finally:
        report['finished']=time.time();atomic_write_json(RECOVERY/'dispatcher.json',report)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['dispatch','worker'])
    parser.add_argument('--task',type=int);parser.add_argument('--hours',type=float,default=12)
    args=parser.parse_args()
    if args.mode=='worker':worker(args.task)
    else:dispatch(args.hours)
