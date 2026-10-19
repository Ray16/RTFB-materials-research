#!/usr/bin/env python
"""Start the second conservative-gradient Hessian shard after the bounded spin queue."""
import json,os,signal,socket,subprocess,sys,time
from pathlib import Path
from redox.core.common import ROOT
BASE=ROOT/'results/validation/physical_thermal_spin_2026-10-10'
GATE='/nfs/lambda_stor_01/homes/rzhu/bin/gpu_reserve'


def main():
    report=dict(status='waiting_for_small_spin_checks',host=socket.gethostname(),pid=os.getpid(),
                started=time.time(),max_wait_seconds=21600,gpu=1,cpus=[2,3])
    path=BASE/'full_response_second_shard_queue.json'
    def save():
        tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(report,indent=2));tmp.replace(path)
    child=None
    def stop(sig,frame):raise KeyboardInterrupt(f'termination {sig}')
    signal.signal(signal.SIGTERM,stop)
    try:
        save();deadline=time.time()+21600
        while time.time()<deadline:
            p=BASE/'all_small_spin_checks/worker_status.json'
            if p.exists():
                current=json.loads(p.read_text())
                if current['status'] in ('electronic_checks_completed','incomplete'):break
            time.sleep(15)
        else:raise RuntimeError('spin dependency deadline expired')
        while time.time()<deadline:
            free=subprocess.check_output([GATE,'list','--hosts','lambda1'],text=True)
            if 'lambda1 1' in free.splitlines():break
            time.sleep(15)
        else:raise RuntimeError('reservation deadline expired')
        if os.getloadavg()[0]+2>(os.cpu_count() or 1):raise RuntimeError('CPU capacity changed')
        argv=[GATE,'run','1','--label','redox-full-response-hessian1-20261011','--','taskset','-c','2,3',
              'timeout','--signal=TERM','--kill-after=30','3h',sys.executable,'-u',
              str(ROOT/'scripts/validation/numerical_hessian_shard.py'),'--shard','1','--grid-response']
        report.update(status='running',argv=argv);save()
        with (BASE/'full_response_second_shard.log').open('w') as stream:
            child=subprocess.Popen(argv,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
            report['gate_pid']=child.pid;save();rc=child.wait()
        report.update(status='completed' if rc==0 else 'failed',returncode=rc,finished=time.time());save()
    except (Exception,KeyboardInterrupt) as exc:
        report.update(status='failed',error=repr(exc),finished=time.time());save();raise
    finally:
        if child is not None and child.poll() is None:
            os.killpg(child.pid,signal.SIGTERM)
            try:child.wait(timeout=30)
            except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()


if __name__=='__main__':main()
