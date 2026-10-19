#!/usr/bin/env python
"""Bounded UMA-only pilot for the two registered motifs after physical control releases GPU1."""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import socket
import signal
import subprocess
import sys
import time

from redox.core.common import ROOT, UMA, read_manifest

IDS=['exp_viologen_methoxyethyl','exp_phthalimide_control']
OUT=ROOT/'results/discovery/uma_pilot_2026-10-11'
GATE='/nfs/lambda_stor_01/homes/rzhu/bin/gpu_reserve'


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    status=dict(status='waiting_for_physical_control',host=socket.gethostname(),pid=os.getpid(),
                process_start=Path('/proc/self/stat').read_text().split()[21],started=time.time(),ids=IDS,
                gpu=1,threads=2,worker_cpus=[2,3],max_wait_seconds=21600,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),steps=[])
    def save():
        p=OUT/'status.json';t=p.with_suffix('.tmp');t.write_text(json.dumps(status,indent=2));t.replace(p)
    handles=[]
    child=None
    def terminate(signum,frame):raise KeyboardInterrupt(f'termination signal {signum}')
    signal.signal(signal.SIGTERM,terminate)
    try:
        save();deadline=time.time()+21600
        required=[ROOT/'results/validation/physical_reorg_2026-10-10/conformer_controls'/mode/'result.json'
                  for mode in ('local_unconstrained','linker_fixed')]
        while time.time()<deadline:
            if all(p.exists() and json.loads(p.read_text()).get('converged_smd') is True for p in required):
                break
            if any(p.with_name('failure.json').exists() for p in required):
                raise RuntimeError('Physical control failed; inspect before starting dependent UMA pilot')
            time.sleep(15)
        else:raise RuntimeError('physical control deadline expired')
        for gid in IDS:
            h=(ROOT/'logs/pipeline_locks'/f'{gid}.lock').open('a+')
            fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB);handles.append(h)
        while time.time()<deadline:
            free=subprocess.check_output([GATE,'list','--hosts','lambda1'],text=True)
            if 'lambda1 1' in free.splitlines():break
            time.sleep(15)
        else:raise RuntimeError('reservation deadline expired')
        if os.getloadavg()[0]+2>(os.cpu_count() or 1):raise RuntimeError('insufficient current CPU capacity')
        status['status']='running';save()
        # Claim GPU once for the complete small pilot so another worker cannot interleave
        # conflicting state writes. CUDA_VISIBLE_DEVICES is set only by the gate.
        command=[GATE,'run','1','--label','redox-new-motif-uma-20261011','--','taskset','-c','2,3',
                 'timeout','--signal=TERM','--kill-after=60','3h','bash','-c',
                 'set -e; "$1" -u -m redox.qm.uma --only exp_viologen_methoxyethyl --device cuda; '
                 '"$1" -u -m redox.qm.uma --only exp_phthalimide_control --device cuda',
                 'uma-pilot',sys.executable]
        status['command']=command;save()
        with (OUT/'worker.log').open('w') as stream:
            child=subprocess.Popen(command,cwd=ROOT,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True)
            status['gate_pid']=child.pid;save();rc=child.wait()
        status['returncode']=rc
        if rc:raise RuntimeError(f'UMA pilot exited {rc}; inspect worker.log')
        states=[]
        for row in read_manifest():
            if row['id'] not in IDS:continue
            sd=UMA/row['id']/row['state'];p=sd/'result.json'
            if not p.exists() or not (sd/'relaxed.xyz').exists():raise RuntimeError(f'missing UMA output {sd}')
            rec=json.loads(p.read_text())
            if rec.get('converged') is not True or rec.get('charge')!=int(row['charge']):
                raise RuntimeError(f'bad UMA convergence/charge {sd}')
            states.append(dict(id=row['id'],state=row['state'],chosen_mult=rec['chosen_mult'],
                               energy_eV=rec['energy_eV'],spin_gap_eV=rec.get('spin_gap_eV'),
                               semantics='gas-phase preoptimization only; no validated redox or kinetics'))
        if len(states)!=5:raise RuntimeError('unexpected UMA state count')
        status.update(status='warm_start_outputs_verified',states=states,finished=time.time(),
                      production_properties_computed=False);save()
    except (Exception,KeyboardInterrupt) as exc:
        status.update(status='failed',error=repr(exc),finished=time.time());save();raise
    finally:
        if child is not None and child.poll() is None:
            os.killpg(child.pid,signal.SIGTERM)
            try:child.wait(timeout=30)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid,signal.SIGKILL);child.wait()
        for h in handles:h.close()


if __name__=='__main__':main()
