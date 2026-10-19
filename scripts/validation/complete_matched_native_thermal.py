#!/usr/bin/env python
"""Bounded same-definition AQ/Fc native thermochemistry closure after metal control."""
import hashlib,json,os,signal,socket,subprocess,sys,time,site
from pathlib import Path
from redox.core.common import ROOT,read_result
from redox.qm.native_thermal import NATIVE_THERMAL,free_energy_native
BASE=ROOT/'results/validation/physical_thermal_spin_2026-10-10'
GATE='/nfs/lambda_stor_01/homes/rzhu/bin/gpu_reserve'


def main():
    p=BASE/'matched_native_thermal_status.json'
    status=dict(status='waiting_for_Fc_implementation_control',host=socket.gethostname(),pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().split()[21],started=time.time(),max_wait_seconds=21600,
        states=['ferrocene:neu','ferrocene:ox','aq_benzylamino:neu','aq_benzylamino:red1','aq_benzylamino:red2'],
        gpu=0,compute_cpu_affinity=[0,1],controller_cpu_affinity=sorted(os.sched_getaffinity(0)),
        nuclear_protocol=NATIVE_THERMAL,legacy_production_modified=False)
    def save():
        t=p.with_suffix('.tmp');t.write_text(json.dumps(status,indent=2));t.replace(p)
    save();child=None
    env=dict(os.environ,PYTHONPATH=str(ROOT/'src'))
    env['PATH']=str(Path(sys.executable).parent)+os.pathsep+env.get('PATH','')
    libs=[str(p) for base in site.getsitepackages() for p in Path(base).glob('nvidia/*/lib')]
    env['LD_LIBRARY_PATH']=':'.join(libs+[env.get('LD_LIBRARY_PATH','')])
    def stop(sig,frame):raise KeyboardInterrupt(f'termination {sig}')
    signal.signal(signal.SIGTERM,stop)
    try:
        deadline=time.time()+21600
        while time.time()<deadline:
            fp=BASE/'uks_fc_full_response_control/status.json'
            if fp.exists():
                fc=json.loads(fp.read_text())
                if fc['status']=='failed' or fc['status']=='uks_full_response_failed':raise RuntimeError('Fc metal-doublet Hessian control failed')
                proc=Path(f'/proc/{fc["pid"]}/stat')
                alive=proc.exists() and proc.read_text().split()[21]==fc['process_start']
                if fc['status']=='uks_full_response_verified' and not alive:break
            time.sleep(15)
        else:raise RuntimeError('Fc control deadline exceeded')
        with (BASE/'native_control_adoption.log').open('w') as log:
            result=subprocess.run([sys.executable,str(ROOT/'scripts/validation/adopt_verified_native_controls.py')],cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT)
        if result.returncode:raise RuntimeError('Verified control adoption failed')
        while time.time()<deadline:
            slots=subprocess.check_output([GATE,'list','--hosts','lambda1'],text=True)
            if 'lambda1 0' in slots.splitlines():break
            time.sleep(15)
        else:raise RuntimeError('Native cohort reservation deadline exceeded')
        if os.getloadavg()[0]+2>os.cpu_count():raise RuntimeError('Native cohort CPU capacity changed')
        # The native writer independently refuses incomplete/stale geometry and bad derivatives.
        argv=[GATE,'run','0','--label','redox-matched-aq-fc-native-thermal-20261011','--','taskset','-c','0,1',
            'timeout','--signal=TERM','--kill-after=30','2h',sys.executable,'-u','-m','redox.qm.native_thermal',
            '--only','ferrocene,aq_benzylamino','--threads','2']
        status.update(status='computing_matched_native_thermal',argv=argv);save()
        for name in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[name]='2'
        with (BASE/'matched_native_thermal.log').open('w') as log:
            child=subprocess.Popen(argv,cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            status['gate_pid']=child.pid;save();rc=child.wait()
        if rc:raise RuntimeError('Matched native thermochemistry failed; inspect exact state/derivative error')
        G={}
        for key in status['states']:
            gid,state=key.split(':');G[key],qc=free_energy_native(gid,state)
            if G[key] is None:raise RuntimeError(f'{key}: native record validation failed {qc}')
        fc=G['ferrocene:ox']-G['ferrocene:neu']
        first=G['aq_benzylamino:neu']-G['aq_benzylamino:red1']-fc
        second=G['aq_benzylamino:red1']-G['aq_benzylamino:red2']-fc
        spacing=G['aq_benzylamino:neu']+G['aq_benzylamino:red2']-2*G['aq_benzylamino:red1']
        status.update(status='matched_native_predictions_verified',free_energies_eV=G,Fc_absolute_V=fc,
            E1_V_vs_Fc=first,E2_V_vs_Fc=second,wave_spacing_eV=spacing,
            paired_arithmetic_error_eV=spacing-(first-second),finished=time.time(),
            model_experimental_accuracy_validated=False)
        if abs(status['paired_arithmetic_error_eV'])>1e-8:raise RuntimeError('Matched native spacing identity failed')
        save()
    except (Exception,KeyboardInterrupt) as exc:
        status.update(status='failed',error=repr(exc),finished=time.time());save();raise
    finally:
        if child is not None and child.poll() is None:
            os.killpg(child.pid,signal.SIGTERM)
            try:child.wait(timeout=30)
            except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()


if __name__=='__main__':main()
