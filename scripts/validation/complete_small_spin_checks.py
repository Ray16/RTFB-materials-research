#!/usr/bin/env python
"""Level-matched DFT checks of all remaining small UMA spin gaps, one reserved GPU.

Alternative geometries/records remain separate. Production states are not replaced.
Only electronic SMD gaps are certified; thermochemical gaps require their own validation.
"""
import fcntl,hashlib,json,os,socket,subprocess,sys,time
from pathlib import Path
from redox.core.common import ROOT,DFT,read_result,smiles_to_xyz_map,read_manifest
from redox.core.protocol import load_record,geom_sha1
from redox.qm.dft import dft_smd,_opt_protocol_mismatch
from redox.qm.sp import ensure

BASE=ROOT/'results/validation/physical_thermal_spin_2026-10-10'
OUT=BASE/'all_small_spin_checks'
GATE='/nfs/lambda_stor_01/homes/rzhu/bin/gpu_reserve'
IDS=['d3_80izcq','aq_benzyloxy','pmdi','mophquinone','dmophquinone','d3_05stji']


def worker():
    status=dict(status='running',host=socket.gethostname(),pid=os.getpid(),
                process_start=Path('/proc/self/stat').read_text().split()[21],
                gpu=os.getenv('CUDA_VISIBLE_DEVICES'),cpus=sorted(os.sched_getaffinity(0)),
                started=time.time(),expected_states=IDS,items=[])
    def save():
        p=OUT/'worker_status.json';tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(status,indent=2));tmp.replace(p)
    OUT.mkdir(parents=True,exist_ok=True);save();failures=0
    for gid in IDS:
        state='red2';chosen=read_result(gid,state);q=int(chosen['charge']);m=3 if chosen['mult']==1 else 1
        sd=ROOT/'calcs/spincheck'/gid/f'{state}_m{m}'
        item=dict(id=gid,state=state,chosen_mult=chosen['mult'],alternative_mult=m,
                  source_geometry_sha1=geom_sha1(DFT/gid/state/'opt.xyz'),started=time.time())
        status['items'].append(item);save()
        lock=(ROOT/'logs/pipeline_locks'/f'{gid}_spincheck.lock').open('a+')
        try:
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
            sd.mkdir(parents=True,exist_ok=True);rj=sd/'result.json'
            if rj.exists():
                raw=json.loads(rj.read_text())
                if _opt_protocol_mismatch(raw,q) or raw['charge']!=q or raw['mult']!=m:
                    raise RuntimeError('existing alternative has a different geometry protocol; preserve and inspect')
            else:
                if (sd/'opt.xyz').exists():raise RuntimeError('partial alternative geometry needs inspection')
                result=dft_smd(DFT/gid/state/'opt.xyz',q,m,backend='gpu',opt_out=sd/'opt.xyz',do_freq=False)
                if result.get('converged_smd') is not True:raise RuntimeError('alternative SCF did not converge')
                result.update(id=gid,state=f'{state}_m{m}',parent_state=state,
                              source_geometry_sha1=item['source_geometry_sha1'],g_thermal_eV=None,
                              validation='electronic spin comparison; thermal correction pending')
                rj.write_text(json.dumps(result,indent=2)+'\n')
            smi=next(r['smiles'] for r in read_manifest() if r['id']==gid)
            if smiles_to_xyz_map(smi,sd/'opt.xyz') is None:raise RuntimeError('alternative changes molecular connectivity')
            _,record=ensure(gid,f'{state}_m{m}',sd,q,m,backend='gpu')
            if load_record(sd,q,m) is None:raise RuntimeError('alternative active SP invalid')
            gap=record['e_smd_eV']-chosen['e_smd_eV']
            item.update(status='electronic_gap_verified',gap_eV=gap,lower_alternative=(gap<-1e-4),
                        thermochemical_gap_validated=False,finished=time.time(),production_replaced=False)
            failures=0;print('SPIN_CHECK_COMPLETE',gid,m,gap,flush=True)
        except Exception as exc:
            item.update(status='failed',error=repr(exc),finished=time.time());failures+=1
            print('SPIN_CHECK_FAILED',gid,repr(exc),flush=True)
        finally:lock.close();save()
        if failures>=3:
            status['circuit_breaker']='three consecutive failures';break
    status.update(status='electronic_checks_completed' if len(status['items'])==len(IDS)
                  and all(r['status']=='electronic_gap_verified' for r in status['items']) else 'incomplete',
                  finished=time.time());save()
    if status['status']!='electronic_checks_completed':raise SystemExit(1)


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    if '--worker' in sys.argv:worker();return
    status=dict(status='waiting_for_numerical_columns',host=socket.gethostname(),pid=os.getpid(),
                process_start=Path('/proc/self/stat').read_text().split()[21],started=time.time(),
                gpu=1,cpus=[2,3],max_wait_seconds=21600,states=IDS)
    path=OUT/'coordinator_status.json'
    def save():
        tmp=path.with_suffix('.tmp');tmp.write_text(json.dumps(status,indent=2));tmp.replace(path)
    save();deadline=time.time()+21600
    try:
        shard=BASE/'numerical_solution_hessian_aq_red2_tight/shard_1_status.json'
        while time.time()<deadline:
            if shard.exists():
                raw=json.loads(shard.read_text())
                if raw['status']=='failed':raise RuntimeError('numerical dependency failed')
                if raw['status']=='columns_verified':break
            time.sleep(15)
        else:raise RuntimeError('dependency deadline expired')
        while time.time()<deadline:
            slots=subprocess.check_output([GATE,'list','--hosts','lambda1'],text=True)
            if 'lambda1 1' in slots.splitlines():break
            time.sleep(15)
        else:raise RuntimeError('reservation deadline expired')
        if os.getloadavg()[0]+2>(os.cpu_count() or 1):raise RuntimeError('insufficient current CPU capacity')
        argv=[GATE,'run','1','--label','redox-small-spin-checks-20261011','--','taskset','-c','2,3',
              'timeout','--signal=TERM','--kill-after=60','6h',sys.executable,'-u',str(Path(__file__).resolve()),'--worker']
        status.update(status='running',argv=argv);save()
        with (OUT/'worker.log').open('w') as log:
            p=subprocess.Popen(argv,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
            status['gate_pid']=p.pid;save();rc=p.wait()
        status.update(status='completed' if rc==0 else 'failed',returncode=rc,finished=time.time());save()
    except Exception as exc:
        status.update(status='failed',error=repr(exc),finished=time.time());save();raise


if __name__=='__main__':main()
