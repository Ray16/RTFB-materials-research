#!/usr/bin/env python
"""Finite electronic-structure pilot for polymer folding and segment association.

No Hessian or association free energy is claimed at this stage. Larger models
depend on an H100 energy/gradient memory pilot before parallel execution.
"""
import argparse
import gc
import hashlib
import json
import math
import os
from pathlib import Path
import shlex
import signal
import socket
import subprocess
import sys
import threading
import time
import traceback
from contextlib import redirect_stdout,redirect_stderr

from redox.core.common import ROOT
from redox.core.io import atomic_write_json

BASE=ROOT/'results/polymer_models/campaigns/2026-10-19_pi_stacking'
SCRIPT=Path(__file__).resolve()
PYTHON=Path('/nfs/lambda_stor_01/homes/rzhu/miniforge3/envs/redox/bin/python')
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify(mapping):
    for p,h in mapping.items():
        if sha(p)!=h:raise RuntimeError('Frozen source/input changed: '+p)


def freeze():
    assert not (BASE/'plan.json').exists(),'Preserve the existing frozen campaign'
    from rdkit import Chem
    cases=json.loads((BASE/'cases.json').read_text())['cases']
    tasks=[]
    for case in cases:
        count=Chem.AddHs(Chem.MolFromSmiles(case['smiles'])).GetNumAtoms()
        for target in (('monomer',) if case['repeats']==1 else ('open','stacked')):
            tasks.append(dict(key=case['key']+'__'+target,candidate_id=case['candidate_id'],
                model_key=case['key'],model_type=case['model_type'],repeats=case['repeats'],
                stereoisomer=case['stereoisomer'],target=target,smiles=case['smiles'],
                fragment=case['fragment'],atoms=count,charge=0,mult=1,
                role='interaction_method_control; not synthesis admission'))
        if case['repeats']==1:
            for orientation in ('parallel','antiparallel'):
                model=case['candidate_id']+'__segments__'+orientation
                tasks.append(dict(key=model+'__'+orientation,candidate_id=case['candidate_id'],
                    model_key=model,model_type='two_nonbonded_segments',repeats=2,
                    stereoisomer='two achiral segments',target=orientation,
                    smiles=case['smiles']+'.'+case['smiles'],fragment=case['fragment'],
                    atoms=2*count,charge=0,mult=1,role='interaction_method_control; not synthesis admission'))
    assert len(tasks)==27
    files=list((ROOT/'src/redox').rglob('*.py'))+list((ROOT/'config').glob('*.py'))+list((ROOT/'config').glob('*.json'))
    files += [SCRIPT,ROOT/'scripts/pipeline/collect_pi_stacking.py',ROOT/'scripts/pipeline/prepare_pi_stacking.py']
    plan=dict(schema=1,created=time.time(),tasks=tasks,pilot_key=max(tasks,key=lambda t:t['atoms'])['key'],
              input_sha256={str(BASE/'cases.json'):sha(BASE/'cases.json')},
              producer_sha256={str(p):sha(p) for p in sorted(set(files))},
              raw_root=str(ROOT/'calcs/polymer_models'/BASE.name),solvent='acetonitrile',
              electronic_protocol='Production r2SCAN-D4/SVP geometry and wB97M-V/TZVPD/SMD energy',
              property_scope='Resting-state electronic folding/association screen; not deltaG or population',
              Hessians=False,maximum_attempts_per_item=1,consecutive_failure_limit=3,
              large_model_memory_pilot_required=True,maximum_memory_pilot_MiB=80*1024,
              small_RKS_geometry_SP_max_atoms=60,small_gpu_minimum_MiB=31*1024)
    for name in ('claims','tasks','done','failures','workers','launches','pilot','live_inputs'):(BASE/name).mkdir(exist_ok=True)
    atomic_write_json(BASE/'plan.json',plan)
    stage_inputs(plan)
    print('Frozen',len(tasks),'neutral electronic-structure cases; largest',max(t['atoms'] for t in tasks),'atoms')


def stage_inputs(plan):
    """Bind an immutable warm start only when its upstream calculation is ready."""
    from prepare_pi_stacking import choose
    full=BASE/'preparation.json'
    if full.exists():
        available=json.loads(full.read_text())['tasks']
        candidates=[(t,Path(t['input_xyz']),full) for t in available]
    else:
        candidates=[]
        for task in plan['tasks']:
            if task['target']!='monomer':continue
            source=BASE/'seeds'/task['model_key']/'ensemble.json'
            if source.exists():
                record,fallback=choose(json.loads(source.read_text()),'monomer')
                candidate=dict(task,pi_groups=record['pi_groups'],seed_stage='converged_topology_checked_xTB_ALPB')
                candidates.append((candidate,Path(record['xyz']),source))
    roster={t['key']:t for t in plan['tasks']}
    for task,xyz,upstream in candidates:
        expected=roster[task['key']]
        if any(task.get(k)!=v for k,v in expected.items()):raise RuntimeError('Prepared task differs from frozen roster: '+task['key'])
        directory=BASE/'live_inputs'/task['key'];receipt=directory/'task.json'
        if receipt.exists():
            old=json.loads(receipt.read_text());verify(old['input_sha256'])
            if sha(xyz)!=sha(old['input_xyz']):raise RuntimeError('Later preparation changed an already submitted warm start')
            continue
        directory.mkdir(exist_ok=True)
        import shutil
        destination=directory/'input.xyz';shutil.copyfile(xyz,destination)
        item=dict(task,input_xyz=str(destination),input_sha256={str(destination):sha(destination),str(upstream):sha(upstream),str(xyz):sha(xyz)})
        atomic_write_json(receipt,item)


def prepared_task(task):
    path=BASE/'live_inputs'/task['key']/'task.json'
    if not path.exists():return None
    ready=json.loads(path.read_text())
    if any(ready.get(k)!=v for k,v in task.items()):raise RuntimeError('Ready task does not match frozen roster')
    verify(ready['input_sha256'])
    return ready


class MemoryMonitor:
    def __init__(self):self.peak=0;self.stop=threading.Event();self.thread=None
    def __enter__(self):
        def poll():
            while not self.stop.is_set():
                try:
                    text=subprocess.check_output(['nvidia-smi','--query-compute-apps=pid,used_gpu_memory','--format=csv,noheader,nounits'],text=True,timeout=10)
                    for line in text.splitlines():
                        pid,memory=map(str.strip,line.split(','))
                        if int(pid)==os.getpid():self.peak=max(self.peak,int(memory))
                except (OSError,ValueError,subprocess.SubprocessError):pass
                self.stop.wait(1)
        self.thread=threading.Thread(target=poll,daemon=True);self.thread.start();return self
    def __exit__(self,*exc):self.stop.set();self.thread.join(timeout=12)


def gradient(xyz):
    from ase.io import read
    import numpy as np
    from pyscf import gto
    from redox.qm.dft import _build_smd_mf,_free_gpu_pool
    atoms=read(xyz)
    mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions.tolist())),unit='Angstrom',
              charge=0,spin=0,basis='def2-svp',verbose=0,max_memory=6000)
    mf=_build_smd_mf(mol,'r2scan','acetonitrile','gpu',disp='d4')
    mf.conv_tol=1e-11;mf.max_cycle=200;mf.kernel()
    if not mf.converged:raise RuntimeError('Independent geometry-gradient SCF failed')
    derivative=mf.nuc_grad_method();derivative.grid_response=True
    raw=derivative.kernel()
    value=np.asarray(raw.get() if hasattr(raw,'get') else raw)
    if not np.isfinite(value).all():raise RuntimeError('Nonfinite gradient')
    maximum=float(abs(value).max())
    del derivative,mf;_free_gpu_pool()
    return maximum


def calculate(task,plan,pilot=False):
    from redox.workflows.resources import verify_gate_owner
    import pyscf.lib
    from redox.qm.dft import dft_smd
    from redox.qm.sp import ensure
    from redox.core.protocol import load_record
    from redox.validation.protonation import hydrogen_topology_status
    from redox.build.polymer_contacts import contact_metrics
    from ase.io import read
    verify_gate_owner();pyscf.lib.num_threads(2)
    verify(task['input_sha256'])
    directory=BASE/'pilot' if pilot else BASE/'tasks'/task['key']
    if pilot:
        assert not (directory/'status.json').exists(),'Prior resource pilot must be reviewed'
    else:directory.mkdir(exist_ok=False)
    sd=Path(plan['raw_root'])/task['key']/'rest'
    if not pilot:sd.mkdir(parents=True,exist_ok=False)
    report=dict(task=task,status='computing',stage='memory_pilot' if pilot else 'dft_geometry',
                started=time.time(),host=socket.gethostname(),pid=os.getpid(),
                free_energy_computed=False,minimum_Hessian_tested=False,synthesis_admission=False)
    def save():atomic_write_json(directory/'status.json',report)
    save()
    memory=None
    try:
        with MemoryMonitor() as memory, (directory/'calculation.log').open('w') as log,redirect_stdout(log),redirect_stderr(log):
            if pilot:
                result=dft_smd(Path(task['input_xyz']),0,1,do_opt=False,do_freq=False,backend='gpu')
                if not (result['converged_smd'] and result['converged_gas']):raise RuntimeError('Memory-pilot energy SCF failed')
                gmax=gradient(task['input_xyz'])
                atomic_write_json(directory/'energy_and_gradient.json',dict(result=result,gradient_max_Ha_Bohr=gmax,
                    geometry='starting geometry only; resource qualification, not physical minimum'))
            else:
                result=dft_smd(Path(task['input_xyz']),0,1,do_opt=True,do_freq=False,backend='gpu',
                    opt_out=sd/'opt.xyz',conv_params=dict(convergence_gmax=1e-4,convergence_grms=6e-5,
                    convergence_energy=1e-7,convergence_dmax=.002,convergence_drms=.001))
                result.update(id=task['key'],state='rest',n_e=0)
                atomic_write_json(sd/'result.json',result)
                topology=hydrogen_topology_status(task['smiles'],sd/'opt.xyz')
                if topology!='ok':raise RuntimeError('Topology/hydrogen change: '+topology)
                ensure(task['key'],'rest',sd,0,1,backend='gpu')
                record=load_record(sd,0,1)
                if record is None:raise RuntimeError('Active electronic SP not admitted')
                gmax=gradient(sd/'opt.xyz')
                if gmax>.00045:raise RuntimeError('Independent solution gradient is not stationary')
                xyz=read(sd/'opt.xyz')
                metrics=contact_metrics(xyz.positions,task['pi_groups'])
                report.update(E_smd_eV=record['e_smd_eV'],gradient_max_Ha_Bohr=gmax,
                              contact_metrics=metrics,stacked_pairs=sum(m['stack_geometry'] for m in metrics),
                              state_dir=str(sd),topology_QC=topology)
        if memory.peak<=0:raise RuntimeError('GPU memory observation unavailable')
        report['observed_peak_GPU_MiB']=memory.peak
        if pilot and memory.peak>plan['maximum_memory_pilot_MiB']:
            raise RuntimeError('Large-model memory pilot exceeds admitted safety budget')
        verify(plan['producer_sha256']);verify(plan['input_sha256']);verify(task['input_sha256'])
        if pilot:report.update(status='memory_pilot_passed',qualified_scope='Neutral RKS energies and gradients up to this case; no Hessian qualification',
                              output_sha256={str(directory/'energy_and_gradient.json'):sha(directory/'energy_and_gradient.json')})
        else:
            report.update(status='state_records_verified',output_sha256={str(p):sha(p) for p in sd.rglob('*') if p.is_file() and p.suffix in ('.json','.xyz')})
    except BaseException as exc:
        report.update(status='failed_preserved',error=repr(exc),traceback=traceback.format_exc());raise
    finally:
        if memory is not None:report['observed_peak_GPU_MiB']=memory.peak
        report['finished']=time.time();save()


def worker(role):
    from redox.workflows.resources import verify_gate_owner
    from redox.workflows.checkpoints import verify_state_checkpoint
    from redox.workflows.processes import stop_owned_children
    from redox.qm.dft import _free_gpu_pool
    verify_gate_owner()
    assert len(os.sched_getaffinity(0))==2 and all(os.environ.get(k)=='2' for k in THREADS)
    plan=json.loads((BASE/'plan.json').read_text());verify(plan['producer_sha256']);verify(plan['input_sha256'])
    owner=dict(host=socket.gethostname(),pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
               gpu=int(os.environ['CUDA_VISIBLE_DEVICES']),cpus=sorted(os.sched_getaffinity(0)),role=role)
    path=BASE/'workers'/f"{owner['host']}_gpu{owner['gpu']}_{role}.json"
    report=dict(owner,status='starting',items=[],started=time.time());atomic_write_json(path,report)
    def stop(signum,frame):raise KeyboardInterrupt(f'signal {signum}')
    signal.signal(signal.SIGTERM,stop)
    consecutive=0
    try:
        if role=='pilot':
            task=prepared_task(next(t for t in plan['tasks'] if t['key']==plan['pilot_key']))
            if task is None:raise RuntimeError('Memory-pilot upstream geometry not ready')
            report.update(status='computing',active_key=task['key']);atomic_write_json(path,report)
            calculate(task,plan,pilot=True)
        else:
            if role=='queue_large':
                pilot=json.loads((BASE/'pilot/status.json').read_text())
                assert pilot['status']=='memory_pilot_passed'
                verify(pilot['output_sha256'])
            for task in plan['tasks']:
                key=task['key'];done=BASE/'done'/f'{key}.json'
                if role=='queue_small' and task['atoms']>plan['small_RKS_geometry_SP_max_atoms']:continue
                ready=prepared_task(task)
                if ready is None:continue
                task=ready
                if done.exists():verify_state_checkpoint(done);continue
                if (BASE/'failures'/f'{key}.json').exists():continue
                claim=BASE/'claims'/key
                try:claim.mkdir()
                except FileExistsError:continue
                atomic_write_json(claim/'owner.json',owner)
                try:
                    if (BASE/'tasks'/key).exists():raise RuntimeError('Prior incomplete attempt requires owner review')
                    report.update(status='computing',active_key=key);atomic_write_json(path,report)
                    print('Actual electronic-structure start',key,flush=True)
                    try:calculate(task,plan)
                    except Exception as exc:
                        if not (BASE/'tasks'/key/'status.json').exists():raise
                        atomic_write_json(BASE/'failures'/f'{key}.json',dict(key=key,owner=owner,error=repr(exc),attempts=1,finished=time.time()))
                        consecutive+=1;report['items'].append(dict(key=key,status='failed_preserved'))
                    else:
                        status=BASE/'tasks'/key/'status.json'
                        atomic_write_json(done,dict(key=key,owner=owner,status_path=str(status),status_sha256=sha(status),finished=time.time()))
                        verify_state_checkpoint(done);consecutive=0;report['items'].append(dict(key=key,status='state_records_verified'))
                finally:
                    stop_owned_children();gc.collect();_free_gpu_pool()
                    (claim/'owner.json').unlink();claim.rmdir();atomic_write_json(path,report)
                env=dict(os.environ,**{k:'1' for k in THREADS})
                with (BASE/'collector.log').open('a') as log:
                    proc=subprocess.run([sys.executable,str(ROOT/'scripts/pipeline/collect_pi_stacking.py')],env=env,stdout=log,stderr=subprocess.STDOUT)
                if proc.returncode:print('Collector needs review; completed records preserved',flush=True)
                if consecutive>=3:raise RuntimeError('Three consecutive actual failures; worker stopped')
        report.update(status='no_unattempted_suitable_work',active_key=None)
    except BaseException as exc:report.update(status='stopped_preserving_outputs',error=repr(exc));raise
    finally:stop_owned_children();report['finished']=time.time();atomic_write_json(path,report)


def launch(role):
    from redox.workflows.resources import allocation_guard,free_slots,gpu_suitability,host_capacity,cpu_allocation,GATE
    host=socket.gethostname().split('.')[0]
    preflight=json.loads((BASE/f'{host}_preflight.json').read_text());assert preflight['ok'] and preflight['python']==sys.executable
    plan=json.loads((BASE/'plan.json').read_text());verify(plan['producer_sha256'])
    record=BASE/'launches'/f'{host}_{role}.json'
    previous=json.loads(record.read_text()) if record.exists() else None
    if previous and role!='queue_large':raise RuntimeError('Inspect existing launch before restart')
    import site
    with allocation_guard():
        suitable=[]
        for _,gpu in free_slots([host]):
            if previous and any(w['gpu']==gpu for w in previous['workers']):continue
            try:device=gpu_suitability(gpu)
            except RuntimeError:continue
            large=device['memory_MiB']>=80*1024
            if (role!='queue_small' and not large) or (role=='queue_small' and large):continue
            suitable.append((gpu,device))
        if role=='pilot':suitable=suitable[:1]  # Explicit dependency: memory pilot before large-model fan-out.
        if not suitable:return
        resources=host_capacity(2*len(suitable),minimum_ram_gib=32*len(suitable));cpus=cpu_allocation(2*len(suitable))
        env=dict(os.environ,PYTHONPATH=str(ROOT/'src'));env.update({k:'2' for k in THREADS})
        env['PATH']=str(PYTHON.parent)+':'+env.get('PATH','')
        env['LD_LIBRARY_PATH']=':'.join([str(p) for b in site.getsitepackages() for p in Path(b).glob('nvidia/*/lib')]+[env.get('LD_LIBRARY_PATH','')])
        report=previous or dict(host=host,role=role,status='submitting',created=time.time(),resources=resources,workers=[],plan_sha256=sha(BASE/'plan.json'))
        report['latest_resources']=resources
        atomic_write_json(record,report)
        for i,(gpu,device) in enumerate(suitable):
            assigned=cpus[2*i:2*i+2];directory=BASE/f'reserved_{host}_gpu{gpu}_{role}';directory.mkdir()
            argv=['taskset','--cpu-list',','.join(map(str,assigned)),str(GATE),'run',str(gpu),'--label','redox-pi-stacking','--',str(PYTHON),'-m','redox.workflows.gpu','--out',str(directory),'--probe','--',str(PYTHON),'-u',str(SCRIPT),'worker','--role',role]
            log=BASE/f'controller_{host}_gpu{gpu}_{role}.log'
            with log.open('w') as stream:
                proc=subprocess.Popen(argv,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True,preexec_fn=lambda assigned=assigned:os.sched_setaffinity(0,assigned))
            report['workers'].append(dict(gpu=gpu,device=device,cpus=assigned,gate_pid=proc.pid,process_start=Path(f'/proc/{proc.pid}/stat').read_text().rsplit(')',1)[1].split()[19],argv=argv,log=str(log)))
            atomic_write_json(record,report);time.sleep(.5)
        report['status']='submitted_actual_starts_pending';atomic_write_json(record,report)
        print(host,role,len(suitable),'reserved workers submitted')


def coordinate():
    """Bounded resource waiter: never preempt, reclaim claims, or retry failed items."""
    from redox.workflows.resources import free_slots,allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    record=BASE/'coordinator.json';assert not record.exists(),'Original coordinator ownership needs review'
    plan=json.loads((BASE/'plan.json').read_text())
    owner=dict(host=socket.gethostname(),pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],started=time.time(),deadline=time.time()+4*3600)
    atomic_write_json(record,dict(owner,status='waiting_for_suitable_free_resources'))
    def remote(host,role):
        env=['env','PYTHONPATH='+str(ROOT/'src'),*[k+'=1' for k in THREADS],str(PYTHON),str(SCRIPT),'launch','--role',role]
        subprocess.run(['ssh','-n','-o','BatchMode=yes','-o','ConnectTimeout=8',host,shlex.join(env)],check=True,timeout=90)
    try:
        while time.time()<owner['deadline']:
            verify(plan['producer_sha256']);verify(plan['input_sha256'])
            stage_inputs(plan)
            done=len(list((BASE/'done').glob('*.json')));failed=len(list((BASE/'failures').glob('*.json')))
            if done+failed==len(plan['tasks']):
                atomic_write_json(record,dict(owner,status='finite_queue_terminal',completed=done,failed=failed));return
            slots=free_slots(['lambda2','lambda11','lambda13'])
            for host in ('lambda2','lambda11'):
                if any(h==host for h,_ in slots) and not (BASE/'launches'/f'{host}_queue_small.json').exists():
                    remote(host,'queue_small')
            pilot=BASE/'pilot/status.json'
            if not (BASE/'launches/lambda13_pilot.json').exists() and any(h=='lambda13' for h,_ in slots) and (BASE/'live_inputs'/plan['pilot_key']/'task.json').exists():remote('lambda13','pilot')
            if pilot.exists():
                state=json.loads(pilot.read_text())['status']
                if state=='failed_preserved':atomic_write_json(record,dict(owner,status='large_model_memory_dependency_failed'));return
                receipt=BASE/'launches/lambda13_queue_large.json'
                launched={w['gpu'] for w in json.loads(receipt.read_text())['workers']} if receipt.exists() else set()
                if state=='memory_pilot_passed' and any(h=='lambda13' and gpu not in launched for h,gpu in slots):remote('lambda13','queue_large')
            atomic_write_json(record,dict(owner,status='finite_campaign_active_or_resource_deferred',completed=done,failed=failed,
                remaining=len(plan['tasks'])-done-failed,memory_pilot=json.loads(pilot.read_text())['status'] if pilot.exists() else 'awaiting_free_H100'))
            time.sleep(15)
        atomic_write_json(record,dict(owner,status='resource_wait_window_ended_remaining_outputs_preserved'))
    except BaseException as exc:atomic_write_json(record,dict(owner,status='stopped_preserving_outputs',error=repr(exc)));raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['freeze','worker','launch','coordinate'])
    parser.add_argument('--role',choices=['pilot','queue_small','queue_large'])
    args=parser.parse_args()
    globals()[args.action](args.role) if args.action in ('worker','launch') else globals()[args.action]()
