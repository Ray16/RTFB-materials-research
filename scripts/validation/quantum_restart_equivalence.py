#!/usr/bin/env python
"""Real interrupted/resumed quantum controls with unchanged physical admission gates."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import shlex
import socket
import subprocess
import sys
import time

import numpy as np

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

SCRIPT = Path(__file__).resolve()
CAPS = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS',
        'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'XTB_THREADS')
CASES = {
    'formaldehyde': dict(smiles='C=O', mult=1,
        xyz='4\nDistorted rigid singlet method control\nC 0 0 0\nO 1.35 0 0\nH -0.68 0.90 0.10\nH -0.65 -0.85 -0.10\n'),
    'methyl_radical': dict(smiles='[CH3]', mult=2,
        xyz='4\nPyramidal radical method control\nC 0 0 0.20\nH 1.15 0 -0.05\nH -0.58 1.00 -0.05\nH -0.58 -1.00 -0.05\n'),
}
OPTIMIZER = dict(convergence_energy=1e-8, convergence_grms=3e-5,
    convergence_gmax=5e-5, convergence_drms=1e-4, convergence_dmax=2e-4)
TOLERANCES = dict(energy_eV=1e-5, free_energy_eV=1e-4,
                  aligned_RMSD_Angstrom=1e-3, frequency_cm=0.5)


def producer_sources():
    sources = [SCRIPT] + [ROOT / 'src/redox' / p for p in (
        'qm/dft.py', 'qm/sp.py', 'qm/native_thermal.py', 'qm/projected_thermal.py',
        'qm/qrrho.py', 'qm/xtb_kernel.py', 'qm/hessian_memory.py', 'core/common.py',
        'core/protocol.py', 'workflows/geometry_checkpoints.py', 'validation/readiness.py')]
    sources += list((ROOT/'src/redox/qm').glob('*.py'))
    sources += [ROOT/'src/redox'/p for p in ('core/io.py','build/geometry.py',
        'validation/protonation.py','validation/spin_ensemble.py')]
    sources += [ROOT/'config/project.json',ROOT/'config/electrolyte.py',
                *list((ROOT/'config').glob('native_*.json')),
                ROOT/'config/projected_thermal_validation.json']
    return sources


def prepare(campaign):
    from redox.validation.protonation import hydrogen_topology_status
    campaign.mkdir(parents=True, exist_ok=False)
    tasks = []
    for name, case in CASES.items():
        path = campaign / (name + '.xyz'); path.write_text(case['xyz'])
        if hydrogen_topology_status(case['smiles'], path) != 'ok':
            raise ValueError('Control seed does not have the exact declared chemistry')
        for phase in ('baseline', 'interrupted', 'resumed'):
            tasks.append(dict(key=name+'__'+phase, case=name, phase=phase,
                charge=0, mult=case['mult'], smiles=case['smiles'], input=str(path),
                depends_on=name+'__interrupted' if phase=='resumed' else None))
    sources = producer_sources()
    atomic_write_json(campaign/'plan.json', dict(tasks=tasks,
        raw_root=str(ROOT/'calcs/quantum_restart'/campaign.name),
        role='Singlet and radical execution method controls; no candidate registration or synthesis admission',
        optimizer=OPTIMIZER, tolerances=TOLERANCES, interrupt_after_evaluations=3,
        expected_interrupt_returncode=87, native_memory_cap_GiB=6,
        minimum_GPU_MiB=31744, actual_attempts_per_task=1,
        producer_sha256={str(p):api.digest(p) for p in sorted(set(sources))},
        input_sha256={str(p):api.digest(p) for p in campaign.glob('*.xyz')},
        interpretation='Checkpoint coordinates remain unqualified; resumed optimization and all native gates must pass'))


def task_folder(campaign, key):
    plan = api.read_plan(campaign)
    return Path(plan['raw_root']) / key


def run_task(campaign, key):
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner()
    if len(os.sched_getaffinity(0)) != 2 or any(os.environ.get(k) != '2' for k in CAPS):
        raise RuntimeError('Two disjoint CPU slots and native caps=2 required')
    plan = api.read_plan(campaign)
    api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])
    task = next(t for t in plan['tasks'] if t['key']==key)
    folder = task_folder(campaign, key); folder.mkdir(parents=True, exist_ok=False)
    owner = dict(host=socket.gethostname(), pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
        cpus=sorted(os.sched_getaffinity(0)), gpu=os.environ['CUDA_VISIBLE_DEVICES'])
    status = dict(status='computing', task=task, owner=owner, started=time.time())
    atomic_write_json(folder/'status.json', status)
    source = Path(task['input'])
    if task['phase']=='resumed':
        proof = json.loads((campaign/'death_checks'/str(task['depends_on']+'.json')).read_text())
        if proof['owner_dead'] is not True or proof['expected_interrupt_returncode']!=87:
            raise ValueError('The intentionally interrupted owner has not been proved dead')
        source = task_folder(campaign, task['depends_on'])/'geometry_checkpoints/last_evaluated.xyz'
        if api.digest(source)!=proof['checkpoint_sha256']:
            raise ValueError('Interrupted geometry changed after ownership audit')
    from redox.qm import dft, native_thermal as native
    from redox.qm.sp import ensure
    from redox.core.protocol import load_record
    from redox.core.common import apply_thermal_policy, free_energy
    from redox.validation.protonation import hydrogen_topology_status
    from redox.workflows.geometry_checkpoints import checkpoint_callback
    previous_manifest = native.read_manifest
    row = dict(id=task['case'],state='rest',charge=0,mult_hint=task['mult'],
               n_e=0,smiles=task['smiles'])
    sd = folder/task['case']/'rest'; sd.mkdir(parents=True)
    save = checkpoint_callback(folder/'geometry_checkpoints')
    count = 0

    def callback(environment):
        nonlocal count
        save(environment); count += 1
        if task['phase']=='interrupted' and count==plan['interrupt_after_evaluations']:
            status.update(status='intentional_process_exit_after_durable_evaluation',
                          checkpoint_sha256=api.digest(folder/'geometry_checkpoints/last_evaluated.xyz'),
                          physical_admission=False, evaluations=count, finished=time.time())
            atomic_write_json(folder/'status.json',status)
            os._exit(87)

    try:
        native.read_manifest = lambda: [row]
        result = dft.dft_smd(source,0,task['mult'],backend='gpu',do_freq=False,
            opt_out=sd/'opt.xyz',conv_params=dict(plan['optimizer'],callback=callback),max_opt_steps=120)
        if task['phase']=='interrupted':
            raise RuntimeError('Control optimized before the requested real interruption')
        result.update(id=task['case'],state='rest',n_e=0)
        atomic_write_json(sd/'result.json',result)
        if hydrogen_topology_status(task['smiles'],sd/'opt.xyz')!='ok':
            raise ValueError('Quantum optimization changed exact chemical identity')
        ensure(task['case'],'rest',sd,0,task['mult'],backend='gpu')
        native.ensure_native_thermal(task['case'],'rest',state_dir=sd,memory_cap_gib=6)
        native.select_policy(sd)
        energy=load_record(sd,0,task['mult']); thermal=native.load_native_thermal(sd,0,task['mult'])
        if energy is None or thermal is None:
            raise ValueError('Active energy or normal native admission failed')
        G,qc=free_energy(apply_thermal_policy(dict(result,e_smd_eV=energy['e_smd_eV']),sd))
        if G is None or qc!='ok':
            raise ValueError('Normal free-energy QC failed')
        api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])
        status.update(status='quantum_control_verified',G_eV=G,state_dir=str(sd),
            output_sha256={str(p):api.digest(p) for p in sd.rglob('*') if p.is_file()})
    except BaseException as error:
        status.update(status='failed_preserved',error=repr(error)); raise
    finally:
        native.read_manifest=previous_manifest
        status['finished']=time.time();atomic_write_json(folder/'status.json',status)


def restart_differences(a,b):
    if a[0].get_chemical_symbols()!=b[0].get_chemical_symbols():
        raise ValueError('Restart atom order changed')
    x=a[0].positions-a[0].positions.mean(axis=0);y=b[0].positions-b[0].positions.mean(axis=0)
    u,_,vt=np.linalg.svd(x.T@y);rotation=u@np.diag([1,1,np.linalg.det(u@vt)])@vt
    rms=float(np.sqrt(np.mean(np.sum((x@rotation-y)**2,axis=1))))
    frequency=float(np.max(abs(np.sort(a[2]['frequencies_cm'])-np.sort(b[2]['frequencies_cm']))))
    difference=dict(energy_eV=max(abs(a[1][k]-b[1][k]) for k in ('e_smd_eV','e_gas_eV')),
        free_energy_eV=abs(a[3]['G_eV']-b[3]['G_eV']),aligned_RMSD_Angstrom=rms,frequency_cm=frequency)
    if any(not np.isfinite(v) or v>TOLERANCES[k] for k,v in difference.items()):
        raise ValueError('Resumed quantum output differs beyond the predeclared tolerance: '+repr(difference))
    return difference


def audit(campaign):
    from ase.io import read
    from redox.core.protocol import load_record
    from redox.qm.native_thermal import load_native_thermal
    plan=api.read_plan(campaign); results=[]
    for name,case in CASES.items():
        records=[]
        for phase in ('baseline','resumed'):
            folder=task_folder(campaign,name+'__'+phase)
            status=json.loads((folder/'status.json').read_text())
            if status['status']!='quantum_control_verified':raise ValueError('Both real quantum controls must finish')
            api.verify_sources(status['output_sha256'])
            sd=Path(status['state_dir']); energy=load_record(sd,0,case['mult']);thermal=load_native_thermal(sd,0,case['mult'])
            if energy is None or thermal is None:raise ValueError('Backing energy/native records no longer qualify')
            expected_G=energy['e_smd_eV']+thermal['g_thermal_total_eV']
            if abs(expected_G-status['G_eV'])>1e-9:
                raise ValueError('Stored free energy does not reconstruct from checked backing records')
            records.append((read(sd/'opt.xyz'),energy,thermal,status))
        if plan['tolerances']!=TOLERANCES:raise ValueError('Predeclared equivalence thresholds changed')
        difference=restart_differences(*records)
        results.append(dict(case=name,mult=case['mult'],differences=difference,
            owner_death_check=str(campaign/'death_checks'/(name+'__interrupted.json'))))
    result=dict(status='real_quantum_interruption_equivalence_verified',controls=results,
        tolerances=plan['tolerances'],raw_root=plan['raw_root'],plan_sha256=api.digest(campaign/'plan.json'),
        scientific_model_changed=False,physical_gates_relaxed=False,
        limits='Two small rigid spin controls; not universal recovery, chemical accuracy or flexible-polymer ensemble evidence')
    atomic_write_json(campaign/'summary.json',result);return result


def coordinate(campaign,max_hours):
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity,free_slots
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    with (campaign/'coordinator.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        deadline=time.time()+max_hours*3600; launches={}; terminals=set()
        (campaign/'death_checks').mkdir();(campaign/'launches').mkdir()
        while time.time()<deadline:
            plan=api.read_plan(campaign);api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
            for key,launch in list(launches.items()):
                receipt=Path(launch['execution'])/'execution.json'
                if not receipt.exists():continue
                execution=json.loads(receipt.read_text())
                if execution['status'] not in ('failed','process_completed'):continue
                state=task_folder(campaign,key)/'status.json'
                if not state.exists():
                    if 'occupied or reserved; no calculation started' in execution.get('error',''):
                        del launches[key];continue
                    raise RuntimeError('Worker did not start; inspect preserved wrapper: '+key)
                status=json.loads(state.read_text())
                if status['task']['phase']=='interrupted':
                    if execution.get('returncode')!=87 or status['status']!='intentional_process_exit_after_durable_evaluation':
                        raise RuntimeError('Unexpected quantum interruption outcome')
                    owner=status['owner']
                    script="from pathlib import Path; p=Path('/proc/"+str(owner['pid'])+"/stat'); print(p.read_text().rsplit(')',1)[1].split()[19] if p.exists() else 'absent')"
                    check=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8',owner['host'],
                        shlex.join(['python3','-c',script])],capture_output=True,text=True,check=True,timeout=20)
                    if check.stdout.strip()==owner['process_start']:raise RuntimeError('Interrupted quantum owner is still live')
                    checkpoint=task_folder(campaign,key)/'geometry_checkpoints/last_evaluated.xyz'
                    atomic_write_json(campaign/'death_checks'/(key+'.json'),dict(owner=owner,owner_dead=True,
                        expected_interrupt_returncode=87,checkpoint_sha256=api.digest(checkpoint),execution_sha256=api.digest(receipt)))
                elif execution.get('returncode')!=0 or status['status']!='quantum_control_verified':
                    raise RuntimeError('Actual quantum control failed; preserve and review: '+key)
                terminals.add(key)
            if len(terminals)==len(plan['tasks']):
                audit(campaign)
                atomic_write_json(campaign/'coordinator.json',dict(status='all_quantum_controls_and_restart_comparisons_verified',
                    terminal_tasks=sorted(terminals),finished=time.time()))
                return
            ready=[t for t in plan['tasks'] if t['key'] not in launches and t['key'] not in terminals
                   and (not t['depends_on'] or t['depends_on'] in terminals)]
            available=free_slots() if ready else []
            for host in sorted({h for h,_ in available}):
                p=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8',host,
                    'nvidia-smi --query-gpu=index,memory.total --format=csv,noheader,nounits'],capture_output=True,text=True,timeout=20)
                memory={int(r.split(',')[0]):float(r.split(',')[1]) for r in p.stdout.splitlines()} if p.returncode==0 else {}
                slots=[g for h,g in available if h==host and memory.get(g,0)>=plan['minimum_GPU_MiB']]
                if not slots:continue
                preflight=[sys.executable,str(ROOT.parent/'.codex-personal/skills/lambda-fleet/scripts/preflight_node.py'),host,
                    '--python',sys.executable,'--module','pyscf','--module','gpu4pyscf','--module','ase','--module','rdkit',
                    '--path',str(campaign/'plan.json'),'--min-available-ram-mib',str(8192*min(len(slots),len(ready))),'--timeout','40']
                p=subprocess.run(preflight,capture_output=True,text=True,timeout=50)
                (campaign/(host+'_preflight.json')).write_text(p.stdout)
                if p.returncode:continue
                for gpu in slots:
                    if not ready:break
                    task=ready.pop(0);execution=campaign/'execution'/(task['key']+'_'+str(time.time_ns()))
                    argv=['env','PYTHONPATH='+str(ROOT/'src'),*[k+'=1' for k in CAPS],sys.executable,
                        '-m','redox.workflows.gpu','--gpu',str(gpu),'--threads','2','--out',str(execution),'--',
                        sys.executable,str(SCRIPT),'worker','--campaign',str(campaign),'--key',task['key']]
                    log=campaign/'launches'/(execution.name+'.log')
                    command='nohup '+shlex.join(argv)+' >'+shlex.quote(str(log))+' 2>&1 </dev/null & echo $!'
                    p=subprocess.run(['ssh','-o','BatchMode=yes',host,command],capture_output=True,text=True,check=True,timeout=20)
                    launch=dict(host=host,gpu=gpu,pid=p.stdout.strip(),argv=argv,execution=str(execution),started=time.time())
                    launches[task['key']]=launch;atomic_write_json(campaign/'launches'/(execution.name+'.json'),launch)
            atomic_write_json(campaign/'coordinator.json',dict(status='running_or_waiting_for_suitable_GPUs',
                pid=os.getpid(),host=socket.gethostname(),deadline=deadline,launches=launches,terminal_tasks=sorted(terminals),updated=time.time()))
            time.sleep(10)
        atomic_write_json(campaign/'coordinator.json',dict(status='finite_resource_window_ended_outputs_preserved',terminal_tasks=sorted(terminals)))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=('prepare','worker','coordinate','audit'))
    parser.add_argument('--campaign',type=Path,required=True);parser.add_argument('--key')
    parser.add_argument('--max-hours',type=float,default=24)
    args=parser.parse_args();campaign=args.campaign.resolve()
    if not 0<args.max_hours<=48:parser.error('Finite resource window up to 48 hours required')
    if args.action=='prepare':prepare(campaign)
    elif args.action=='worker':run_task(campaign,args.key)
    elif args.action=='coordinate':
        try:coordinate(campaign,args.max_hours)
        except BaseException as error:
            atomic_write_json(campaign/'coordinator.json',dict(status='coordinator_failed_preserving_actual_workers',
                error=repr(error),finished=time.time()))
            raise
    else:print(json.dumps(audit(campaign)))
