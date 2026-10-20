#!/usr/bin/env python
"""Separate frozen intersegment attraction from deformation using intact fragments.

Six small, fixed-geometry SPs for the three currently checked segment pairs.
Same production SMD(MeCN) electronic protocol; no new capping atoms, optimization,
BSSE correction, free energy or electronic-coupling claim.
"""
import argparse
from contextlib import redirect_stdout,redirect_stderr
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import socket
import subprocess
import sys
import time
import traceback

from redox.core.common import ROOT,EV_KJ
from redox.core.io import atomic_write_json
from redox.workflows.checkpoints import verify_state_checkpoint

BASE=ROOT/'results/polymer_models/campaigns/2026-10-20_frozen_fragments'
SCRIPT=Path(__file__).resolve()
THREADS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify(mapping):
    for path,expected in mapping.items():
        if sha(path)!=expected:raise RuntimeError('Changed input/source: '+path)


def decompose(pair,monomer,a,b):
    interaction=(pair-a-b)*EV_KJ
    deformation=(a+b-2*monomer)*EV_KJ
    association=(pair-2*monomer)*EV_KJ
    assert abs(interaction+deformation-association)<1e-6
    return dict(frozen_interaction_deltaE_kJ_mol=interaction,
                segment_deformation_deltaE_kJ_mol=deformation,
                relaxed_association_deltaE_kJ_mol=association,
                cycle_closure_error_kJ_mol=interaction+deformation-association)


def prepare():
    from ase.io import read,write
    from rdkit import Chem
    from rdkit.Chem import rdDetermineBonds
    from redox.validation.protonation import hydrogen_topology_status
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    assert not (BASE/'plan.json').exists(),'Preserve the existing finite fragment campaign'
    for folder in ('claims','done','tasks','failures','workers','inputs'):(BASE/folder).mkdir(parents=True,exist_ok=True)
    feed=json.loads((ROOT/'results/polymer_models/interaction_checked.json').read_text())
    monomers={r['candidate_id']:r for r in feed['models'] if r['repeats']==1 and r['E_smd_eV'] is not None}
    pairs=[r for r in feed['models'] if r['model_type']=='two_nonbonded_segments' and r['E_smd_eV'] is not None]
    assert len(pairs)==3,'Review expansion beyond the original three-pair scope'
    tasks=[];panel=[];sources={}
    for i,row in enumerate(pairs):
        mono=monomers[row['candidate_id']];atoms=read(row['final_xyz'])
        observed=Chem.MolFromXYZFile(row['final_xyz']);rdDetermineBonds.DetermineConnectivity(observed,charge=0)
        fragments=Chem.GetMolFrags(observed);assert len(fragments)==2
        keys=[]
        for j,indices in enumerate(fragments):
            key=f'frozen_{i}_{j}';path=BASE/'inputs'/f'{key}.xyz';part=atoms[list(indices)]
            assert len(part) in (35,45),'Only separately validated small monomer classes are admitted'
            write(path,part,format='xyz');assert hydrogen_topology_status(mono['smiles'],path)=='ok'
            task=dict(key=key,id=key,atoms=len(part),smiles=mono['smiles'],charge=0,mult=1,
                      source_pair_key=row['key'],source_atom_indices=list(indices),input_xyz=str(path),
                      input_sha256={str(path):sha(path),row['final_xyz']:sha(row['final_xyz'])})
            tasks.append(task);keys.append(key)
        panel.append(dict(pair_key=row['key'],candidate_id=row['candidate_id'],orientation=row['target'],
                          fragment_keys=keys,pair_E_eV=row['E_smd_eV'],relaxed_monomer_E_eV=mono['E_smd_eV']))
        for model in (row,mono):
            campaign=ROOT/(model.get('recovery_campaign') or feed['status']['campaign'])
            receipt=campaign/'done'/f"{model['key']}.json";report=verify_state_checkpoint(receipt)
            assert abs(report['E_smd_eV']-model['E_smd_eV'])<1e-9
            sources.update(report['output_sha256']);sources[str(receipt)]=sha(receipt)
    snapshot=BASE/'reference_snapshot.json';atomic_write_json(snapshot,dict(panel=panel,source_sha256=sources))
    files=list((ROOT/'src/redox').rglob('*.py'))+list((ROOT/'config').glob('*.py'))+list((ROOT/'config').glob('*.json'))+[SCRIPT]
    atomic_write_json(BASE/'plan.json',dict(tasks=tasks,created=time.time(),solvent='acetonitrile',
        protocol='Active wB97M-V/def2-TZVPD/SMD electronic single points at fixed fragment coordinates',
        input_sha256={str(snapshot):sha(snapshot)},producer_sha256={str(p):sha(p) for p in files},
        raw_root=str(ROOT/'calcs/polymer_models'/BASE.name),attempts_per_item=1,
        memory_qualification='35-atom phthalimide and 45-atom PMDI monomer SPs previously completed at 13,568/24,224 MiB observed peaks; admit 32-GB GPUs only',
        scope='Full intact segment interaction/deformation; not pure pendant-core interaction, CP/BSSE, deltaG or rates'))
    print('Prepared six frozen-fragment single points; all original atoms retained without new caps')


def collect():
    import fcntl
    with (BASE/'collection.lock').open('a+') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return
        plan=json.loads((BASE/'plan.json').read_text());verify(plan['input_sha256'])
        snapshot=json.loads((BASE/'reference_snapshot.json').read_text());verify(snapshot['source_sha256'])
        energies={};invalid=[]
        for path in (BASE/'done').glob('*.json'):
            try:
                record=verify_state_checkpoint(path);energies[record['task']['key']]=record['E_smd_eV']
            except (OSError,KeyError,ValueError) as exc:invalid.append(dict(receipt=str(path),error=repr(exc)))
        rows=[]
        for pair in snapshot['panel']:
            keys=pair['fragment_keys'];ready=all(k in energies for k in keys)
            values=decompose(pair['pair_E_eV'],pair['relaxed_monomer_E_eV'],*(energies[k] for k in keys)) if ready else dict(frozen_interaction_deltaE_kJ_mol=None,segment_deformation_deltaE_kJ_mol=None,relaxed_association_deltaE_kJ_mol=(pair['pair_E_eV']-2*pair['relaxed_monomer_E_eV'])*EV_KJ)
            rows.append(dict(pair,**values,status='checked_electronic_decomposition' if ready else 'pending_fragment_SP'))
        status=dict(updated=time.time(),qualified_fragment_SP=len(energies),expected_fragment_SP=6,
                    completed_pair_decompositions=sum(r['status']=='checked_electronic_decomposition' for r in rows),
                    failures=len(list((BASE/'failures').glob('*.json'))),invalid_receipts=invalid,
                    active_claims=[p.name for p in (BASE/'claims').iterdir() if p.is_dir()],
                    raw_root=plan['raw_root'],scope=plan['scope'],solvent='acetonitrile')
        atomic_write_json(ROOT/'results/polymer_models/interaction_decomposition.json',dict(status=status,pairs=rows,
            definitions=dict(interaction='E_AB(R_AB) - E_A(R_A frozen) - E_B(R_B frozen)',
                             deformation='E_A(R_A frozen) + E_B(R_B frozen) - 2 E_monomer(relaxed)',
                             closure='relaxed association = frozen interaction + deformation'),
            BSSE_corrected=False,free_energy_computed=False,electron_transfer_computed=False,
            note='SMD fragment cavity changes are included; this is an operational full-segment energy decomposition, not a unique solvent/pure-core energy decomposition'))
        print(json.dumps(status))


def worker():
    from redox.workflows.resources import verify_gate_owner
    from redox.workflows.processes import stop_owned_children
    from redox.qm.sp import ensure
    from redox.core.protocol import load_record
    from redox.qm.dft import _free_gpu_pool
    import pyscf.lib
    import shutil
    verify_gate_owner();pyscf.lib.num_threads(2)
    assert len(os.sched_getaffinity(0))==2 and all(os.environ.get(k)=='2' for k in THREADS)
    plan=json.loads((BASE/'plan.json').read_text());verify(plan['producer_sha256']);verify(plan['input_sha256'])
    owner=dict(host=socket.gethostname().split('.')[0],pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],gpu=int(os.environ['CUDA_VISIBLE_DEVICES']),cpus=sorted(os.sched_getaffinity(0)))
    report=dict(owner,status='claiming',items=[],started=time.time());path=BASE/'workers'/f"{owner['host']}_gpu{owner['gpu']}.json"
    def stop(signum,frame):raise KeyboardInterrupt(f'signal {signum}')
    signal.signal(signal.SIGTERM,stop);consecutive=0
    try:
        for task in plan['tasks']:
            key=task['key'];receipt=BASE/'done'/f'{key}.json'
            if receipt.exists():verify_state_checkpoint(receipt);continue
            if (BASE/'failures'/f'{key}.json').exists():continue
            claim=BASE/'claims'/key
            try:claim.mkdir()
            except FileExistsError:continue
            atomic_write_json(claim/'owner.json',owner)
            try:
                verify(task['input_sha256']);dest=BASE/'tasks'/key;dest.mkdir(exist_ok=False)
                report.update(status='computing',active_key=key);atomic_write_json(path,report)
                state=dict(status='computing',task=task,owner=owner,started=time.time(),geometry_optimization=False,thermal_computed=False)
                status_path=dest/'status.json';atomic_write_json(status_path,state)
                print('Actual frozen-fragment SP start',key,flush=True)
                try:
                    sd=Path(plan['raw_root'])/key/'rest';sd.mkdir(parents=True,exist_ok=False)
                    shutil.copyfile(task['input_xyz'],sd/'opt.xyz')
                    atomic_write_json(sd/'result.json',dict(id=key,state='rest',charge=0,mult=1,solvent='acetonitrile',optimized=False,geometry_definition='Intact segment frozen at checked pair coordinates'))
                    with (dest/'calculation.log').open('w') as log,redirect_stdout(log),redirect_stderr(log):ensure(key,'rest',sd,0,1,backend='gpu')
                    record=load_record(sd,0,1)
                    if record is None or not math.isfinite(record['e_smd_eV']):raise RuntimeError('SP energy not admitted')
                    verify(plan['producer_sha256']);verify(task['input_sha256'])
                    state.update(status='state_records_verified',E_smd_eV=record['e_smd_eV'],state_dir=str(sd),finished=time.time(),output_sha256={str(p):sha(p) for p in sd.rglob('*') if p.is_file() and p.suffix in ('.xyz','.json')})
                    atomic_write_json(status_path,state);atomic_write_json(receipt,dict(key=key,status_path=str(status_path),status_sha256=sha(status_path),owner=owner));verify_state_checkpoint(receipt)
                    report['items'].append(dict(key=key,status='state_records_verified'));consecutive=0
                except Exception as exc:
                    state.update(status='failed_preserved',error=repr(exc),traceback=traceback.format_exc(),finished=time.time());atomic_write_json(status_path,state)
                    atomic_write_json(BASE/'failures'/f'{key}.json',dict(key=key,owner=owner,error=repr(exc)));report['items'].append(dict(key=key,status='failed_preserved'));consecutive+=1
            finally:
                stop_owned_children();_free_gpu_pool();(claim/'owner.json').unlink();claim.rmdir();atomic_write_json(path,report)
            collect()
            if consecutive>=3:raise RuntimeError('Three consecutive actual failures')
        report.update(status='finite_queue_no_unattempted_work',active_key=None)
    finally:
        stop_owned_children();report['finished']=time.time();atomic_write_json(path,report)


def launch():
    from redox.workflows.resources import allocation_guard,free_slots,gpu_suitability,host_capacity,cpu_allocation,GATE
    import site
    host=socket.gethostname().split('.')[0];plan=json.loads((BASE/'plan.json').read_text())
    preflight=json.loads((BASE/f'{host}_preflight.json').read_text());assert preflight['ok'] and preflight['python']==sys.executable
    verify(plan['producer_sha256']);verify(plan['input_sha256']);receipt=BASE/f'launch_{host}.json';assert not receipt.exists()
    remaining=[t for t in plan['tasks'] if not (BASE/'done'/f"{t['key']}.json").exists() and not (BASE/'failures'/f"{t['key']}.json").exists()]
    with allocation_guard():
        suitable=[]
        for _,gpu in free_slots([host]):
            try:suitable.append((gpu,gpu_suitability(gpu)))
            except RuntimeError:continue
        suitable=suitable[:len(remaining)]
        if not suitable:raise RuntimeError('No suitable unoccupied GPU for pending fragment work')
        capacity=host_capacity(2*len(suitable),minimum_ram_gib=32*len(suitable));cpus=cpu_allocation(2*len(suitable))
        env=dict(os.environ,PYTHONPATH=str(ROOT/'src'));env.update({k:'2' for k in THREADS})
        env['LD_LIBRARY_PATH']=':'.join([str(p) for b in site.getsitepackages() for p in Path(b).glob('nvidia/*/lib')]+[env.get('LD_LIBRARY_PATH','')])
        record=dict(host=host,created=time.time(),resources=capacity,plan_sha256=sha(BASE/'plan.json'),workers=[])
        for i,(gpu,device) in enumerate(suitable):
            assigned=cpus[2*i:2*i+2];destination=BASE/f'reserved_{host}_gpu{gpu}';destination.mkdir()
            argv=['taskset','--cpu-list',','.join(map(str,assigned)),str(GATE),'run',str(gpu),'--label','redox-frozen-fragments','--',sys.executable,'-m','redox.workflows.gpu','--out',str(destination),'--probe','--',sys.executable,'-u',str(SCRIPT),'worker']
            with (BASE/f'controller_{host}_gpu{gpu}.log').open('w') as stream:
                proc=subprocess.Popen(argv,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,stdout=stream,stderr=subprocess.STDOUT,start_new_session=True,preexec_fn=lambda assigned=assigned:os.sched_setaffinity(0,assigned))
            record['workers'].append(dict(gpu=gpu,device=device,cpus=assigned,gate_pid=proc.pid,argv=argv));atomic_write_json(receipt,record);time.sleep(.5)
        print('Submitted',len(suitable),'reserved frozen-fragment workers on',host)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['prepare','launch','worker','collect']);args=parser.parse_args();globals()[args.action]()
