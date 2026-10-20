#!/usr/bin/env python
"""Finite DBU/DIPEA proton-capture pilot; preserve the neutral grafting campaign."""
from __future__ import annotations

import argparse
from contextlib import redirect_stdout
import gc
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import signal
import socket
import subprocess
import sys
import time

from redox.core.common import ROOT, EV_KJ
from redox.core.io import atomic_write_json
from redox.workflows.checkpoints import verify_state_checkpoint

BASE = ROOT/'results/grafting/campaigns/2026-10-19_base_assisted_pilot'
OUT = ROOT/'results/grafting/base_assisted'
SCRIPT = Path(__file__).resolve()
UPSTREAM = ROOT/'results/grafting/campaigns/2026-10-19_mecn_reference'
THREADS = ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS',
           'MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
spec = importlib.util.spec_from_file_location('reference_grafting', ROOT/'scripts/pipeline/grafting_campaign.py')
G = importlib.util.module_from_spec(spec)
spec.loader.exec_module(G)
G.BASE = BASE
G.SCRIPT = SCRIPT

BASES = {
    'dbu': dict(neutral='C1CCC2=NCCCN2CC1', protonated='C1CCC2=[NH+]CCCN2CC1',
                source='https://onlinelibrary.wiley.com/doi/abs/10.1002/pat.1968',
                precedent='AQ carboxylate attachment in NMP; MeCN and exact phenolic designs are proposed'),
    'dipea': dict(neutral='CCN(C(C)C)C(C)C', protonated='CC[NH+](C(C)C)C(C)C',
                  source='https://chemistry-europe.onlinelibrary.wiley.com/doi/10.1002/chem.202200149',
                  precedent='Ferrocene carboxylate attachment to Merrifield beads in DMF; exact phenolic designs are proposed')}
PANEL = ['exp_bq_oxy_methoxy','exp_bq_oxy_trimethyl','exp_aq_1oxy',
         'exp_aq_2oxy_6methoxy','design_ndi_ethyl','design_pmdi_ethyl']


def concentration_correction(delta_n, temperature=298.15):
    """kJ/mol: ideal thermal 1-atm convention to 1 M, once per solute."""
    return delta_n * .00831446261815324 * temperature * math.log(.08205736608096 * temperature)


def corrections(energies, base, hcl, chloride):
    neutral = energies.get(base+'_neutral')
    protonated = energies.get(base+'_protonated')
    pair = energies.get(base+'_chloride_pair')
    free = ((protonated+chloride-neutral-hcl)*EV_KJ
            if neutral is not None and protonated is not None else None)
    paired = ((pair-neutral-hcl)*EV_KJ + concentration_correction(-1)
              if neutral is not None and pair is not None else None)
    association = ((pair-protonated-chloride)*EV_KJ + concentration_correction(-1)
                   if pair is not None and protonated is not None else None)
    return dict(free_ions_kJ_mol=free, contact_pair_kJ_mol=paired,
                ion_pair_association_kJ_mol=association)


def pair_seeds(protonated, directory):
    """Place chloride beyond the declared N-H bond; no overlapping fragment embedding."""
    import numpy as np
    from rdkit import Chem
    from redox.build.structures import conformer_ensemble
    mol, cids, _, _ = conformer_ensemble(Chem.MolFromSmiles(protonated), net_charge=1, top_k=3)
    anchor = next(a for a in mol.GetAtoms() if a.GetSymbol()=='N' and a.GetFormalCharge()==1)
    hydrogen = next(a for a in anchor.GetNeighbors() if a.GetAtomicNum()==1)
    for i,cid in enumerate(cids):
        xyz=np.array(mol.GetConformer(cid).GetPositions())
        axis=xyz[hydrogen.GetIdx()]-xyz[anchor.GetIdx()]
        axis/=np.linalg.norm(axis)
        cl=xyz[hydrogen.GetIdx()]+(2.2+.4*i)*axis
        assert np.linalg.norm(xyz-cl,axis=1).min()>1.5
        lines=[str(len(xyz)+1),'Explicit protonated-base/chloride contact-pair starting pose']
        lines += [f'{a.GetSymbol()} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}' for a,p in zip(mol.GetAtoms(),xyz)]
        lines += [f'Cl {cl[0]:.10f} {cl[1]:.10f} {cl[2]:.10f}']
        (directory/f'conf_{i:02d}.xyz').write_text('\n'.join(lines)+'\n')
    return len(cids)


def prepare():
    from collections import Counter
    from rdkit import Chem
    from redox.build.grafting import atom_inventory
    from redox.build.structures import conformer_ensemble,to_xyz
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    assert not (BASE/'plan.json').exists(), 'Preserve prior campaign'
    with allocation_guard():
        host_capacity(2,minimum_ram_gib=8)
        os.sched_setaffinity(0,cpu_allocation(2))
    for folder in ('claims','done','failures','tasks','workers','workspace/library'):
        (BASE/folder).mkdir(parents=True,exist_ok=True)
    OUT.mkdir(parents=True,exist_ok=True)
    actual=ROOT/'calcs/grafting'/BASE.name
    actual.mkdir(exist_ok=False)
    (BASE/'workspace/calcs').mkdir()
    (BASE/'workspace/calcs/dft').symlink_to(actual,target_is_directory=True)
    refs={}
    for key in ('hydrogen_chloride','chloride','merrifield_local_chloride'):
        receipt=UPSTREAM/'done'/f'{key}.json'
        report=verify_state_checkpoint(receipt)
        refs[key]=dict(G_eV=report['G_eV'],receipt=str(receipt),receipt_sha256=G.digest(receipt))
    cached=json.loads((UPSTREAM/'reused_products.json').read_text())
    roster=json.loads((ROOT/'results/grafting/reactions.json').read_text())['reactions']
    panel=[]
    for reaction in roster:
        if reaction['id'] not in PANEL:continue
        assert reaction['reaction_class']=='NH_OH_substitution'
        prec=UPSTREAM/'done'/f"{reaction['precursor_id']}.json"
        item=dict(id=reaction['id'],name=reaction['name'],family=reaction['family'],
                  reaction=reaction,reference_deltaG_kJ_mol=None,source_sha256={})
        product=cached.get(reaction['id'])
        try:
            report=verify_state_checkpoint(prec)
            if product is None:raise ValueError('No qualified reusable product')
            G.verify(product['source_sha256'])
            value=product['G_eV']+refs['hydrogen_chloride']['G_eV']-report['G_eV']-refs['merrifield_local_chloride']['G_eV']
            item.update(reference_deltaG_kJ_mol=value*EV_KJ,
                        source_sha256={str(prec):G.digest(prec),**report['output_sha256'],**product['source_sha256']})
        except (OSError,ValueError,RuntimeError) as exc:item['reference_hold']=repr(exc)
        panel.append(item)
    assert len(panel)==6
    atomic_write_json(BASE/'reference_snapshot.json',dict(references=refs,panel=panel))
    tasks=[]
    for base,identity in BASES.items():
        assert atom_inventory(Chem.MolFromSmiles(identity['neutral']))+Counter({1:1})==atom_inventory(Chem.MolFromSmiles(identity['protonated']))
        for form in ('neutral','protonated','chloride_pair'):
            key=base+'_'+form
            smiles=identity['neutral'] if form=='neutral' else identity['protonated']+('.[Cl-]' if form=='chloride_pair' else '')
            mol=Chem.MolFromSmiles(smiles);charge=Chem.GetFormalCharge(mol)
            directory=BASE/'workspace/library'/key
            confdir=directory/'conformers';confdir.mkdir(parents=True)
            if form=='chloride_pair':count=pair_seeds(identity['protonated'],confdir)
            else:
                molh,cids,ff,es=conformer_ensemble(mol,net_charge=charge,top_k=3)
                for i,(cid,e) in enumerate(zip(cids,es)):
                    one=Chem.Mol(molh);one.RemoveAllConformers();one.AddConformer(molh.GetConformer(cid),assignId=True)
                    (confdir/f'conf_{i:02d}.xyz').write_text(to_xyz(one,f'{ff} seed E={e}'))
                count=len(cids)
            shutil.copyfile(confdir/'conf_00.xyz',directory/'rest.xyz')
            row=dict(id=key,name=key,family='base_capture_pilot',state='rest',charge=charge,
                     mult_hint=1,n_e=0,n_conf=count,solv_preopt=1,smiles=smiles)
            tasks.append(dict(key=key,id=key,state='rest',role='base_component',atoms=Chem.AddHs(mol).GetNumAtoms(),
                              manifest_row=row,input_sha256=G.hashes(directory),atom=False,preopt='alpb'))
    producers=list((ROOT/'src/redox').rglob('*.py'))+list((ROOT/'config').glob('*.json'))+list((ROOT/'config').glob('*.py'))
    producers += [SCRIPT,ROOT/'scripts/pipeline/grafting_campaign.py',ROOT/'scripts/plotting/grafting/update_base_assisted_figures.py']
    upstream_plan=json.loads((UPSTREAM/'plan.json').read_text())
    plan=dict(schema=1,created=time.time(),tasks=tasks,bases=BASES,solvent='acetonitrile',temperature_K=298.15,
              raw_root=str(actual),attempts_per_item=1,worker_consecutive_failure_limit=3,
              optimizer_stopping_tolerances=upstream_plan['optimizer_stopping_tolerances'],
              producer_sha256={str(p):G.digest(p) for p in producers},
              input_sha256={str(BASE/'reference_snapshot.json'):G.digest(BASE/'reference_snapshot.json')},
              thermal_standard_state='native ideal-gas 1 atm; explicit delta-n correction to 1 M in reaction sums',
              purpose='MeCN thermodynamic sensitivity pilot; no barriers, yields, solid salts or synthesis admission',
              physical_limits=['Approximate contact-pair entropy and incomplete conformer sampling',
                               'Base alkylation and quinone side reactions not calculated',
                               'Actual experimental solvent/temperature may differ from this matched MeCN screen'])
    atomic_write_json(BASE/'plan.json',plan)
    print(json.dumps(dict(prepared_components=len(tasks),reference_panel=len(panel),reference_holds=sum(r['reference_deltaG_kJ_mol'] is None for r in panel))))


def collect():
    import fcntl
    with (OUT/'.collection.lock').open('a+') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:return
        plan=json.loads((BASE/'plan.json').read_text());G.verify(plan['input_sha256'])
        source=json.loads((BASE/'reference_snapshot.json').read_text())
        for ref in source['references'].values():
            assert G.digest(ref['receipt'])==ref['receipt_sha256']
            verify_state_checkpoint(ref['receipt'])
        energies={};invalid=[]
        for p in (BASE/'done').glob('*.json'):
            try:
                r=verify_state_checkpoint(p);energies[r['task']['key']]=r['G_eV']
            except (OSError,ValueError,KeyError) as exc:invalid.append(dict(receipt=str(p),error=repr(exc)))
        shifts={b:corrections(energies,b,source['references']['hydrogen_chloride']['G_eV'],
                              source['references']['chloride']['G_eV']) for b in BASES}
        rows=[]
        for item in source['panel']:
            row=dict(item);reference=row['reference_deltaG_kJ_mol']
            try:G.verify(row['source_sha256'])
            except (OSError,RuntimeError):reference=None;row['reference_hold']='Upstream data changed'
            row['base_models']={b:{k:reference+v if reference is not None and v is not None else None
                                   for k,v in shifts[b].items() if k!='ion_pair_association_kJ_mol'} for b in BASES}
            row['eligible_for_synthesis']=False
            rows.append(row)
        status=dict(updated=time.time(),components_qualified=len(energies),components_expected=6,
                    failures=len(list((BASE/'failures').glob('*.json'))),invalid_receipts=invalid,
                    base_models_ready=sum(s['free_ions_kJ_mol'] is not None for s in shifts.values()),
                    ion_pair_models_ready=sum(s['contact_pair_kJ_mol'] is not None for s in shifts.values()),
                    scope=plan['purpose'],campaign=str(BASE.relative_to(ROOT)),raw_root=plan['raw_root'],
                    active_claims=[p.name for p in (BASE/'claims').iterdir() if p.is_dir()],
                    calculation_start_records=len(list((BASE/'tasks').glob('*/status.json'))),
                    exact_grafting_accuracy_validated=False,synthesis_admission=False)
        atomic_write_json(OUT/'current.json',dict(status=status,corrections=shifts,candidates=rows,
                          bases=BASES,physical_limits=plan['physical_limits'],standard_state=plan['thermal_standard_state']))
        atomic_write_json(OUT/'status.json',status)
        subprocess.run([str(ROOT/'.venv-figures/bin/python'),str(ROOT/'scripts/plotting/grafting/update_base_assisted_figures.py')],check=True)
        print(json.dumps(status))


def worker():
    from redox.workflows.resources import verify_gate_owner
    from redox.workflows.processes import stop_owned_children
    verify_gate_owner()
    plan=json.loads((BASE/'plan.json').read_text());G.verify(plan['producer_sha256']);G.verify(plan['input_sha256'])
    from redox.qm.native_thermal import implementation_controls
    implementation_controls()
    assert len(os.sched_getaffinity(0))==2 and all(os.environ.get(k)=='2' for k in THREADS)
    owner=dict(host=socket.gethostname().split('.')[0],pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
               cpus=sorted(os.sched_getaffinity(0)),gpu=int(os.environ['CUDA_VISIBLE_DEVICES']))
    path=BASE/'workers'/f"{owner['host']}_gpu{owner['gpu']}.json"
    report=dict(owner,status='claiming',started=time.time(),items=[]);consecutive=0
    def stop(signum,frame):raise KeyboardInterrupt(f'signal {signum}')
    signal.signal(signal.SIGTERM,stop)
    try:
        for task in plan['tasks']:
            key=task['key'];done=BASE/'done'/f'{key}.json'
            if done.exists():verify_state_checkpoint(done);continue
            if (BASE/'failures'/f'{key}.json').exists():continue
            claim=BASE/'claims'/key
            try:claim.mkdir()
            except FileExistsError:continue
            atomic_write_json(claim/'owner.json',owner)
            try:
                if done.exists():verify_state_checkpoint(done);continue
                if (BASE/'tasks'/key).exists():raise RuntimeError('Prior attempt needs owner review')
                report.update(status='computing',active_key=key);atomic_write_json(path,report)
                print('Actual calculation start',key,flush=True)
                try:G.run_state(task,plan)
                except Exception as exc:
                    if not (BASE/'tasks'/key/'status.json').exists():raise
                    atomic_write_json(BASE/'failures'/f'{key}.json',dict(key=key,owner=owner,error=repr(exc),attempts=1,finished=time.time()))
                    report['items'].append(dict(key=key,status='failed_preserved'));consecutive+=1
                else:
                    state=BASE/'tasks'/key/'status.json'
                    atomic_write_json(done,dict(key=key,owner=owner,status_path=str(state),status_sha256=G.digest(state),finished=time.time()))
                    verify_state_checkpoint(done);report['items'].append(dict(key=key,status='state_records_verified'));consecutive=0
            finally:
                stop_owned_children();gc.collect()
                from redox.qm.dft import _free_gpu_pool
                _free_gpu_pool();(claim/'owner.json').unlink();claim.rmdir();atomic_write_json(path,report)
            env=dict(os.environ,**{k:'1' for k in THREADS})
            with (BASE/'collector.log').open('a') as stream:
                result=subprocess.run([sys.executable,str(SCRIPT),'collect'],env=env,stdout=stream,stderr=subprocess.STDOUT)
            if result.returncode:print('Collector needs review',flush=True)
            if consecutive>=3:raise RuntimeError('Three consecutive actual failures')
        report.update(status='no_unclaimed_unattempted_tasks',active_key=None)
    except BaseException as exc:
        report.update(status='stopped_preserving_outputs',error=repr(exc));raise
    finally:
        report['finished']=time.time();atomic_write_json(path,report)


def launch():
    from redox.workflows.resources import allocation_guard,cpu_allocation,free_slots,gpu_suitability,host_capacity,GATE
    import site
    plan=json.loads((BASE/'plan.json').read_text());G.verify(plan['producer_sha256']);G.verify(plan['input_sha256'])
    host=socket.gethostname().split('.')[0];preflight=json.loads((BASE/f'{host}_preflight.json').read_text())
    assert preflight['ok'] and preflight['python']==sys.executable
    receipt=BASE/f'launch_{host}.json';assert not receipt.exists(),'Inspect prior launch'
    with allocation_guard():
        remaining=[t for t in plan['tasks'] if not (BASE/'done'/f"{t['key']}.json").exists()
                   and not (BASE/'failures'/f"{t['key']}.json").exists() and not (BASE/'claims'/t['key']).exists()]
        suitable=[]
        for _,gpu in free_slots([host]):
            try:suitable.append((gpu,gpu_suitability(gpu)))
            except RuntimeError:continue
        suitable=suitable[:len(remaining)]
        if not suitable:print('No unclaimed work or suitable free GPU');return
        capacity=host_capacity(2*len(suitable),minimum_ram_gib=32*len(suitable));cpus=cpu_allocation(2*len(suitable))
        libs=[str(p) for b in site.getsitepackages() for p in Path(b).glob('nvidia/*/lib')]
        env=dict(os.environ,PYTHONPATH=str(ROOT/'src'),REDOX_HESSIAN_MEMORY_GIB='24')
        env['PATH']=str(Path(sys.executable).parent)+os.pathsep+env.get('PATH','')
        env['LD_LIBRARY_PATH']=':'.join(libs+[env.get('LD_LIBRARY_PATH','')]);env.update({k:'2' for k in THREADS})
        record=dict(host=host,created=time.time(),status='submitting',resources=capacity,plan_sha256=G.digest(BASE/'plan.json'),workers=[])
        atomic_write_json(receipt,record)
        for i,(gpu,device) in enumerate(suitable):
            assigned=cpus[2*i:2*i+2];reservation=BASE/f'reserved_{host}_gpu{gpu}';reservation.mkdir()
            argv=['taskset','--cpu-list',','.join(map(str,assigned)),str(GATE),'run',str(gpu),'--label','redox-base-pilot','--',
                  sys.executable,'-m','redox.workflows.gpu','--out',str(reservation),'--probe','--',sys.executable,'-u',str(SCRIPT),'worker']
            log=BASE/f'controller_{host}_gpu{gpu}.log'
            with log.open('w') as stream:
                proc=subprocess.Popen(argv,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,stdout=stream,stderr=subprocess.STDOUT,
                                      start_new_session=True,preexec_fn=lambda assigned=assigned:os.sched_setaffinity(0,assigned))
            record['workers'].append(dict(gpu=gpu,device=device,cpus=assigned,gate_pid=proc.pid,
                process_start=Path(f'/proc/{proc.pid}/stat').read_text().rsplit(')',1)[1].split()[19],argv=argv,log=str(log)))
            atomic_write_json(receipt,record);time.sleep(.5)
        record['status']='submitted_actual_start_verification_pending';atomic_write_json(receipt,record)
        print(json.dumps(dict(host=host,submitted_workers=len(suitable))))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['prepare','worker','launch','collect'])
    globals()[parser.parse_args().action]()
