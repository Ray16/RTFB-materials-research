#!/usr/bin/env python
"""Prospective 16-GB-card energy pilots for reviewed 19/38-atom local controls.

Frozen single points only. This is not the large-polymer geometry/Hessian
execution class, and its resource observations cannot qualify that workflow.
"""
import argparse
from contextlib import redirect_stdout,redirect_stderr
import gc
import fcntl
import json
import math
import os
from pathlib import Path
import socket
import subprocess
import sys
import time

from redox.core.common import ROOT,HARTREE_EV
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

PROFILE=ROOT/'environments/cutensor22/runtime.json'
SETTINGS=dict(xc='wb97m-v',basis='def2-tzvpd',nlc='vv10',disp=None,density_fit=True,
    only_dfj=False,grid_level=3,nlc_grid_level=3,SCF_conv_tol=1e-10,SCF_conv_tol_grad=1e-7,
    SCF_max_cycle=200,solvent_model='unchanged GPU4PySCF legacy SMD',SMD_surface_points=590,
    media=['gas','acetonitrile'],geometry_optimized=False,Hessian_computed=False)


def device(gpu):
    r=subprocess.run(['nvidia-smi',f'--id={gpu}','--query-gpu=index,uuid,name,memory.total,compute_cap',
        '--format=csv,noheader,nounits'],capture_output=True,text=True,timeout=15,check=True)
    fields=[s.strip() for s in r.stdout.strip().split(',')]
    if len(fields)!=5 or int(fields[0])!=gpu or float(fields[3])<15360:
        raise RuntimeError('Local-only prospective energy pilot needs at least 15360 MiB')
    return dict(index=gpu,uuid=fields[1],name=fields[2],memory_MiB=float(fields[3]),compute_capability=fields[4])


def prepare(construction,review_path,out):
    witness=json.loads(review_path.read_text());api.verify_sources(witness['source_sha256'])
    if witness['status']!='fixed_parent_local_motif_construction_independently_checked':raise ValueError('Checked frozen local constructions required')
    manifest=json.loads((construction/'manifest.json').read_text())
    if witness['source_sha256'].get(str(construction/'manifest.json'))!=api.digest(construction/'manifest.json'):
        raise ValueError('Construction review belongs to another manifest')
    row=next(r for r in manifest['entries'] if r['key']=='occ_2_3__pose0__cap109')
    out.mkdir(parents=True,exist_ok=False);entries=[]
    from rdkit import Chem
    for label in ('A','AB'):
        item=next(r for r in row['molecules'] if r['id']==label)
        mol=Chem.MolFromMolFile(item['molecular_graph'],removeHs=False)
        if item['charge']!=0 or item['atoms'] not in (19,38) or mol is None or Chem.GetFormalCharge(mol)!=0:
            raise ValueError('These first pilots are neutral N-methylphthalimide monomer/pair controls')
        expected=Chem.MolToSmiles(Chem.MolFromSmiles('CN1C(=O)c2ccccc2C1=O'))
        if any(Chem.MolToSmiles(Chem.RemoveHs(part))!=expected for part in Chem.GetMolFrags(mol,asMols=True)):
            raise ValueError('Unmeasured motif cannot use the first local resource pilots')
        case=out/label;case.mkdir();xyz=case/'input.xyz';xyz.write_bytes(Path(item['geometry']).read_bytes())
        inputs={str(p):api.digest(p) for p in (construction/'manifest.json',review_path,Path(item['geometry']),Path(item['molecular_graph']),xyz,PROFILE)}
        inputs.update(witness['source_sha256'])
        sources=[Path(__file__).resolve(),ROOT/'src/redox/qm/dft.py',ROOT/'src/redox/core/common.py',ROOT/'src/redox/workflows/resources.py',
            ROOT/'scripts/pipeline/run_pi_stacking.py',ROOT/'scripts/fleet/tensor_runtime_worker.py',
            ROOT/'scripts/fleet/exclusive_python_worker.py',ROOT/'scripts/fleet/bounded_host_python_worker.py']
        plan=dict(task=dict(key=label,atoms=item['atoms'],charge=0,mult=1,source_parent=row['parent_key'],
            role='frozen_local_motif_electronic_resource_control'),settings=SETTINGS,
            input_xyz=str(xyz),geometry_sha256=api.digest(xyz),input_sha256=inputs,
            producer_sha256={str(p):api.digest(p) for p in sources},execution_driver=str(Path(__file__).resolve()),
            runtime_profile=str(PROFILE),raw_root=str(ROOT/'calcs/local_motif_interactions'/out.name/label),
            maximum_host_RSS_GiB=32,minimum_free_RAM_GiB=32,maximum_GPU_MiB=14336,
            minimum_GPU_MiB=15360,resource_class='prospective_19_or_38_atom_frozen_energy_only',
            primary_candidate_admission=False,thermochemistry_computed=False,
            experimental_smooth_surface_used=False,created=time.time())
        atomic_write_json(case/'plan.json',plan);entries.append(dict(key=label,campaign=str(case),plan_sha256=api.digest(case/'plan.json')))
    atomic_write_json(out/'campaign.json',dict(status='two_local_energy_resource_pilots_prepared',entries=entries,
        source_construction=str(construction),source_review=str(review_path),primary_candidate_admission=False,
        scientific_energies_transferred=False,large_polymer_or_Hessian_resource_qualification=False))


def run(campaign):
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner();plan=api.read_plan(campaign);api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    if plan['settings']!=SETTINGS or plan['primary_candidate_admission'] or plan['task']['atoms'] not in (19,38):raise ValueError('Wrong local energy scope')
    caps=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
    if len(os.sched_getaffinity(0))!=2 or set(os.environ.get(k) for k in caps)!={'2'}:raise ValueError('Two disjoint CPUs and thread caps=2 required')
    actual=device(int(os.environ['CUDA_VISIBLE_DEVICES']));raw=Path(plan['raw_root']);raw.mkdir(parents=True,exist_ok=False)
    stage=campaign/'electronic';stage.mkdir(exist_ok=False);status=stage/'status.json'
    report=dict(status='frozen_local_energy_resource_pilot_running',started=time.time(),raw_root=str(raw),
        owner=dict(host=socket.gethostname().split('.')[0],pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),
        actual_device=actual,primary_candidate_admission=False,thermochemistry_computed=False,results={})
    atomic_write_json(status,report)
    try:
        from ase.io import read
        from pyscf import gto
        from redox.qm import dft
        from gpu4pyscf import dft as gdft
        from gpu4pyscf.solvent import smd
        if smd.get_cds_legacy.__module__!='gpu4pyscf.solvent.smd':raise ValueError('Local control requires the unchanged legacy surface model')
        atoms=read(plan['input_xyz'])
        if len(atoms)!=plan['task']['atoms'] or set(atoms.get_chemical_symbols())-{'C','H','N','O'}:raise ValueError('Wrong local molecular resource class')
        monitor=api.load_module('local_point_memory_observation',ROOT/'scripts/pipeline/run_pi_stacking.py')
        with monitor.MemoryMonitor() as memory,(raw/'calculation.log').open('x') as log,redirect_stdout(log),redirect_stderr(log):
            for medium in SETTINGS['media']:
                gc.collect();dft._free_gpu_pool();folder=raw/medium;folder.mkdir()
                mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions.tolist())),unit='Angstrom',
                    basis=SETTINGS['basis'],charge=0,spin=0,verbose=0,max_memory=6000)
                if medium=='gas':
                    mf=gdft.RKS(mol).density_fit(only_dfj=False);mf.xc=SETTINGS['xc'];mf.nlc=SETTINGS['nlc']
                else:
                    mf=dft._build_smd_mf(mol,SETTINGS['xc'],medium,'gpu',nlc=SETTINGS['nlc']);mf.with_solvent.sasa_ng=590
                mf.conv_tol=SETTINGS['SCF_conv_tol'];mf.conv_tol_grad=SETTINGS['SCF_conv_tol_grad'];mf.max_cycle=200
                mf.grids.level=3;mf.nlcgrids.level=3;mf.chkfile=str(folder/'SCF.h5')
                mf.callback=lambda env:atomic_write_json(folder/'progress.json',dict(cycle=int(env['cycle'])+1,E_Ha=float(env['e_tot']),checked=time.time(),qualified=False))
                atomic_write_json(folder/'settings.json',dict(SETTINGS,medium=medium,geometry_sha256=plan['geometry_sha256']))
                mf.kernel();diag=dft._orbital_diag(mf,0)
                if not mf.converged or not math.isfinite(mf.e_tot) or abs(diag.get('spin_contam',math.inf))>.1:raise ValueError('Local frozen SCF/spin gate failed')
                point=dict(medium=medium,E_Ha=float(mf.e_tot),E_eV=float(mf.e_tot)*HARTREE_EV,
                    converged=True,orbital_diagnostics=diag,geometry_sha256=plan['geometry_sha256'],settings=SETTINGS,
                    charge=0,mult=1,geometry_optimized=False,thermochemistry_computed=False,independent_review_pending=True)
                atomic_write_json(folder/'energy.json',point);report['results'][medium]=point;atomic_write_json(status,report)
                del mf,mol;gc.collect();dft._free_gpu_pool()
        if not 0<memory.peak<=plan['maximum_GPU_MiB']:raise ValueError('Prospective local GPU ceiling exceeded')
        api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
        report.update(status='frozen_local_energy_finished_independent_review_pending',observed_peak_GPU_MiB=memory.peak,
            output_sha256={str(p):api.digest(p) for p in raw.rglob('*') if p.is_file() and p.suffix in ('.json','.h5')})
    except BaseException as error:report.update(status='frozen_local_energy_failed_preserved',error=repr(error));raise
    finally:report['finished']=time.time();atomic_write_json(status,report)


def launch(campaign,gpu,preflight_path):
    with (campaign/'dispatch.lock').open('a+') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        return _launch(campaign,gpu,preflight_path)


def _launch(campaign,gpu,preflight_path):
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation,free_slots,GATE
    plan=api.read_plan(campaign);api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    declared=json.loads((campaign.parent/'campaign.json').read_text())
    entries=[r for r in declared['entries'] if r['campaign']==str(campaign)]
    if len(entries)!=1 or entries[0]['plan_sha256']!=api.digest(campaign/'plan.json'):
        raise ValueError('Local pilot plan differs from its exact declared campaign')
    profile=json.loads(PROFILE.read_text());preflight=json.loads(preflight_path.read_text());host=socket.gethostname().split('.')[0]
    if (str(Path(sys.executable))!=profile['python'] or not preflight['ok'] or preflight['python']!=profile['python']
            or not preflight['GPU_imports_blocked_in_CPU_preflight'] or preflight['GPU_modules_loaded']):raise ValueError('CPU-safe exact-interpreter preflight required')
    if any((campaign/p).exists() for p in ['launch.json','electronic','residency_guard']):raise ValueError('Preserve an existing scientific attempt')
    placement=api.load_module('local_point_shared_host_promises',ROOT/'scripts/pipeline/coordinate_geometry_batch.py')
    with allocation_guard():
        capacity=host_capacity(2,minimum_ram_gib=64);promise=placement.remaining_promises()
        if not placement.can_fit(capacity['available_RAM_bytes'],promise,bound=32*2**30):return dict(status='resource_capacity_unavailable')
        cpus=cpu_allocation(2);actual=device(gpu)
        if (host,gpu) not in free_slots([host]):return dict(status='atomic_placement_unavailable')
        argv=['taskset','--cpu-list',','.join(map(str,cpus)),str(GATE),'run',str(gpu),'--label','redox-local-motif-electronic-resource-pilot','--','env',
            *[k+'=2' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')],
            'PYTHONPATH='+str(ROOT/'src'),'LD_LIBRARY_PATH='+profile['library_path'],profile['python'],str(ROOT/'scripts/fleet/exclusive_python_worker.py'),
            '--guard-out',str(campaign/'residency_guard'),str(ROOT/'scripts/fleet/bounded_host_python_worker.py'),'--out',str(campaign/'host_guard'),
            '--maximum-rss-gib','32','--minimum-free-gib','32',str(ROOT/'scripts/fleet/tensor_runtime_worker.py'),'--profile',str(PROFILE),
            '--out',str(campaign/'tensor_runtime'),str(Path(__file__).resolve()),'run','--campaign',str(campaign)]
        with (campaign/'worker.log').open('x') as log:process=subprocess.Popen(argv,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        time.sleep(.1);owner=dict(host=host,pid=process.pid,process_start=Path(f'/proc/{process.pid}/stat').read_text().rsplit(')',1)[1].split()[19])
        row=dict(status='gate_started_scientific_start_pending',owner=owner,actual_device=actual,CPU_affinity=cpus,argv=argv,
            plan_sha256=api.digest(campaign/'plan.json'),preflight_sha256=api.digest(preflight_path),started=time.time())
        atomic_write_json(campaign/'launch.json',row);return row


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','run','launch'))
    for name in ('campaign','construction','review','preflight'):p.add_argument('--'+name,type=Path)
    p.add_argument('--gpu',type=int);a=p.parse_args()
    if a.action=='prepare':
        from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
        with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
        prepare(a.construction.resolve(),a.review.resolve(),a.campaign.resolve())
    elif a.action=='run':run(a.campaign.resolve())
    else:print(json.dumps(launch(a.campaign.resolve(),a.gpu,a.preflight.resolve())))


if __name__=='__main__':main()
