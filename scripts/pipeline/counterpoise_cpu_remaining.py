#!/usr/bin/env python3
"""Checkpoint one remaining CPU counterpoise term after matched pilot gates.

The three permitted terms are AB, B-own and B-ghost. Every term uses the exact
saved A-own grid and full-complex RI basis. A completed term remains unqualified
until independent checkpoint and orbital-stationarity review.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import threading
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

CORE=ROOT/'scripts/validation/local_counterpoise.py'
FACTORY=ROOT/'scripts/validation/run_cpu_counterpoise_reference.py'
CONTROL=ROOT/'results/reviews/2026-10-26/counterpoise_cpu_pilot/remaining_water_check.json'
CAPS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
EXPECTED={'AB':(168,1026,0),'B_own':(84,513,0),'B_ghost':(84,1026,19)}


def inputs(case,geometry_review,own,own_stationarity,ghost,ghost_stationarity):
    if any(not p.is_relative_to(ROOT.resolve()) for p in (case,geometry_review,own,own_stationarity,ghost,ghost_stationarity)):
        raise ValueError('Inputs must be on shared project storage')
    helper=api.load_module('remaining_source_geometry_binding',ROOT/'scripts/pipeline/counterpoise_cpu_pilot.py')
    plan,sources=helper.source_hashes(case,geometry_review)
    control=json.loads(CONTROL.read_text());api.verify_sources(control['source_sha256'])
    if (control['status']!='remaining_water_counterpoise_terms_stationary' or not control['passed']
            or len(control['terms'])!=3 or {x['role'] for x in control['terms']}!=set(EXPECTED)
            or not all(x['passed'] for x in control['terms'])):
        raise ValueError('All three source-bound water controls are required')
    records=[]
    for raw,role,gradient_path,expected_status,expected_gradient in (
        (own,'A_own',own_stationarity,'CPU_pilot_converged_independent_review_pending','CPU_pilot_energy_and_orbital_stationarity_checked'),
        (ghost,'A_ghost',ghost_stationarity,'CPU_ghost_pilot_converged_independent_review_pending','CPU_ghost_energy_and_orbital_stationarity_checked')):
        status=json.loads((raw/'status.json').read_text())
        if status.get('status')!=expected_status:raise ValueError('Both matched A terms must complete before remaining terms')
        state=json.loads((raw/role/'energy.json').read_text());launch=json.loads((raw/'launch.json').read_text());gradient=json.loads(gradient_path.read_text())
        if (state!=status or not state['converged'] or state['qualified'] or state['source_sha256']!=launch['source_sha256']
                or gradient['status']!=expected_gradient or not gradient['passed'] or gradient['full_counterpoise_computed']
                or gradient['solvent_counterpoise_computed'] or abs(gradient['reference_energy_Ha']-state['E_Ha'])>1e-9):
            raise ValueError('Matched A-term energy or stationarity gate differs')
        api.verify_sources(state['source_sha256']);api.verify_sources(gradient['source_sha256'])
        bound=[raw/'launch.json',raw/'status.json',raw/role/'energy.json',raw/role/'SCF.h5',raw/'grid_coordinates_AU.npy',raw/'grid_weights.npy']
        if any(gradient['source_sha256'].get(str(p))!=api.digest(p) for p in bound):
            raise ValueError('Stationarity result does not bind the exact A-term checkpoint and grid')
        sources.update({str(p):api.digest(p) for p in [*bound,gradient_path]});records.append((state,launch))
    if (set(records[0][0]['common_grid_hashes'])!={'grid_coordinates_AU.npy','grid_weights.npy'}
            or records[0][0]['common_grid_hashes']!=records[1][0]['common_grid_hashes']
            or records[0][1]['versions']!=records[1][1]['versions']):
        raise ValueError('Matched A terms use different common grids or numerical environments')
    for name,digest in records[0][0]['common_grid_hashes'].items():
        if api.digest(own/name)!=digest or api.digest(ghost/name)!=digest:raise ValueError('Matched common grid changed')
    sources.update({str(p):api.digest(p) for p in (CORE,FACTORY,CONTROL,Path(__file__).resolve())})
    return plan,sources,records[0][0],records[0][1]['versions']


def run(role,case,geometry_review,own,own_stationarity,ghost,ghost_stationarity,raw,threads):
    import numpy as np
    import psutil
    import pyscf
    from pyscf import dft,df,lib
    launch=json.loads((raw/'launch.json').read_text());plan,sources,own_state,versions=inputs(case,geometry_review,own,own_stationarity,ghost,ghost_stationarity)
    if (role not in EXPECTED or launch['role']!=role or versions!={'pyscf':pyscf.__version__,'numpy':np.__version__}
            or launch['versions']!=versions or sources!=launch['source_sha256'] or sorted(os.sched_getaffinity(0))!=launch['cpu_affinity']
            or any(os.environ.get(k)!=str(threads) for k in CAPS)):
        raise ValueError('Term role, numerical environment, source or CPU allocation differs')
    core=api.load_module('remaining_term_counterpoise_protocol',CORE);fixed=api.load_module('remaining_term_fixed_auxiliary',FACTORY)
    if plan['settings']!=core.SETTINGS:raise ValueError('Counterpoise method differs')
    process=psutil.Process();limit=32*1024**3;peak=[process.memory_info().rss];stop=threading.Event()
    def watch():
        while not stop.wait(2):
            rss=process.memory_info().rss;peak[0]=max(peak[0],rss)
            if rss>limit:os.kill(os.getpid(),signal.SIGTERM);return
    def interrupted(signum,frame):raise MemoryError('Term stopped by sampled RSS guard or termination')
    signal.signal(signal.SIGTERM,interrupted);owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]);status=raw/'status.json'
    atomic_write_json(status,dict(status='CPU_remaining_term_initializing_unqualified',role=role,owner=owner,source_sha256=sources,started=time.time(),qualified=False))
    watcher=threading.Thread(target=watch,daemon=True);watcher.start()
    try:
        geometry=json.loads(Path(plan['input_geometry']).read_text());mols=core.molecules(geometry['symbols'],geometry['coordinates_A'],geometry['A_indices']);full=mols['AB'];mol=mols[role]
        electrons,basis,ghosts=EXPECTED[role]
        if (mol.nelectron!=electrons or mol.nao_nr()!=basis or int(sum(mol.atom_charges()==0))!=ghosts):
            raise ValueError('Remaining term nuclei, electrons or AO basis differ')
        grids=[]
        for name in ('grid_coordinates_AU.npy','grid_weights.npy'):
            target=raw/name;shutil.copyfile(own/name,target)
            if api.digest(target)!=own_state['common_grid_hashes'][name]:raise ValueError('Common grid differs from admitted A terms')
            grids.append(np.load(target))
        coords,weights=grids;aux=df.addons.make_auxmol(full,None);folder=raw/role;folder.mkdir(exist_ok=False)
        mf=dft.RKS(mol).density_fit(only_dfj=False);mf.with_df.auxmol=aux.copy();mf.xc=core.SETTINGS['xc'];mf.nlc=core.SETTINGS['nlc'];mf.small_rho_cutoff=0.
        mf.conv_tol=core.SETTINGS['SCF_conv_tol'];mf.conv_tol_grad=core.SETTINGS['SCF_conv_tol_grad'];mf.max_cycle=core.SETTINGS['SCF_max_cycle'];mf.grids.level=3;mf.nlcgrids.level=3
        for grid in (mf.grids,mf.nlcgrids):grid.coords=coords.copy();grid.weights=weights.copy()
        mf.chkfile=str(folder/'SCF.h5');lib.num_threads(threads)
        def callback(env):
            rss=process.memory_info().rss;peak[0]=max(peak[0],rss)
            atomic_write_json(folder/'progress.json',dict(cycle=int(env['cycle'])+1,E_Ha=float(env['e_tot']),host_RSS_bytes=rss,peak_RSS_bytes=peak[0],checked=time.time(),qualified=False))
            if rss>limit:raise MemoryError('Remaining term exceeds sampled 32 GiB RSS guard')
        mf.callback=callback;previous=launch.get('resume_from');initial=None
        if previous:
            checkpoint=Path(previous['checkpoint'])
            if api.digest(checkpoint)!=previous['checkpoint_sha256']:raise ValueError('Restart checkpoint changed')
            saved=lib.chkfile.load_mol(str(checkpoint))
            if (saved.nao_nr()!=mol.nao_nr() or saved.nelectron!=mol.nelectron or not np.array_equal(saved.atom_charges(),mol.atom_charges())
                    or not np.allclose(saved.atom_coords(),mol.atom_coords(),rtol=0,atol=1e-10)):
                raise ValueError('Restart checkpoint molecular identity differs')
            initial=mf.from_chk(str(checkpoint),project=False)
        atomic_write_json(status,dict(status='CPU_remaining_term_SCF_running_unqualified',role=role,owner=owner,started=time.time(),qualified=False))
        started=time.time()
        with fixed.fixed_complex_auxiliary_factory(df.addons,full,mols.values()):energy=float(mf.kernel(dm0=initial))
        actual=mf.with_df.auxmol
        if (not mf.converged or not np.isfinite(energy) or actual.nao_nr()!=aux.nao_nr()
                or not np.array_equal(actual.atom_charges(),aux.atom_charges())
                or not np.allclose(actual.atom_coords(),aux.atom_coords(),rtol=0,atol=1e-10)):
            raise ValueError('Term SCF failed or full-complex auxiliary basis changed')
        _,after,_,_=inputs(case,geometry_review,own,own_stationarity,ghost,ghost_stationarity)
        if after!=sources:raise ValueError('Term sources changed during calculation')
        result=dict(status='CPU_remaining_term_converged_independent_review_pending',role=role,E_Ha=energy,
            converged=True,electrons=electrons,basis_functions=basis,ghost_centers=ghosts,auxiliary_functions=int(actual.nao_nr()),
            grid_points=len(weights),common_grid_hashes=own_state['common_grid_hashes'],settings=core.SETTINGS,seconds=time.time()-started,
            peak_RSS_bytes=peak[0],owner=owner,checkpoint_sha256=api.digest(folder/'SCF.h5'),source_sha256=sources,
            primary_candidate_admission=False,full_counterpoise_computed=False,qualified=False)
        atomic_write_json(folder/'energy.json',result);atomic_write_json(status,result)
    except BaseException as error:
        atomic_write_json(status,dict(status='CPU_remaining_term_interrupted_or_failed_preserved',role=role,error=repr(error),owner=owner,
            peak_RSS_bytes=peak[0],source_sha256=sources,qualified=False,checkpoint_preserved=(raw/role/'SCF.h5').exists()));raise
    finally:
        stop.set();watcher.join(timeout=3)


def launch(role,case,geometry_review,own,own_stationarity,ghost,ghost_stationarity,raw,threads,
           resume_from=None,exclude_cpus=(),allowed_cpus=None,minimum_ram_gib=80):
    import numpy as np
    import pyscf
    if role not in EXPECTED or threads<1 or pyscf.__version__!='2.14.0':raise ValueError('Supported term and qualified numerical environment required')
    if not raw.is_relative_to((ROOT/'calcs').resolve()):raise ValueError('Use durable project checkpoint storage')
    if raw.exists():raise FileExistsError('Preserve prior term attempts')
    plan,sources,_,versions=inputs(case,geometry_review,own,own_stationarity,ghost,ghost_stationarity)
    if versions!={'pyscf':pyscf.__version__,'numpy':np.__version__}:raise ValueError('Term and admitted references use different numerical environments')
    previous=None
    if resume_from is not None:
        old=resume_from.resolve()
        if not old.is_relative_to((ROOT/'calcs').resolve()):raise ValueError('Restart must use durable project storage')
        state=json.loads((old/'status.json').read_text());prior=json.loads((old/'launch.json').read_text())
        if (state.get('status')!='CPU_remaining_term_interrupted_or_failed_preserved' or state.get('role')!=role
                or not state.get('checkpoint_preserved') or state.get('source_sha256')!=sources
                or prior.get('source_sha256')!=sources or prior.get('versions')!=versions):
            raise ValueError('Restart requires a stopped, source-matched checkpoint for this role')
        owner=state['owner']
        if owner['host']!=os.uname().nodename:raise ValueError('Cannot establish remote restart-owner status')
        proc=Path('/proc')/str(owner['pid'])/'stat'
        if proc.exists() and proc.read_text().rsplit(')',1)[1].split()[19]==owner['process_start']:raise ValueError('Previous term owner is still live')
        checkpoint=old/role/'SCF.h5';previous=dict(checkpoint=str(checkpoint),checkpoint_sha256=api.digest(checkpoint),prior_attempt=str(old),prior_status_sha256=api.digest(old/'status.json'))
    with allocation_guard():
        original=sorted(os.sched_getaffinity(0))
        try:
            if allowed_cpus is not None:os.sched_setaffinity(0,set(allowed_cpus))
            capacity=host_capacity(threads,minimum_ram_gib=minimum_ram_gib);choices=cpu_allocation(threads+len(set(exclude_cpus)))
            affinity=[cpu for cpu in choices if cpu not in set(exclude_cpus)][:threads]
            if len(affinity)!=threads:raise RuntimeError('Insufficient disjoint CPU allocation')
        finally:os.sched_setaffinity(0,set(original))
        raw.mkdir(parents=True);(raw/'scratch').mkdir();env=dict(os.environ,PYSCF_TMPDIR=str(raw/'scratch'))
        for key in CAPS:env[key]=str(threads)
        command=[sys.executable,str(Path(__file__).resolve()),'run','--role',role,'--case',str(case),'--geometry-review',str(geometry_review),
            '--own',str(own),'--own-stationarity',str(own_stationarity),'--ghost',str(ghost),'--ghost-stationarity',str(ghost_stationarity),'--raw',str(raw),'--threads',str(threads)]
        atomic_write_json(raw/'launch.json',dict(status='CPU_remaining_term_launch_prepared',role=role,case=str(case),threads=threads,cpu_affinity=affinity,
            capacity=capacity,source_sha256=sources,interpreter=sys.executable,command=command,versions=versions,resume_from=previous,
            maximum_host_RSS_GiB=32,full_counterpoise_computed=False,qualified=False))
        with (raw/'worker.log').open('x') as log:
            process=subprocess.Popen(command,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,env=env,cwd=ROOT,
                start_new_session=True,preexec_fn=lambda:os.sched_setaffinity(0,set(affinity)))
        atomic_write_json(raw/'launcher.json',dict(pid=process.pid,host=os.uname().nodename,launched=time.time(),cpu_affinity=affinity))
    return process.pid


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('launch','run'));p.add_argument('--role',choices=EXPECTED,required=True)
    for name in ('case','geometry-review','own','own-stationarity','ghost','ghost-stationarity','raw'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--threads',type=int,default=8);p.add_argument('--resume-from',type=Path);a=p.parse_args()
    values=[getattr(a,k).resolve() for k in ('case','geometry_review','own','own_stationarity','ghost','ghost_stationarity','raw')]
    if a.action=='launch':print(json.dumps(dict(pid=launch(a.role,*values,a.threads,a.resume_from),raw=str(values[-1]))))
    else:run(a.role,*values,a.threads)
