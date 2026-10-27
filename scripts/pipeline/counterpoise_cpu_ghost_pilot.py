#!/usr/bin/env python3
"""Checkpointed A-ghost CPU resource pilot, gated on a stationary A-own term.

This is one term of a five-term gas counterpoise calculation. It is not a
complete interaction, solvent correction or candidate-admission certificate.
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
GHOST_CONTROL=ROOT/'results/reviews/2026-10-26/counterpoise_cpu_pilot/ghost_water_check.json'
CAPS=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')


def inputs(case,geometry_review,own,stationarity):
    if any(not p.is_relative_to(ROOT.resolve()) for p in (case,geometry_review,own,stationarity)):
        raise ValueError('Inputs must be on shared project storage')
    helper=api.load_module('ghost_source_geometry_binding',ROOT/'scripts/pipeline/counterpoise_cpu_pilot.py')
    plan,sources=helper.source_hashes(case,geometry_review)
    status=json.loads((own/'status.json').read_text())
    if status.get('status')!='CPU_pilot_converged_independent_review_pending':
        raise ValueError('A completed A-own term is required before ghost pilot')
    state=json.loads((own/'A_own/energy.json').read_text());launch=json.loads((own/'launch.json').read_text())
    gradient=json.loads(stationarity.read_text());control=json.loads(GHOST_CONTROL.read_text())
    if (state!=status or state['status']!='CPU_pilot_converged_independent_review_pending' or not state['converged']
            or state['source_sha256']!=launch['source_sha256'] or state['qualified']
            or gradient['status']!='CPU_pilot_energy_and_orbital_stationarity_checked' or not gradient['passed']
            or gradient['full_counterpoise_computed'] or gradient['solvent_counterpoise_computed']
            or abs(gradient['reference_energy_Ha']-state['E_Ha'])>1e-9):
        raise ValueError('A matching, independently stationary A-own term is required')
    if (control['status']!='water_ghost_energy_and_orbital_stationarity_checked' or not control['passed']
            or control['role']!='A_ghost' or control['orbital_gradient_norm_au']>control['orbital_gradient_threshold_au']):
        raise ValueError('Source-bound water ghost control is required')
    api.verify_sources(state['source_sha256']);api.verify_sources(gradient['source_sha256']);api.verify_sources(control['source_sha256'])
    grids=[own/'grid_coordinates_AU.npy',own/'grid_weights.npy']
    bound=[own/'launch.json',own/'status.json',own/'A_own/energy.json',own/'A_own/SCF.h5',*grids]
    if any(gradient['source_sha256'].get(str(p))!=api.digest(p) for p in bound):
        raise ValueError('Stationarity result does not bind this exact A-own checkpoint and grid')
    if set(state['common_grid_hashes'])!={p.name for p in grids} or any(api.digest(p)!=state['common_grid_hashes'][p.name] for p in grids):
        raise ValueError('A-own common-grid files changed')
    paths=[own/'launch.json',own/'status.json',own/'A_own/energy.json',own/'A_own/SCF.h5',*grids,
           stationarity,GHOST_CONTROL,CORE,FACTORY,Path(__file__).resolve()]
    sources.update({str(p):api.digest(p) for p in paths})
    return plan,sources,state


def run(case,geometry_review,own,stationarity,raw,threads):
    import numpy as np
    import psutil
    import pyscf
    from pyscf import dft,df,lib
    launch=json.loads((raw/'launch.json').read_text());plan,sources,own_state=inputs(case,geometry_review,own,stationarity)
    if (pyscf.__version__!='2.14.0' or launch['versions']!={'pyscf':pyscf.__version__,'numpy':np.__version__}
            or sources!=launch['source_sha256'] or sorted(os.sched_getaffinity(0))!=launch['cpu_affinity']
            or any(os.environ.get(k)!=str(threads) for k in CAPS)):
        raise ValueError('Ghost pilot software, source or CPU allocation differs')
    core=api.load_module('cpu_ghost_counterpoise_protocol',CORE);fixed=api.load_module('cpu_ghost_fixed_auxiliary',FACTORY)
    if plan['settings']!=core.SETTINGS:raise ValueError('Counterpoise protocol differs')
    process=psutil.Process();limit=32*1024**3;peak=[process.memory_info().rss];stop=threading.Event()
    def watch():
        while not stop.wait(2):
            current=process.memory_info().rss;peak[0]=max(peak[0],current)
            if current>limit:os.kill(os.getpid(),signal.SIGTERM);return
    def interrupted(signum,frame):raise MemoryError('Ghost pilot stopped by sampled RSS guard or termination')
    signal.signal(signal.SIGTERM,interrupted)
    owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    status=raw/'status.json';atomic_write_json(status,dict(status='CPU_ghost_pilot_initializing_unqualified',owner=owner,
        source_sha256=sources,started=time.time(),qualified=False))
    watcher=threading.Thread(target=watch,daemon=True);watcher.start()
    try:
        geometry=json.loads(Path(plan['input_geometry']).read_text());mols=core.molecules(geometry['symbols'],geometry['coordinates_A'],geometry['A_indices']);full=mols['AB'];mol=mols['A_ghost']
        if mol.nao_nr()!=1026 or full.nao_nr()!=1026 or mol.nelectron!=84 or int(sum(mol.atom_charges()==0))!=19:
            raise ValueError('Ghost nuclei, electron count or AO basis differs')
        grids=[]
        for name in ('grid_coordinates_AU.npy','grid_weights.npy'):
            source=own/name;target=raw/name;shutil.copyfile(source,target)
            if api.digest(target)!=own_state['common_grid_hashes'][name]:raise ValueError('Common grid differs from A-own term')
            grids.append(np.load(target))
        coords,weights=grids;aux=df.addons.make_auxmol(full,None)
        role=raw/'A_ghost';role.mkdir(exist_ok=False);mf=dft.RKS(mol).density_fit(only_dfj=False);mf.with_df.auxmol=aux.copy()
        mf.xc=core.SETTINGS['xc'];mf.nlc=core.SETTINGS['nlc'];mf.small_rho_cutoff=0.
        mf.conv_tol=core.SETTINGS['SCF_conv_tol'];mf.conv_tol_grad=core.SETTINGS['SCF_conv_tol_grad'];mf.max_cycle=core.SETTINGS['SCF_max_cycle']
        mf.grids.level=3;mf.nlcgrids.level=3
        for grid in (mf.grids,mf.nlcgrids):grid.coords=coords.copy();grid.weights=weights.copy()
        mf.chkfile=str(role/'SCF.h5');lib.num_threads(threads)
        def callback(env):
            current=process.memory_info().rss;peak[0]=max(peak[0],current)
            atomic_write_json(role/'progress.json',dict(cycle=int(env['cycle'])+1,E_Ha=float(env['e_tot']),
                host_RSS_bytes=current,peak_RSS_bytes=peak[0],checked=time.time(),qualified=False))
            if current>limit:raise MemoryError('Ghost pilot exceeds sampled 32 GiB RSS guard')
        mf.callback=callback
        previous=launch.get('resume_from');initial=None
        if previous:
            checkpoint=Path(previous['checkpoint'])
            if api.digest(checkpoint)!=previous['checkpoint_sha256']:raise ValueError('Ghost restart checkpoint changed')
            saved=lib.chkfile.load_mol(str(checkpoint))
            if (saved.nao_nr()!=mol.nao_nr() or saved.nelectron!=mol.nelectron
                    or not np.array_equal(saved.atom_charges(),mol.atom_charges())
                    or not np.allclose(saved.atom_coords(),mol.atom_coords(),rtol=0,atol=1e-10)):
                raise ValueError('Ghost restart molecular identity differs')
            initial=mf.from_chk(str(checkpoint),project=False)
        atomic_write_json(status,dict(status='CPU_ghost_pilot_SCF_running_unqualified',owner=owner,
            role='A_ghost',started=time.time(),qualified=False))
        started=time.time()
        with fixed.fixed_complex_auxiliary_factory(df.addons,full,mols.values()):energy=float(mf.kernel(dm0=initial))
        actual=mf.with_df.auxmol
        if (not mf.converged or not np.isfinite(energy) or actual.nao_nr()!=aux.nao_nr()
                or not np.array_equal(actual.atom_charges(),aux.atom_charges())
                or not np.allclose(actual.atom_coords(),aux.atom_coords(),rtol=0,atol=1e-10)):
            raise ValueError('Ghost SCF failed or full-complex auxiliary basis changed')
        _,after,_=inputs(case,geometry_review,own,stationarity)
        if after!=sources:raise ValueError('Ghost pilot sources changed during computation')
        result=dict(status='CPU_ghost_pilot_converged_independent_review_pending',role='A_ghost',E_Ha=energy,
            converged=True,electrons=int(mol.nelectron),basis_functions=int(mol.nao_nr()),
            ghost_centers=int(sum(mol.atom_charges()==0)),auxiliary_functions=int(actual.nao_nr()),grid_points=len(weights),
            common_grid_hashes=own_state['common_grid_hashes'],settings=core.SETTINGS,seconds=time.time()-started,
            peak_RSS_bytes=peak[0],owner=owner,checkpoint_sha256=api.digest(role/'SCF.h5'),source_sha256=sources,
            primary_candidate_admission=False,full_counterpoise_computed=False,qualified=False)
        atomic_write_json(role/'energy.json',result);atomic_write_json(status,result)
    except BaseException as error:
        atomic_write_json(status,dict(status='CPU_ghost_pilot_interrupted_or_failed_preserved',error=repr(error),owner=owner,
            peak_RSS_bytes=peak[0],source_sha256=sources,qualified=False,checkpoint_preserved=(raw/'A_ghost/SCF.h5').exists()))
        raise
    finally:
        stop.set();watcher.join(timeout=3)


def launch(case,geometry_review,own,stationarity,raw,threads,resume_from=None,exclude_cpus=(),allowed_cpus=None):
    import numpy as np
    import pyscf
    if threads<1 or pyscf.__version__!='2.14.0':raise ValueError('Qualified PySCF and positive CPU count required')
    if not raw.is_relative_to((ROOT/'calcs').resolve()):raise ValueError('Use durable project checkpoint storage')
    if raw.exists():raise FileExistsError('Preserve prior ghost-pilot attempts')
    plan,sources,own_state=inputs(case,geometry_review,own,stationarity)
    versions={'pyscf':pyscf.__version__,'numpy':np.__version__};previous=None
    if json.loads((own/'launch.json').read_text())['versions']!=versions:
        raise ValueError('Ghost and A-own numerical environments differ')
    if resume_from is not None:
        old=resume_from.resolve()
        if not old.is_relative_to((ROOT/'calcs').resolve()):raise ValueError('Restart must use durable project storage')
        old_status=json.loads((old/'status.json').read_text());old_launch=json.loads((old/'launch.json').read_text())
        if (old_status.get('status')!='CPU_ghost_pilot_interrupted_or_failed_preserved' or not old_status.get('checkpoint_preserved')
                or old_status.get('source_sha256')!=sources or old_launch.get('source_sha256')!=sources
                or old_launch.get('versions')!=versions):
            raise ValueError('Restart requires a stopped, source-matched ghost checkpoint')
        owner=old_status['owner']
        if owner['host']!=os.uname().nodename:raise ValueError('Cannot establish remote restart-owner status')
        proc=Path('/proc')/str(owner['pid'])/'stat'
        if proc.exists() and proc.read_text().rsplit(')',1)[1].split()[19]==owner['process_start']:
            raise ValueError('Previous ghost-pilot owner is still live')
        checkpoint=old/'A_ghost/SCF.h5'
        previous=dict(checkpoint=str(checkpoint),checkpoint_sha256=api.digest(checkpoint),prior_attempt=str(old),
                      prior_status_sha256=api.digest(old/'status.json'))
    with allocation_guard():
        original=sorted(os.sched_getaffinity(0))
        try:
            if allowed_cpus is not None:os.sched_setaffinity(0,set(allowed_cpus))
            capacity=host_capacity(threads,minimum_ram_gib=80);choices=cpu_allocation(threads+len(set(exclude_cpus)))
            affinity=[cpu for cpu in choices if cpu not in set(exclude_cpus)][:threads]
            if len(affinity)!=threads:raise RuntimeError('Insufficient disjoint CPU allocation')
        finally:
            os.sched_setaffinity(0,set(original))
        raw.mkdir(parents=True);(raw/'scratch').mkdir();env=dict(os.environ,PYSCF_TMPDIR=str(raw/'scratch'))
        for name in CAPS:env[name]=str(threads)
        command=[sys.executable,str(Path(__file__).resolve()),'run','--case',str(case),'--geometry-review',str(geometry_review),
                 '--own',str(own),'--stationarity',str(stationarity),'--raw',str(raw),'--threads',str(threads)]
        atomic_write_json(raw/'launch.json',dict(status='CPU_ghost_pilot_launch_prepared',case=str(case),role='A_ghost',threads=threads,
            cpu_affinity=affinity,capacity=capacity,source_sha256=sources,interpreter=sys.executable,command=command,
            versions=versions,resume_from=previous,maximum_host_RSS_GiB=32,full_counterpoise_computed=False,qualified=False))
        with (raw/'worker.log').open('x') as log:
            process=subprocess.Popen(command,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,env=env,cwd=ROOT,
                start_new_session=True,preexec_fn=lambda:os.sched_setaffinity(0,set(affinity)))
        atomic_write_json(raw/'launcher.json',dict(pid=process.pid,host=os.uname().nodename,launched=time.time(),cpu_affinity=affinity))
    return process.pid


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('launch','run'))
    for name in ('case','geometry-review','own','stationarity','raw'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--threads',type=int,default=8);p.add_argument('--resume-from',type=Path);a=p.parse_args()
    case=a.case.resolve();review=a.geometry_review.resolve();own=a.own.resolve();stationarity=a.stationarity.resolve();raw=a.raw.resolve()
    if a.action=='launch':print(json.dumps(dict(pid=launch(case,review,own,stationarity,raw,a.threads,a.resume_from),raw=str(raw))))
    else:run(case,review,own,stationarity,raw,a.threads)
