#!/usr/bin/env python3
"""CPU-only, checkpointed A-own counterpoise resource pilot.

No full counterpoise value follows from this single term. A future five-term
run requires independent resource and physical review of the completed pilot.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import threading
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity

CORE = ROOT/'scripts/validation/local_counterpoise.py'
CPU_REFERENCE = ROOT/'scripts/validation/run_cpu_counterpoise_reference.py'
CAPS = ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')


def source_hashes(case, review):
    plan=api.read_plan(case);api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    checked=json.loads(review.read_text());api.verify_sources(checked['source_sha256'])
    if (checked['status']!='prepared_motif_counterpoise_geometry_independently_checked_no_scientific_start'
            or Path(checked['campaign'])!=case or checked['parent_key']!=plan['source_construction_key']
            or checked['atoms']!=38 or checked['basis_functions']['A_own']!=513
            or checked['basis_functions']['AB']!=1026):
        raise ValueError('Matching independently reviewed geometry required')
    paths=(case/'plan.json',Path(plan['input_geometry']),review,CORE,CPU_REFERENCE,Path(__file__).resolve())
    return plan,{str(p):api.digest(p) for p in paths}


def run(case, review, raw, threads):
    import numpy as np
    import psutil
    import pyscf
    from pyscf import dft,df,lib
    from pyscf.dft import gen_grid
    if pyscf.__version__!='2.14.0' or len(os.sched_getaffinity(0))!=threads or any(os.environ.get(k)!=str(threads) for k in CAPS):
        raise ValueError('Validated PySCF and exclusive declared CPU/thread settings required')
    launch=json.loads((raw/'launch.json').read_text());plan,sources=source_hashes(case,review)
    if (sources!=launch['source_sha256'] or sorted(os.sched_getaffinity(0))!=launch['cpu_affinity']
            or launch['versions']!={'pyscf':pyscf.__version__,'numpy':np.__version__}):
        raise ValueError('Pilot identity or allocation differs from launch')
    core=api.load_module('cpu_counterpoise_original_protocol',CORE)
    if core.SETTINGS!=plan['settings']:
        raise ValueError('Counterpoise protocol differs from reviewed control')
    process=psutil.Process();limit=32*1024**3;peak=[process.memory_info().rss];stop=threading.Event()
    def watch():
        while not stop.wait(2):
            current=process.memory_info().rss;peak[0]=max(peak[0],current)
            if current>limit:
                os.kill(os.getpid(),signal.SIGTERM);return
    def interrupted(signum,frame):
        raise MemoryError('CPU pilot stopped by memory ceiling or termination signal')
    signal.signal(signal.SIGTERM,interrupted)
    owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    status=raw/'status.json';atomic_write_json(status,dict(status='CPU_pilot_initializing_unqualified',owner=owner,source_sha256=sources,started=time.time(),qualified=False))
    watcher=threading.Thread(target=watch,daemon=True);watcher.start()
    try:
        geometry=json.loads(Path(plan['input_geometry']).read_text());mols=core.molecules(geometry['symbols'],geometry['coordinates_A'],geometry['A_indices']);full=mols['AB'];mol=mols['A_own']
        grid=gen_grid.Grids(full);grid.level=3;grid.build(with_non0tab=False)
        auxiliary=df.addons.make_auxmol(full,None)
        if not (mol.nao_nr()==513 and full.nao_nr()==1026 and mol.nelectron==84):raise ValueError('Pilot geometry/electron/basis identity differs')
        np.save(raw/'grid_coordinates_AU.npy',grid.coords);np.save(raw/'grid_weights.npy',grid.weights)
        grid_hashes={name:api.digest(raw/name) for name in ('grid_coordinates_AU.npy','grid_weights.npy')}
        previous=launch.get('resume_from')
        if previous and previous['grid_hashes']!=grid_hashes:
            raise ValueError('Restart quadrature differs from the saved checkpoint')
        role=raw/'A_own';role.mkdir(exist_ok=False);mf=dft.RKS(mol).density_fit(only_dfj=False);mf.with_df.auxmol=auxiliary.copy()
        mf.xc=core.SETTINGS['xc'];mf.nlc=core.SETTINGS['nlc'];mf.small_rho_cutoff=0.
        mf.conv_tol=core.SETTINGS['SCF_conv_tol'];mf.conv_tol_grad=core.SETTINGS['SCF_conv_tol_grad'];mf.max_cycle=core.SETTINGS['SCF_max_cycle']
        mf.grids.level=3;mf.nlcgrids.level=3
        for grids in (mf.grids,mf.nlcgrids):grids.coords=grid.coords.copy();grids.weights=grid.weights.copy()
        mf.chkfile=str(role/'SCF.h5');lib.num_threads(threads)
        def callback(env):
            current=process.memory_info().rss;peak[0]=max(peak[0],current)
            atomic_write_json(role/'progress.json',dict(cycle=int(env['cycle'])+1,E_Ha=float(env['e_tot']),
                host_RSS_bytes=current,peak_RSS_bytes=peak[0],checked=time.time(),qualified=False))
            if current>limit:raise MemoryError('CPU pilot exceeds 32 GiB host RSS ceiling')
        mf.callback=callback
        atomic_write_json(status,dict(status='CPU_pilot_SCF_running_unqualified',owner=owner,role='A_own',started=time.time(),qualified=False))
        initial=None
        if previous:
            checkpoint=Path(previous['checkpoint'])
            if api.digest(checkpoint)!=previous['checkpoint_sha256']:
                raise ValueError('Restart checkpoint changed')
            saved=lib.chkfile.load_mol(str(checkpoint))
            if (saved.nao_nr()!=mol.nao_nr() or saved.nelectron!=mol.nelectron
                    or not np.array_equal(saved.atom_charges(),mol.atom_charges())
                    or not np.allclose(saved.atom_coords(),mol.atom_coords(),rtol=0,atol=1e-10)):
                raise ValueError('Restart checkpoint molecular identity differs')
            initial=mf.from_chk(str(checkpoint),project=False)
        fixed=api.load_module('cpu_counterpoise_fixed_auxiliary',CPU_REFERENCE)
        started=time.time()
        with fixed.fixed_complex_auxiliary_factory(df.addons,full,mols.values()):
            energy=float(mf.kernel(dm0=initial))
        if not mf.converged or not np.isfinite(energy):raise ValueError('CPU pilot SCF did not converge')
        actual_aux=mf.with_df.auxmol
        if (actual_aux.nao_nr()!=auxiliary.nao_nr()
                or not np.array_equal(actual_aux.atom_charges(),auxiliary.atom_charges())
                or not np.allclose(actual_aux.atom_coords(),auxiliary.atom_coords(),rtol=0,atol=1e-10)):
            raise ValueError('CPU SCF did not retain the full-complex auxiliary basis')
        _,current_sources=source_hashes(case,review)
        if current_sources!=sources:raise ValueError('Pilot sources changed during calculation')
        result=dict(status='CPU_pilot_converged_independent_review_pending',role='A_own',E_Ha=energy,
            converged=True,electrons=int(mol.nelectron),basis_functions=int(mol.nao_nr()),
            auxiliary_functions=int(mf.with_df.auxmol.nao_nr()),grid_points=len(grid.weights),common_grid_hashes=grid_hashes,
            settings=core.SETTINGS,seconds=time.time()-started,peak_RSS_bytes=peak[0],owner=owner,
            checkpoint_sha256=api.digest(role/'SCF.h5'),source_sha256=sources,
            primary_candidate_admission=False,full_counterpoise_computed=False,qualified=False)
        atomic_write_json(role/'energy.json',result);atomic_write_json(status,result)
    except BaseException as error:
        atomic_write_json(status,dict(status='CPU_pilot_interrupted_or_failed_preserved',error=repr(error),owner=owner,
            peak_RSS_bytes=peak[0],source_sha256=sources,qualified=False,checkpoint_preserved=(raw/'A_own/SCF.h5').exists()))
        raise
    finally:
        stop.set();watcher.join(timeout=3)


def launch(case, review, raw, threads, resume_from=None):
    import numpy as np
    import psutil
    import pyscf
    if threads < 1 or pyscf.__version__ != '2.14.0':
        raise ValueError('Positive thread count and qualified CPU PySCF version required')
    if not raw.is_relative_to((ROOT/'calcs').resolve()):
        raise ValueError('Pilot checkpoints must use durable project storage')
    if raw.exists():raise FileExistsError('Preserve existing pilot attempt and checkpoints')
    plan,sources=source_hashes(case,review)
    versions={'pyscf':pyscf.__version__,'numpy':np.__version__};previous=None
    if resume_from is not None:
        old=resume_from.resolve()
        if not old.is_relative_to((ROOT/'calcs').resolve()):raise ValueError('Restart must use a durable project checkpoint')
        old_status=json.loads((old/'status.json').read_text());old_launch=json.loads((old/'launch.json').read_text())
        if (old_status.get('status')!='CPU_pilot_interrupted_or_failed_preserved'
                or not old_status.get('checkpoint_preserved') or old_status.get('source_sha256')!=sources
                or old_launch.get('source_sha256')!=sources or old_launch.get('versions')!=versions):
            raise ValueError('Restart requires a stopped, source-matched checkpoint')
        owner=old_status['owner']
        if owner['host']!=os.uname().nodename:raise ValueError('Cannot establish remote restart-owner status')
        proc=Path('/proc')/str(owner['pid'])/'stat'
        if proc.exists() and proc.read_text().rsplit(')',1)[1].split()[19]==owner['process_start']:
            raise ValueError('Previous pilot owner is still live')
        checkpoint=old/'A_own/SCF.h5'
        grid_hashes={name:api.digest(old/name) for name in ('grid_coordinates_AU.npy','grid_weights.npy')}
        previous=dict(checkpoint=str(checkpoint),checkpoint_sha256=api.digest(checkpoint),grid_hashes=grid_hashes,
                      prior_attempt=str(old),prior_status_sha256=api.digest(old/'status.json'))
    with allocation_guard():
        capacity=host_capacity(threads,minimum_ram_gib=80);affinity=cpu_allocation(threads)
        raw.mkdir(parents=True);(raw/'scratch').mkdir()
        env=dict(os.environ,PYSCF_TMPDIR=str(raw/'scratch'))
        for name in CAPS:env[name]=str(threads)
        command=[sys.executable,str(Path(__file__).resolve()),'run','--case',str(case),'--review',str(review),'--raw',str(raw),'--threads',str(threads)]
        manifest=dict(status='CPU_pilot_launch_prepared',case=str(case),role='A_own',threads=threads,cpu_affinity=affinity,
            capacity=capacity,source_sha256=sources,interpreter=sys.executable,command=command,
            versions=versions,resume_from=previous,maximum_host_RSS_GiB=32,full_counterpoise_computed=False,qualified=False)
        atomic_write_json(raw/'launch.json',manifest)
        with (raw/'worker.log').open('x') as log:
            process=subprocess.Popen(command,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,env=env,
                cwd=ROOT,start_new_session=True,preexec_fn=lambda:os.sched_setaffinity(0,set(affinity)))
        atomic_write_json(raw/'launcher.json',dict(pid=process.pid,host=os.uname().nodename,launched=time.time(),cpu_affinity=affinity))
    return process.pid


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('launch','run'))
    for key in ('case','review','raw'):p.add_argument('--'+key,type=Path,required=True)
    p.add_argument('--threads',type=int,default=8);p.add_argument('--resume-from',type=Path);a=p.parse_args()
    case=a.case.resolve();review=a.review.resolve();raw=a.raw.resolve()
    if a.action=='launch':print(json.dumps(dict(pid=launch(case,review,raw,a.threads,a.resume_from),raw=str(raw))))
    else:run(case,review,raw,a.threads)
