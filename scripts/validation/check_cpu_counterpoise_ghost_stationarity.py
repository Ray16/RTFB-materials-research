#!/usr/bin/env python3
"""Recompute A-ghost energy and orbital stationarity for a reviewed CPU pilot."""
import argparse
import json
import math
import os
from pathlib import Path
import signal
import threading
import time

import numpy as np
import psutil
import pyscf
from pyscf import dft,df,lib
from pyscf.scf.chkfile import load_scf

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

CORE=ROOT/'scripts/validation/local_counterpoise.py'
FACTORY=ROOT/'scripts/validation/run_cpu_counterpoise_reference.py'
CONTROL=ROOT/'results/reviews/2026-10-26/counterpoise_cpu_pilot/ghost_water_check.json'


def evaluate(geometry,checkpoint,coords_path,weights_path,reference_energy,expected_ghost_centers):
    core=api.load_module('ghost_stationarity_protocol',CORE);fixed=api.load_module('ghost_stationarity_fixed_auxiliary',FACTORY)
    mols=core.molecules(geometry['symbols'],geometry['coordinates_A'],geometry['A_indices']);mol=mols['A_ghost'];full=mols['AB'];saved,state=load_scf(str(checkpoint))
    identity=api.load_module('ghost_stationarity_identity',ROOT/'scripts/validation/check_cpu_counterpoise_ghost_control.py')
    identity.verify_ghost_checkpoint(saved,mol,full,expected_ghost_centers=expected_ghost_centers)
    coords=np.load(coords_path);weights=np.load(weights_path)
    if coords.shape!=(len(weights),3) or not np.isfinite(coords).all() or not np.isfinite(weights).all():raise ValueError('Invalid common grid')
    mf=dft.RKS(mol).density_fit(only_dfj=False);aux=df.addons.make_auxmol(full,None);mf.with_df.auxmol=aux.copy()
    mf.xc=core.SETTINGS['xc'];mf.nlc=core.SETTINGS['nlc'];mf.small_rho_cutoff=0.;mf.grids.level=3;mf.nlcgrids.level=3
    for grid in (mf.grids,mf.nlcgrids):grid.coords=coords.copy();grid.weights=weights.copy()
    dm=mf.make_rdm1(state['mo_coeff'],state['mo_occ']);hcore=mf.get_hcore();overlap=mf.get_ovlp()
    with fixed.fixed_complex_auxiliary_factory(df.addons,full,mols.values()):
        veff=mf.get_veff(mol,dm);fock=mf.get_fock(h1e=hcore,s1e=overlap,vhf=veff,dm=dm)
        energy=float(mf.energy_tot(dm=dm,h1e=hcore,vhf=veff));gradient=np.asarray(mf.get_grad(state['mo_coeff'],state['mo_occ'],fock))
    actual=mf.with_df.auxmol
    if (actual.nao_nr()!=aux.nao_nr() or not np.array_equal(actual.atom_charges(),aux.atom_charges())
            or not np.allclose(actual.atom_coords(),aux.atom_coords(),rtol=0,atol=1e-10)):
        raise ValueError('Full-complex auxiliary basis was not retained')
    if not math.isfinite(energy) or not np.isfinite(gradient).all():raise ValueError('Nonfinite ghost energy or gradient')
    norm=float(np.linalg.norm(gradient));error=abs(energy-reference_energy)
    return dict(energy_Ha=energy,reference_energy_Ha=reference_energy,energy_error_Ha=error,
        checkpoint_energy_error_Ha=abs(energy-float(state['e_tot'])),orbital_gradient_norm_au=norm,
        orbital_gradient_threshold_au=core.SETTINGS['SCF_conv_tol_grad'],electrons=int(mol.nelectron),
        basis_functions=int(mol.nao_nr()),ghost_centers=int(sum(mol.atom_charges()==0)),auxiliary_functions=int(actual.nao_nr()))


def check(raw,review_path,out):
    raw=raw.resolve();review_path=review_path.resolve();out=out.resolve()
    if out.exists():raise FileExistsError('Preserve previous ghost stationarity attempts')
    if not raw.is_relative_to((ROOT/'calcs').resolve()) or not out.is_relative_to((ROOT/'calcs').resolve()):
        raise ValueError('Use durable project storage')
    scratch=Path(os.environ.get('PYSCF_TMPDIR','/nonexistent')).resolve()
    if not scratch.is_dir() or not scratch.is_relative_to((ROOT/'calcs').resolve()):raise ValueError('Durable PySCF scratch required')
    review=json.loads(review_path.read_text())
    if (review['status']!='CPU_ghost_pilot_checkpoint_and_sampled_resources_reviewed_gradient_pending'
            or review['physical_gradient_recomputed'] or review['full_counterpoise_computed']):
        raise ValueError('Completed independent ghost checkpoint review required')
    api.verify_sources(review['source_sha256']);control=json.loads(CONTROL.read_text());api.verify_sources(control['source_sha256'])
    if control['status']!='water_ghost_energy_and_orbital_stationarity_checked' or not control['passed']:
        raise ValueError('Source-bound small ghost control required')
    launch=json.loads((raw/'launch.json').read_text());record=json.loads((raw/'A_ghost/energy.json').read_text())
    if (launch['versions']!={'pyscf':pyscf.__version__,'numpy':np.__version__}
            or record['status']!='CPU_ghost_pilot_converged_independent_review_pending'
            or abs(record['E_Ha']-review['energy_Ha'])>1e-9):
        raise ValueError('Ghost numerical environment or reviewed energy differs')
    api.verify_sources(record['source_sha256'])
    checker=api.load_module('ghost_stationarity_prior_review',ROOT/'scripts/validation/review_counterpoise_cpu_pilot.py')
    checker.verify_owner_stopped(record['owner'])
    plan=api.read_plan(Path(launch['case']));geometry=json.loads(Path(plan['input_geometry']).read_text())
    checkpoint=raw/'A_ghost/SCF.h5';coords_path=raw/'grid_coordinates_AU.npy';weights_path=raw/'grid_weights.npy'
    if (api.digest(checkpoint)!=record['checkpoint_sha256'] or api.digest(coords_path)!=record['common_grid_hashes']['grid_coordinates_AU.npy']
            or api.digest(weights_path)!=record['common_grid_hashes']['grid_weights.npy']):
        raise ValueError('Ghost checkpoint or common grid changed')
    bound=[raw/'launch.json',raw/'status.json',raw/'A_ghost/energy.json',checkpoint,coords_path,weights_path]
    if any(review['source_sha256'].get(str(p))!=api.digest(p) for p in bound):
        raise ValueError('Independent review does not bind this exact ghost checkpoint')
    core=api.load_module('ghost_stationarity_protocol_gate',CORE)
    if core.SETTINGS!=record['settings'] or plan['settings']!=core.SETTINGS:raise ValueError('Ghost method differs')
    paths=bound+[review_path,Path(launch['case'])/'plan.json',Path(plan['input_geometry']),
           ROOT/'scripts/validation/check_cpu_counterpoise_ghost_control.py',CORE,FACTORY,CONTROL,Path(__file__).resolve()]
    source={str(p):api.digest(p) for p in paths};out.mkdir(parents=True);started=time.time()
    process=psutil.Process();peak=[process.memory_info().rss];stop=threading.Event();limit=32*1024**3
    def monitor():
        while not stop.wait(2):
            rss=process.memory_info().rss;peak[0]=max(peak[0],rss)
            if rss>limit:os.kill(os.getpid(),signal.SIGTERM);return
    def interrupted(signum,frame):raise MemoryError('Ghost stationarity stopped by sampled RSS guard or termination')
    signal.signal(signal.SIGTERM,interrupted);watcher=threading.Thread(target=monitor,daemon=True);watcher.start()
    try:
        computed=evaluate(geometry,checkpoint,coords_path,weights_path,record['E_Ha'],expected_ghost_centers=19)
    except BaseException as error:
        atomic_write_json(out/'failure.json',dict(status='CPU_ghost_stationarity_failed_preserved',error=repr(error),
            sampled_peak_RSS_bytes=peak[0],source_sha256=source,qualified=False));raise
    finally:
        stop.set();watcher.join(timeout=3)
    api.verify_sources(source);api.verify_sources(record['source_sha256']);api.verify_sources(review['source_sha256'])
    passed=(computed['energy_error_Ha']<=1e-8 and computed['checkpoint_energy_error_Ha']<=1e-8
            and computed['orbital_gradient_norm_au']<=core.SETTINGS['SCF_conv_tol_grad'])
    result=dict(status='CPU_ghost_energy_and_orbital_stationarity_checked' if passed else 'CPU_ghost_stationarity_failed_preserved',
        passed=passed,**computed,sampled_peak_RSS_bytes=peak[0],seconds=time.time()-started,source_sha256=source,
        full_counterpoise_computed=False,solvent_counterpoise_computed=False,
        scope='A-ghost frozen electronic state; no geometry gradient or complete interaction')
    atomic_write_json(out/'summary.json',result)
    if not passed:raise ValueError('Ghost checkpoint fails the declared stationarity threshold')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('raw','review','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--threads',type=int,default=8);a=p.parse_args()
    if a.threads<1 or pyscf.__version__!='2.14.0':raise ValueError('Qualified PySCF and positive thread count required')
    with allocation_guard():host_capacity(a.threads,minimum_ram_gib=80);os.sched_setaffinity(0,cpu_allocation(a.threads))
    if any(os.environ.get(k)!=str(a.threads) for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')):
        raise ValueError('Declared CPU threads and capped native pools required')
    lib.num_threads(a.threads);result=check(a.raw,a.review,a.out);print(json.dumps(result))
