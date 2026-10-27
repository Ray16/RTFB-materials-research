#!/usr/bin/env python3
"""Recompute energy and orbital gradient from a saved CPU counterpoise state.

This evaluates the stored density at the original common grid and full-complex
auxiliary basis; it is not a geometry gradient or a full counterpoise result.
"""
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


def evaluate(symbols,coordinates,indices,grid_coordinates,grid_weights,checkpoint,reference_energy):
    core=api.load_module('stationarity_counterpoise_protocol',CORE)
    fixed=api.load_module('stationarity_fixed_auxiliary',FACTORY)
    mols=core.molecules(symbols,coordinates,indices);mol=mols['A_own'];full=mols['AB']
    saved,state=load_scf(str(checkpoint))
    if (saved.nao_nr()!=mol.nao_nr() or saved.nelectron!=mol.nelectron
            or not np.array_equal(saved.atom_charges(),mol.atom_charges())
            or not np.allclose(saved.atom_coords(),mol.atom_coords(),rtol=0,atol=1e-10)):
        raise ValueError('Checkpoint molecule differs from declared A-own fragment')
    coords=np.load(grid_coordinates);weights=np.load(grid_weights)
    if coords.shape!=(len(weights),3) or not np.isfinite(coords).all() or not np.isfinite(weights).all():
        raise ValueError('Invalid common-grid arrays')
    mf=dft.RKS(mol).density_fit(only_dfj=False);aux=df.addons.make_auxmol(full,None);mf.with_df.auxmol=aux.copy()
    mf.xc=core.SETTINGS['xc'];mf.nlc=core.SETTINGS['nlc'];mf.small_rho_cutoff=0.
    mf.grids.level=3;mf.nlcgrids.level=3
    for grids in (mf.grids,mf.nlcgrids):grids.coords=coords.copy();grids.weights=weights.copy()
    dm=mf.make_rdm1(state['mo_coeff'],state['mo_occ']);hcore=mf.get_hcore();overlap=mf.get_ovlp()
    with fixed.fixed_complex_auxiliary_factory(df.addons,full,mols.values()):
        veff=mf.get_veff(mol,dm)
        fock=mf.get_fock(h1e=hcore,s1e=overlap,vhf=veff,dm=dm)
        energy=float(mf.energy_tot(dm=dm,h1e=hcore,vhf=veff))
        gradient=np.asarray(mf.get_grad(state['mo_coeff'],state['mo_occ'],fock))
    actual=mf.with_df.auxmol
    if (actual.nao_nr()!=aux.nao_nr() or not np.array_equal(actual.atom_charges(),aux.atom_charges())
            or not np.allclose(actual.atom_coords(),aux.atom_coords(),rtol=0,atol=1e-10)):
        raise ValueError('Full-complex auxiliary basis was not retained')
    if not math.isfinite(energy) or not np.isfinite(gradient).all():raise ValueError('Nonfinite state evaluation')
    return dict(recomputed_energy_Ha=energy,checkpoint_energy_Ha=float(state['e_tot']),reference_energy_Ha=float(reference_energy),
        reference_energy_error_Ha=abs(energy-reference_energy),checkpoint_energy_error_Ha=abs(energy-float(state['e_tot'])),
        orbital_gradient_norm_au=float(np.linalg.norm(gradient)),orbital_gradient_max_au=float(np.max(np.abs(gradient))),
        orbital_gradient_threshold_au=core.SETTINGS['SCF_conv_tol_grad'],grid_points=len(weights),
        auxiliary_functions=int(actual.nao_nr()),electrons=int(mol.nelectron),basis_functions=int(mol.nao_nr()))


def check_water(out):
    if out.exists():raise FileExistsError('Preserve earlier stationarity checks')
    review_path=ROOT/'results/reviews/2026-10-26/counterpoise_control/independent_review.json'
    review=json.loads(review_path.read_text());api.verify_sources(review['source_sha256'])
    if review['status']!='gas_counterpoise_CPU_GPU_control_independently_checked':raise ValueError('Reviewed control required')
    base=ROOT/'calcs/validation/2026-10-26_counterpoise/CPU_fixed_auxiliary_factory'
    summary_path=base/'summary.json';energy_path=base/'A_own/energy.json';checkpoint=base/'A_own/SCF.h5'
    paths=[summary_path,energy_path,checkpoint,base/'grid_coordinates_AU.npy',base/'grid_weights.npy',CORE,FACTORY,review_path,Path(__file__).resolve()]
    for p in paths[:5]:
        if review['source_sha256'].get(str(p))!=api.digest(p):raise ValueError('Reviewed source changed')
    before={str(p):api.digest(p) for p in paths};summary=json.loads(summary_path.read_text());record=json.loads(energy_path.read_text())
    out.mkdir(parents=True);start=time.time()
    result=evaluate(summary['source_symbols'],summary['source_coordinates_A'],summary['A_indices'],
                    paths[3],paths[4],checkpoint,record['E_Ha'])
    api.verify_sources(before)
    passed=(result['reference_energy_error_Ha']<=1e-8 and result['checkpoint_energy_error_Ha']<=1e-8
            and result['orbital_gradient_norm_au']<=result['orbital_gradient_threshold_au'])
    result.update(status='water_counterpoise_stationarity_recomputed' if passed else 'water_counterpoise_stationarity_failed_preserved',
        passed=passed,seconds=time.time()-start,source_sha256=before,
        scope='CPU water-dimer A-own electronic state only; no large-pilot or five-term qualification')
    atomic_write_json(out/'summary.json',result)
    if not passed:raise ValueError('Saved CPU control is not stationary at the declared threshold')
    return result


def check_pilot(raw,review_path,out):
    raw=raw.resolve();review_path=review_path.resolve();out=out.resolve()
    if out.exists():raise FileExistsError('Preserve earlier stationarity checks')
    if not raw.is_relative_to((ROOT/'calcs').resolve()) or not out.is_relative_to((ROOT/'calcs').resolve()):
        raise ValueError('Pilot inputs and output must be on durable project storage')
    scratch=Path(os.environ.get('PYSCF_TMPDIR','/nonexistent')).resolve()
    if not scratch.is_dir() or not scratch.is_relative_to((ROOT/'calcs').resolve()):
        raise ValueError('A durable PySCF scratch directory is required')
    review=json.loads(review_path.read_text())
    if (review['status']!='CPU_pilot_checkpoint_and_sampled_resources_reviewed_gradient_pending'
            or review['physical_gradient_recomputed'] or review['full_counterpoise_computed']):
        raise ValueError('A matching independent pilot checkpoint review is required')
    api.verify_sources(review['source_sha256'])
    record=json.loads((raw/'A_own/energy.json').read_text());launch=json.loads((raw/'launch.json').read_text())
    if launch['versions']!={'pyscf':pyscf.__version__,'numpy':np.__version__}:
        raise ValueError('Pilot and stationarity software versions differ')
    if record['status']!='CPU_pilot_converged_independent_review_pending' or abs(record['E_Ha']-review['energy_Ha'])>1e-9:
        raise ValueError('Pilot result differs from the independent review')
    reviewer=api.load_module('stationarity_pilot_owner_review',ROOT/'scripts/validation/review_counterpoise_cpu_pilot.py')
    reviewer.verify_owner_stopped(record['owner'])
    plan=api.read_plan(Path(launch['case']));geometry=json.loads(Path(plan['input_geometry']).read_text())
    checkpoint=raw/'A_own/SCF.h5';grid_coordinates=raw/'grid_coordinates_AU.npy';grid_weights=raw/'grid_weights.npy'
    if (api.digest(checkpoint)!=record['checkpoint_sha256'] or api.digest(grid_coordinates)!=record['common_grid_hashes']['grid_coordinates_AU.npy']
            or api.digest(grid_weights)!=record['common_grid_hashes']['grid_weights.npy']):
        raise ValueError('Pilot checkpoint or integration grid changed')
    paths=[raw/'launch.json',raw/'status.json',raw/'A_own/energy.json',checkpoint,grid_coordinates,grid_weights,
           review_path,CORE,FACTORY,Path(__file__).resolve()]
    before={str(p):api.digest(p) for p in paths};out.mkdir(parents=True);start=time.time()
    process=psutil.Process();peak=[process.memory_info().rss];stop=threading.Event();limit=32*1024**3
    def monitor():
        while not stop.wait(2):
            rss=process.memory_info().rss;peak[0]=max(peak[0],rss)
            if rss>limit:os.kill(os.getpid(),signal.SIGTERM);return
    def interrupted(signum,frame):raise MemoryError('Stationarity process stopped by sampled RSS guard or termination')
    signal.signal(signal.SIGTERM,interrupted);watcher=threading.Thread(target=monitor,daemon=True);watcher.start()
    try:
        result=evaluate(geometry['symbols'],geometry['coordinates_A'],geometry['A_indices'],
                        grid_coordinates,grid_weights,checkpoint,record['E_Ha'])
    except BaseException as error:
        atomic_write_json(out/'failure.json',dict(status='CPU_pilot_stationarity_failed_preserved',error=repr(error),
            sampled_peak_RSS_bytes=peak[0],source_sha256=before,qualified=False));raise
    finally:
        stop.set();watcher.join(timeout=3)
    api.verify_sources(before)
    passed=(result['reference_energy_error_Ha']<=1e-8 and result['checkpoint_energy_error_Ha']<=1e-8
            and result['orbital_gradient_norm_au']<=result['orbital_gradient_threshold_au'])
    result.update(status='CPU_pilot_energy_and_orbital_stationarity_checked' if passed else 'CPU_pilot_stationarity_failed_preserved',
        passed=passed,seconds=time.time()-start,source_sha256=before,sampled_peak_RSS_bytes=peak[0],
        full_counterpoise_computed=False,solvent_counterpoise_computed=False,
        scope='A-own frozen-fragment electronic state; no geometry gradient or complete interaction')
    atomic_write_json(out/'summary.json',result)
    if not passed:raise ValueError('Pilot checkpoint is not stationary at the declared threshold')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--water-out',type=Path)
    p.add_argument('--pilot-raw',type=Path);p.add_argument('--pilot-review',type=Path);p.add_argument('--pilot-out',type=Path)
    p.add_argument('--threads',type=int,default=2);a=p.parse_args()
    water=a.water_out is not None
    if water==bool(a.pilot_raw) or (not water and not all((a.pilot_raw,a.pilot_review,a.pilot_out))):
        p.error('Choose a water control or all three pilot paths')
    if a.threads<1:raise ValueError('Positive CPU thread count required')
    if pyscf.__version__!='2.14.0':raise ValueError('Qualified PySCF version required')
    with allocation_guard():host_capacity(a.threads,minimum_ram_gib=8 if water else 80);os.sched_setaffinity(0,cpu_allocation(a.threads))
    if any(os.environ.get(k)!=str(a.threads) for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')):
        raise ValueError('Declared CPU threads and capped native pools required')
    lib.num_threads(a.threads)
    result=check_water(a.water_out.resolve()) if water else check_pilot(a.pilot_raw,a.pilot_review,a.pilot_out)
    print(json.dumps(result))
