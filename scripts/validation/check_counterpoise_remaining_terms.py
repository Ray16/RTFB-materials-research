#!/usr/bin/env python3
"""Recompute saved CPU counterpoise term energies and orbital gradients.

The first use checks the remaining three terms of a reviewed water-dimer
control. Frozen-geometry orbital stationarity is not geometry stationarity,
solution thermochemistry or a polymer material assessment.
"""
import argparse
import json
import math
import os
from pathlib import Path
import time

import numpy as np
import pyscf
from pyscf import dft,df,lib
from pyscf.scf.chkfile import load_scf

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

CORE=ROOT/'scripts/validation/local_counterpoise.py'
FACTORY=ROOT/'scripts/validation/run_cpu_counterpoise_reference.py'
ROLES=('AB','B_own','B_ghost')


def evaluate(geometry,role,checkpoint,coords_path,weights_path,reference_energy):
    if role not in ROLES:raise ValueError('Unsupported remaining counterpoise role')
    core=api.load_module('remaining_counterpoise_protocol',CORE);fixed=api.load_module('remaining_counterpoise_auxiliary',FACTORY)
    mols=core.molecules(geometry['symbols'],geometry['coordinates_A'],geometry['A_indices']);full=mols['AB'];mol=mols[role]
    saved,state=load_scf(str(checkpoint))
    if (saved.nao_nr()!=mol.nao_nr() or saved.nelectron!=mol.nelectron or saved.natm!=mol.natm
            or not np.array_equal(saved.atom_charges(),mol.atom_charges())
            or not np.allclose(saved.atom_coords(),mol.atom_coords(),rtol=0,atol=1e-10)):
        raise ValueError('Checkpoint nuclei, basis or geometry differs from requested role')
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
    if not math.isfinite(energy) or not np.isfinite(gradient).all():raise ValueError('Nonfinite energy or orbital gradient')
    return dict(role=role,energy_Ha=energy,reference_energy_Ha=reference_energy,
        energy_error_Ha=abs(energy-reference_energy),checkpoint_energy_error_Ha=abs(energy-float(state['e_tot'])),
        orbital_gradient_norm_au=float(np.linalg.norm(gradient)),orbital_gradient_threshold_au=core.SETTINGS['SCF_conv_tol_grad'],
        electrons=int(mol.nelectron),basis_functions=int(mol.nao_nr()),ghost_centers=int(sum(mol.atom_charges()==0)),
        auxiliary_functions=int(actual.nao_nr()))


def check_water(out):
    if out.exists():raise FileExistsError('Preserve previous water-control attempts')
    review_path=ROOT/'results/reviews/2026-10-26/counterpoise_control/independent_review.json'
    review=json.loads(review_path.read_text());api.verify_sources(review['source_sha256'])
    if review['status']!='gas_counterpoise_CPU_GPU_control_independently_checked':raise ValueError('Reviewed water control required')
    core=api.load_module('remaining_control_settings',CORE)
    if review['settings']!=core.SETTINGS:raise ValueError('Water-control protocol differs')
    base=ROOT/'calcs/validation/2026-10-26_counterpoise/CPU_fixed_auxiliary_factory';summary_path=base/'summary.json'
    paths=[summary_path,base/'grid_coordinates_AU.npy',base/'grid_weights.npy',CORE,FACTORY,review_path,Path(__file__).resolve()]
    for role in ROLES:paths.extend([base/role/'energy.json',base/role/'SCF.h5'])
    for path in paths[:3]+paths[7:]:
        if review['source_sha256'].get(str(path))!=api.digest(path):raise ValueError('Reviewed water input changed')
    sources={str(p):api.digest(p) for p in paths};baseline=json.loads(summary_path.read_text())
    geometry=dict(symbols=baseline['source_symbols'],coordinates_A=baseline['source_coordinates_A'],A_indices=baseline['A_indices'])
    out.mkdir(parents=True);results=[];started=time.time()
    for role in ROLES:
        record=json.loads((base/role/'energy.json').read_text())
        result=evaluate(geometry,role,base/role/'SCF.h5',base/'grid_coordinates_AU.npy',base/'grid_weights.npy',record['E_Ha'])
        passed=(result['energy_error_Ha']<=1e-8 and result['checkpoint_energy_error_Ha']<=1e-8
                and result['orbital_gradient_norm_au']<=result['orbital_gradient_threshold_au'])
        result['passed']=passed;results.append(result);atomic_write_json(out/(role+'.json'),result)
        if not passed:raise ValueError('Water term fails energy or stationarity: '+role)
    api.verify_sources(sources)
    report=dict(status='remaining_water_counterpoise_terms_stationary',passed=True,terms=results,seconds=time.time()-started,
        source_sha256=sources,scope='Three water-dimer electronic terms; no large-motif resource or polymer qualification')
    atomic_write_json(out/'summary.json',report);return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--water-out',type=Path,required=True);a=p.parse_args()
    if pyscf.__version__!='2.14.0':raise ValueError('Qualified PySCF required')
    with allocation_guard():host_capacity(2);os.sched_setaffinity(0,cpu_allocation(2))
    if any(os.environ.get(k)!='2' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')):
        raise ValueError('Two capped CPU threads required')
    lib.num_threads(2);result=check_water(a.water_out.resolve());print(json.dumps(result))
