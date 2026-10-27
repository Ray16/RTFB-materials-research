#!/usr/bin/env python3
"""Recompute the saved water-dimer ghost-basis energy and orbital gradient."""
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


def verify_ghost_checkpoint(saved,mol,full,expected_ghost_centers):
    if (saved.nao_nr()!=full.nao_nr() or saved.nao_nr()!=mol.nao_nr() or saved.nelectron!=mol.nelectron
            or sum(saved.atom_charges()==0)!=expected_ghost_centers
            or not np.array_equal(saved.atom_charges(),mol.atom_charges())
            or not np.allclose(saved.atom_coords(),mol.atom_coords(),rtol=0,atol=1e-10)):
        raise ValueError('Ghost checkpoint nuclei, basis or geometry differs')


def check(out):
    if out.exists():raise FileExistsError('Preserve earlier ghost controls')
    review_path=ROOT/'results/reviews/2026-10-26/counterpoise_control/independent_review.json'
    review=json.loads(review_path.read_text());api.verify_sources(review['source_sha256'])
    if review['status']!='gas_counterpoise_CPU_GPU_control_independently_checked':raise ValueError('Reviewed control required')
    core_path=ROOT/'scripts/validation/local_counterpoise.py';factory_path=ROOT/'scripts/validation/run_cpu_counterpoise_reference.py'
    core=api.load_module('ghost_water_counterpoise_protocol',core_path);fixed=api.load_module('ghost_water_fixed_auxiliary',factory_path)
    if core.SETTINGS!=review['settings']:raise ValueError('Control protocol differs')
    base=ROOT/'calcs/validation/2026-10-26_counterpoise/CPU_fixed_auxiliary_factory'
    summary_path=base/'summary.json';energy_path=base/'A_ghost/energy.json';checkpoint=base/'A_ghost/SCF.h5'
    paths=[summary_path,energy_path,checkpoint,base/'grid_coordinates_AU.npy',base/'grid_weights.npy',core_path,factory_path,review_path,Path(__file__).resolve()]
    for p in paths[:5]:
        if review['source_sha256'].get(str(p))!=api.digest(p):raise ValueError('Reviewed source changed')
    source={str(p):api.digest(p) for p in paths};baseline=json.loads(summary_path.read_text());record=json.loads(energy_path.read_text())
    mols=core.molecules(baseline['source_symbols'],baseline['source_coordinates_A'],baseline['A_indices']);full=mols['AB'];mol=mols['A_ghost']
    saved,state=load_scf(str(checkpoint))
    verify_ghost_checkpoint(saved,mol,full,expected_ghost_centers=3)
    coords=np.load(paths[3]);weights=np.load(paths[4]);mf=dft.RKS(mol).density_fit(only_dfj=False);aux=df.addons.make_auxmol(full,None)
    mf.with_df.auxmol=aux.copy();mf.xc=core.SETTINGS['xc'];mf.nlc=core.SETTINGS['nlc'];mf.small_rho_cutoff=0.;mf.grids.level=3;mf.nlcgrids.level=3
    for grids in (mf.grids,mf.nlcgrids):grids.coords=coords.copy();grids.weights=weights.copy()
    out.mkdir(parents=True);started=time.time();dm=mf.make_rdm1(state['mo_coeff'],state['mo_occ']);h=mf.get_hcore();overlap=mf.get_ovlp()
    with fixed.fixed_complex_auxiliary_factory(df.addons,full,mols.values()):
        veff=mf.get_veff(mol,dm);fock=mf.get_fock(h1e=h,s1e=overlap,vhf=veff,dm=dm)
        energy=float(mf.energy_tot(dm=dm,h1e=h,vhf=veff));gradient=np.asarray(mf.get_grad(state['mo_coeff'],state['mo_occ'],fock))
    actual=mf.with_df.auxmol
    if (actual.nao_nr()!=aux.nao_nr() or not np.array_equal(actual.atom_charges(),aux.atom_charges())
            or not np.allclose(actual.atom_coords(),aux.atom_coords(),rtol=0,atol=1e-10)):
        raise ValueError('Full-complex auxiliary basis was not retained')
    api.verify_sources(source);norm=float(np.linalg.norm(gradient));error=abs(energy-record['E_Ha'])
    passed=math.isfinite(energy) and np.isfinite(gradient).all() and error<=1e-8 and norm<=core.SETTINGS['SCF_conv_tol_grad']
    result=dict(status='water_ghost_energy_and_orbital_stationarity_checked' if passed else 'water_ghost_control_failed_preserved',
        passed=bool(passed),role='A_ghost',energy_Ha=energy,reference_energy_Ha=record['E_Ha'],energy_error_Ha=error,
        orbital_gradient_norm_au=norm,orbital_gradient_threshold_au=core.SETTINGS['SCF_conv_tol_grad'],
        electrons=int(mol.nelectron),basis_functions=int(mol.nao_nr()),ghost_centers=int(sum(mol.atom_charges()==0)),
        auxiliary_functions=int(actual.nao_nr()),seconds=time.time()-started,source_sha256=source,
        scope='CPU water-dimer ghost-basis state only; no large-motif qualification')
    atomic_write_json(out/'summary.json',result)
    if not passed:raise ValueError('Water ghost control failed stationarity')
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    with allocation_guard():host_capacity(2);os.sched_setaffinity(0,cpu_allocation(2))
    if pyscf.__version__!='2.14.0' or any(os.environ.get(k)!='2' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')):
        raise ValueError('Qualified PySCF with two capped CPU threads required')
    lib.num_threads(2);print(json.dumps(check(a.out.resolve())))
