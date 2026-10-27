#!/usr/bin/env python3
"""Independently inspect a completed CPU counterpoise resource-pilot checkpoint.

This checks identity, stored SCF content and sampled resource evidence; it does
not recompute the Fock gradient or qualify a five-term interaction.
"""
import argparse
import json
import math
import os
from pathlib import Path
import time

import numpy as np
from pyscf import df
from pyscf.scf.chkfile import load_scf

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity


def verify_orbitals(coefficients,occupations,orbital_energies,overlap,electrons,tolerance=1e-7):
    coeff=np.asarray(coefficients);occ=np.asarray(occupations);energy=np.asarray(orbital_energies);overlap=np.asarray(overlap)
    if overlap.ndim!=2:raise ValueError('Overlap matrix must be two dimensional')
    n=overlap.shape[0]
    if coeff.ndim!=2:raise ValueError('Checkpoint orbital matrix must be two dimensional')
    m=coeff.shape[1]
    if (overlap.shape!=(n,n) or coeff.shape[0]!=n or m<electrons//2 or occ.shape!=(m,) or energy.shape!=(m,)
            or not all(np.isfinite(a).all() for a in (coeff,occ,energy,overlap))
            or not np.all(np.isin(occ,[0.,2.])) or abs(float(occ.sum())-electrons)>1e-8):
        raise ValueError('Nonfinite or inconsistent closed-shell checkpoint orbitals')
    error=float(np.max(np.abs(coeff.T@overlap@coeff-np.eye(m))))
    if error>tolerance:raise ValueError('Checkpoint orbitals are not overlap-orthonormal')
    return error


def verify_owner_stopped(owner):
    if owner['host']!=os.uname().nodename:raise ValueError('Cannot establish remote owner status')
    stat=Path('/proc')/str(owner['pid'])/'stat'
    if stat.exists():
        fields=stat.read_text().rsplit(')',1)[1].split()
        if fields[19]==owner['process_start'] and fields[0] not in ('Z','X'):
            raise ValueError('Pilot process is still live')


def review(raw,out):
    raw=raw.resolve();out=out.resolve()
    if out.exists():raise FileExistsError('Preserve prior reviews')
    if not raw.is_relative_to((ROOT/'calcs').resolve()):raise ValueError('Pilot must reside on project storage')
    launch=json.loads((raw/'launch.json').read_text());status=json.loads((raw/'status.json').read_text())
    if status.get('status')!='CPU_pilot_converged_independent_review_pending':
        raise ValueError('Completed but unqualified one-term pilot required')
    role=raw/'A_own';record=json.loads((role/'energy.json').read_text())
    if (status!=record or record['status']!='CPU_pilot_converged_independent_review_pending'
            or record['role']!='A_own' or not record['converged'] or record['qualified']
            or record['full_counterpoise_computed'] or record['primary_candidate_admission']):
        raise ValueError('Completed but unqualified one-term pilot required')
    verify_owner_stopped(record['owner'])
    if record['source_sha256']!=launch['source_sha256']:
        raise ValueError('Pilot source identity changed')
    api.verify_sources(record['source_sha256'])
    case=Path(launch['case']);plan=api.read_plan(case);geometry=json.loads(Path(plan['input_geometry']).read_text())
    core=api.load_module('review_cpu_pilot_counterpoise',ROOT/'scripts/validation/local_counterpoise.py')
    if record['settings']!=plan['settings'] or plan['settings']!=core.SETTINGS:
        raise ValueError('Electronic settings differ from reviewed plan')
    if (record['peak_RSS_bytes']<=0 or record['peak_RSS_bytes']>32*1024**3
            or record['electrons']!=84 or record['basis_functions']!=513):
        raise ValueError('Resource or electronic identity is out of scope')
    mols=core.molecules(geometry['symbols'],geometry['coordinates_A'],geometry['A_indices'])
    mol=mols['A_own'];aux=df.addons.make_auxmol(mols['AB'],None)
    if aux.nao_nr()!=record['auxiliary_functions']:
        raise ValueError('Auxiliary basis size differs from full-complex basis')
    if set(record['common_grid_hashes'])!={'grid_coordinates_AU.npy','grid_weights.npy'}:
        raise ValueError('Saved integration grid inventory differs')
    for name,expected in record['common_grid_hashes'].items():
        if api.digest(raw/name)!=expected:
            raise ValueError('Saved integration grid changed')
    coords=np.load(raw/'grid_coordinates_AU.npy',mmap_mode='r');weights=np.load(raw/'grid_weights.npy',mmap_mode='r')
    if (coords.shape!=(record['grid_points'],3) or weights.shape!=(record['grid_points'],)
            or not np.isfinite(coords).all() or not np.isfinite(weights).all()):
        raise ValueError('Saved integration grid has invalid dimensions or values')
    checkpoint=role/'SCF.h5'
    if api.digest(checkpoint)!=record['checkpoint_sha256']:
        raise ValueError('Checkpoint hash differs')
    saved,scf=load_scf(str(checkpoint))
    if (saved.nao_nr()!=mol.nao_nr() or saved.nelectron!=mol.nelectron
            or not np.array_equal(saved.atom_charges(),mol.atom_charges())
            or not np.allclose(saved.atom_coords(),mol.atom_coords(),rtol=0,atol=1e-10)):
        raise ValueError('Checkpoint molecule differs from planned fragment')
    energy=float(scf['e_tot'])
    if not math.isfinite(energy) or abs(energy-record['E_Ha'])>1e-9:
        raise ValueError('Checkpoint and recorded energies differ')
    error=verify_orbitals(scf['mo_coeff'],scf['mo_occ'],scf['mo_energy'],mol.intor_symmetric('int1e_ovlp'),mol.nelectron)
    paths=[raw/'launch.json',raw/'status.json',role/'energy.json',checkpoint,
           raw/'grid_coordinates_AU.npy',raw/'grid_weights.npy',Path(__file__).resolve()]
    result=dict(status='CPU_pilot_checkpoint_and_sampled_resources_reviewed_gradient_pending',checked=time.time(),
        owner=record['owner'],role='A_own',energy_Ha=energy,electrons=mol.nelectron,basis_functions=mol.nao_nr(),
        auxiliary_functions=aux.nao_nr(),grid_points=record['grid_points'],orbital_orthonormality_max_error=error,
        sampled_peak_RSS_bytes=record['peak_RSS_bytes'],source_sha256={str(p):api.digest(p) for p in paths},
        physical_gradient_recomputed=False,full_counterpoise_computed=False,primary_candidate_admission=False,
        scope='Checkpoint identity and sampled resources only; no independently recomputed SCF gradient or full interaction')
    out.parent.mkdir(parents=True,exist_ok=True);atomic_write_json(out,result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    result=review(a.raw,a.out);print(json.dumps(dict(status=result['status'],energy_Ha=result['energy_Ha'])))
