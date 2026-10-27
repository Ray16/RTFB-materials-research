#!/usr/bin/env python3
"""Review a completed remaining CPU counterpoise term without promoting it."""
import argparse
import json
import math
import os
from pathlib import Path
import time

import numpy as np
import pyscf
from pyscf import df
from pyscf.scf.chkfile import load_scf

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

EXPECTED={'AB':(168,1026,0),'B_own':(84,513,0),'B_ghost':(84,1026,19)}


def review(raw,out):
    raw=raw.resolve();out=out.resolve()
    if out.exists():raise FileExistsError('Preserve prior reviews')
    if not raw.is_relative_to((ROOT/'calcs').resolve()):raise ValueError('Term must reside on durable project storage')
    launch=json.loads((raw/'launch.json').read_text());status=json.loads((raw/'status.json').read_text());role=launch['role']
    if launch['versions']!={'pyscf':pyscf.__version__,'numpy':np.__version__}:raise ValueError('Term and review environments differ')
    if role not in EXPECTED or status.get('status')!='CPU_remaining_term_converged_independent_review_pending':
        raise ValueError('Completed unqualified remaining term required')
    folder=raw/role;record=json.loads((folder/'energy.json').read_text())
    if (record!=status or record['role']!=role or not record['converged'] or record['qualified']
            or record['full_counterpoise_computed'] or record['primary_candidate_admission']
            or record['source_sha256']!=launch['source_sha256']):raise ValueError('Term result or source identity differs')
    helper=api.load_module('remaining_review_common_checks',ROOT/'scripts/validation/review_counterpoise_cpu_pilot.py')
    helper.verify_owner_stopped(record['owner']);api.verify_sources(record['source_sha256'])
    plan=api.read_plan(Path(launch['case']));geometry=json.loads(Path(plan['input_geometry']).read_text())
    core=api.load_module('remaining_review_protocol',ROOT/'scripts/validation/local_counterpoise.py')
    if record['settings']!=plan['settings'] or plan['settings']!=core.SETTINGS:raise ValueError('Term electronic method differs')
    electrons,basis,ghosts=EXPECTED[role]
    if (record['peak_RSS_bytes']<=0 or record['peak_RSS_bytes']>32*1024**3 or record['electrons']!=electrons
            or record['basis_functions']!=basis or record['ghost_centers']!=ghosts):raise ValueError('Term identity or sampled resources differ')
    binder=api.load_module('remaining_review_launch_arguments',ROOT/'scripts/validation/review_counterpoise_cpu_ghost.py')
    own=binder.argument(launch['command'],'--own')
    if set(record['common_grid_hashes'])!={'grid_coordinates_AU.npy','grid_weights.npy'}:raise ValueError('Common-grid inventory differs')
    for name,digest in record['common_grid_hashes'].items():
        if api.digest(raw/name)!=digest or api.digest(own/name)!=digest:raise ValueError('Term grid differs from A-own reference')
    coords=np.load(raw/'grid_coordinates_AU.npy',mmap_mode='r');weights=np.load(raw/'grid_weights.npy',mmap_mode='r')
    if (coords.shape!=(record['grid_points'],3) or weights.shape!=(record['grid_points'],)
            or not np.isfinite(coords).all() or not np.isfinite(weights).all()):raise ValueError('Invalid common integration grid')
    checkpoint=folder/'SCF.h5'
    if api.digest(checkpoint)!=record['checkpoint_sha256']:raise ValueError('Checkpoint hash differs')
    mols=core.molecules(geometry['symbols'],geometry['coordinates_A'],geometry['A_indices']);mol=mols[role];full=mols['AB'];saved,scf=load_scf(str(checkpoint))
    if (saved.nao_nr()!=mol.nao_nr() or saved.nelectron!=mol.nelectron or saved.natm!=mol.natm
            or not np.array_equal(saved.atom_charges(),mol.atom_charges())
            or not np.allclose(saved.atom_coords(),mol.atom_coords(),rtol=0,atol=1e-10)):
        raise ValueError('Checkpoint nuclei, basis or geometry differs')
    aux=df.addons.make_auxmol(full,None)
    if aux.nao_nr()!=record['auxiliary_functions']:raise ValueError('Full-complex auxiliary basis differs')
    energy=float(scf['e_tot'])
    if not math.isfinite(energy) or abs(energy-record['E_Ha'])>1e-9:raise ValueError('Checkpoint and recorded energies differ')
    error=helper.verify_orbitals(scf['mo_coeff'],scf['mo_occ'],scf['mo_energy'],mol.intor_symmetric('int1e_ovlp'),mol.nelectron)
    paths=[raw/'launch.json',raw/'status.json',folder/'energy.json',checkpoint,raw/'grid_coordinates_AU.npy',raw/'grid_weights.npy',
           Path(__file__).resolve(),ROOT/'scripts/validation/review_counterpoise_cpu_pilot.py',ROOT/'scripts/validation/review_counterpoise_cpu_ghost.py']
    result=dict(status='CPU_remaining_term_checkpoint_and_sampled_resources_reviewed_gradient_pending',checked=time.time(),owner=record['owner'],
        role=role,energy_Ha=energy,electrons=electrons,basis_functions=basis,ghost_centers=ghosts,auxiliary_functions=aux.nao_nr(),
        grid_points=record['grid_points'],orbital_orthonormality_max_error=error,sampled_peak_RSS_bytes=record['peak_RSS_bytes'],
        source_sha256={str(p):api.digest(p) for p in paths},physical_gradient_recomputed=False,full_counterpoise_computed=False,
        primary_candidate_admission=False,scope='Checkpoint identity and sampled resources only; no recomputed SCF gradient or full interaction')
    out.parent.mkdir(parents=True,exist_ok=True);atomic_write_json(out,result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    result=review(a.raw,a.out);print(json.dumps(dict(status=result['status'],role=result['role'],energy_Ha=result['energy_Ha'])))
