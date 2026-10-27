#!/usr/bin/env python3
"""Review a finished CPU A-ghost pilot without claiming its stationarity."""
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


def argument(command,flag):
    if command.count(flag)!=1 or command.index(flag)+1>=len(command):raise ValueError('Missing or duplicated launch argument')
    return Path(command[command.index(flag)+1]).resolve()


def review(raw,out):
    raw=raw.resolve();out=out.resolve()
    if out.exists():raise FileExistsError('Preserve prior reviews')
    if not raw.is_relative_to((ROOT/'calcs').resolve()):raise ValueError('Ghost pilot must use shared project storage')
    launch=json.loads((raw/'launch.json').read_text());status=json.loads((raw/'status.json').read_text())
    if launch['versions']!={'pyscf':pyscf.__version__,'numpy':np.__version__}:
        raise ValueError('Ghost review and calculation environments differ')
    if status.get('status')!='CPU_ghost_pilot_converged_independent_review_pending':
        raise ValueError('Completed but unqualified ghost pilot required')
    role=raw/'A_ghost';record=json.loads((role/'energy.json').read_text())
    if (record!=status or not record['converged'] or record['role']!='A_ghost' or record['qualified']
            or record['full_counterpoise_computed'] or record['primary_candidate_admission']
            or record['source_sha256']!=launch['source_sha256']):
        raise ValueError('Ghost pilot result or source identity differs')
    own_reviewer=api.load_module('ghost_review_common_checks',ROOT/'scripts/validation/review_counterpoise_cpu_pilot.py')
    own_reviewer.verify_owner_stopped(record['owner']);api.verify_sources(record['source_sha256'])
    case=Path(launch['case']);plan=api.read_plan(case);geometry=json.loads(Path(plan['input_geometry']).read_text())
    core=api.load_module('ghost_review_counterpoise_protocol',ROOT/'scripts/validation/local_counterpoise.py')
    if record['settings']!=plan['settings'] or plan['settings']!=core.SETTINGS:
        raise ValueError('Ghost electronic settings differ')
    if (record['peak_RSS_bytes']<=0 or record['peak_RSS_bytes']>32*1024**3
            or record['electrons']!=84 or record['basis_functions']!=1026 or record['ghost_centers']!=19):
        raise ValueError('Ghost pilot identity or sampled resources differ')
    own=argument(launch['command'],'--own');stationarity=argument(launch['command'],'--stationarity')
    gradient=json.loads(stationarity.read_text())
    if (gradient['status']!='CPU_pilot_energy_and_orbital_stationarity_checked' or not gradient['passed']
            or gradient['source_sha256'].get(str(own/'A_own/energy.json'))!=api.digest(own/'A_own/energy.json')):
        raise ValueError('Ghost term is not bound to the stationary A-own reference')
    api.verify_sources(gradient['source_sha256'])
    if set(record['common_grid_hashes'])!={'grid_coordinates_AU.npy','grid_weights.npy'}:
        raise ValueError('Ghost common-grid inventory differs')
    for name,digest in record['common_grid_hashes'].items():
        if api.digest(raw/name)!=digest or api.digest(own/name)!=digest:
            raise ValueError('Ghost grid does not match the A-own reference')
    coords=np.load(raw/'grid_coordinates_AU.npy',mmap_mode='r');weights=np.load(raw/'grid_weights.npy',mmap_mode='r')
    if (coords.shape!=(record['grid_points'],3) or weights.shape!=(record['grid_points'],)
            or not np.isfinite(coords).all() or not np.isfinite(weights).all()):raise ValueError('Invalid ghost integration grid')
    checkpoint=role/'SCF.h5'
    if api.digest(checkpoint)!=record['checkpoint_sha256']:raise ValueError('Ghost checkpoint hash differs')
    mols=core.molecules(geometry['symbols'],geometry['coordinates_A'],geometry['A_indices']);mol=mols['A_ghost'];full=mols['AB']
    saved,scf=load_scf(str(checkpoint));ghost_control=api.load_module('ghost_review_ghost_identity',ROOT/'scripts/validation/check_cpu_counterpoise_ghost_control.py')
    ghost_control.verify_ghost_checkpoint(saved,mol,full,expected_ghost_centers=19)
    aux=df.addons.make_auxmol(full,None)
    if aux.nao_nr()!=record['auxiliary_functions']:raise ValueError('Full-complex auxiliary basis size differs')
    energy=float(scf['e_tot'])
    if not math.isfinite(energy) or abs(energy-record['E_Ha'])>1e-9:raise ValueError('Checkpoint and recorded energies differ')
    error=own_reviewer.verify_orbitals(scf['mo_coeff'],scf['mo_occ'],scf['mo_energy'],mol.intor_symmetric('int1e_ovlp'),mol.nelectron)
    paths=[raw/'launch.json',raw/'status.json',role/'energy.json',checkpoint,raw/'grid_coordinates_AU.npy',raw/'grid_weights.npy',
           stationarity,Path(__file__).resolve(),ROOT/'scripts/validation/review_counterpoise_cpu_pilot.py',ROOT/'scripts/validation/check_cpu_counterpoise_ghost_control.py']
    result=dict(status='CPU_ghost_pilot_checkpoint_and_sampled_resources_reviewed_gradient_pending',checked=time.time(),
        owner=record['owner'],role='A_ghost',energy_Ha=energy,electrons=mol.nelectron,basis_functions=mol.nao_nr(),
        ghost_centers=19,auxiliary_functions=aux.nao_nr(),grid_points=record['grid_points'],orbital_orthonormality_max_error=error,
        sampled_peak_RSS_bytes=record['peak_RSS_bytes'],source_sha256={str(p):api.digest(p) for p in paths},
        physical_gradient_recomputed=False,full_counterpoise_computed=False,primary_candidate_admission=False,
        scope='Ghost checkpoint identity and sampled resources only; no independently recomputed SCF gradient or full interaction')
    out.parent.mkdir(parents=True,exist_ok=True);atomic_write_json(out,result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--raw',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    result=review(a.raw,a.out);print(json.dumps(dict(status=result['status'],energy_Ha=result['energy_Ha'])))
