#!/usr/bin/env python
"""Review completed local single points and their limited resource evidence."""
import argparse
import json
import math
from pathlib import Path
import time

import numpy as np

from redox.core.common import ROOT,HARTREE_EV
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def review(campaign,out):
    plan=api.read_plan(campaign);api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    launch=json.loads((campaign/'launch.json').read_text())
    if launch['plan_sha256']!=api.digest(campaign/'plan.json'):raise ValueError('Executed local plan changed')
    path=campaign/'electronic/status.json';row=json.loads(path.read_text())
    ownership=api.load_module('local_electronic_stopped_owner',ROOT/'scripts/pipeline/run_polymer_spacing.py')
    if row['status']!='frozen_local_energy_finished_independent_review_pending' or ownership.owner_alive(row['owner']):
        raise ValueError('Completed stopped local electronic control required')
    api.verify_sources(row['output_sha256'])
    if row['primary_candidate_admission'] or row['thermochemistry_computed'] or plan['primary_candidate_admission']:
        raise ValueError('Local frozen energies cannot qualify a candidate or nuclear partition')
    raw=Path(row['raw_root']);symbols=[];coords=[]
    for line in Path(plan['input_xyz']).read_text().splitlines()[2:]:
        fields=line.split();symbols.append(fields[0]);coords.append([float(v) for v in fields[1:]])
    coords=np.asarray(coords);records={};checkpoints={}
    from pyscf.lib import chkfile
    for medium in ('gas','acetonitrile'):
        point=json.loads((raw/medium/'energy.json').read_text());settings=json.loads((raw/medium/'settings.json').read_text())
        if (not point['converged'] or point['medium']!=medium or point['charge']!=0 or point['mult']!=1
                or point['geometry_sha256']!=plan['geometry_sha256'] or point['settings']!=plan['settings']
                or settings['medium']!=medium or any(settings[k]!=v for k,v in plan['settings'].items())
                or point['geometry_optimized'] or point['thermochemistry_computed']
                or not math.isfinite(point['E_Ha']) or not math.isfinite(point['E_eV']) or abs(point['E_eV']-point['E_Ha']*HARTREE_EV)>1e-8
                or abs(point['orbital_diagnostics'].get('spin_contam',math.inf))>.1):
            raise ValueError('Local energy settings, charge/spin, geometry or scalar binding differs')
        checkpoint=raw/medium/'SCF.h5';mol=chkfile.load_mol(str(checkpoint));stored=chkfile.load(str(checkpoint),'scf')
        if (mol.charge!=0 or mol.spin!=0 or str(mol.basis).lower()!='def2-tzvpd'
                or [mol.atom_symbol(i) for i in range(mol.natm)]!=symbols
                or not np.allclose(mol.atom_coords(unit='Angstrom'),coords,atol=1e-10,rtol=0)
                or not math.isfinite(float(stored['e_tot'])) or abs(float(stored['e_tot'])-point['E_Ha'])>1e-10
                or abs(float(np.sum(stored['mo_occ']))-mol.nelectron)>1e-10):
            raise ValueError('Actual SCF checkpoint geometry, state, basis or energy differs')
        records[medium]=point;checkpoints[medium]=dict(actual_electrons=int(mol.nelectron),basis_functions=int(mol.nao_nr()),
            checkpoint_energy_error_Ha=abs(float(stored['e_tot'])-point['E_Ha']))
    guard=campaign/'residency_guard';owner=json.loads((guard/'owner.json').read_text());state=json.loads((guard/'state.json').read_text())
    host_path=campaign/'host_guard/status.json';host=json.loads(host_path.read_text());runtime_path=campaign/'tensor_runtime/status.json';rt=json.loads(runtime_path.read_text());profile=json.loads(Path(plan['runtime_profile']).read_text())
    if (not owner['guard_stopped'] or owner['pid']!=row['owner']['pid'] or not state['own_residency_observed']
            or row['actual_device']['uuid']!=launch['actual_device']['uuid']
            or state['foreign_resident_PIDs'] or (guard/'interruption.json').exists()
            or host['pid']!=row['owner']['pid'] or host['status']!='program_returned_review_required'
            or not 0<host['peak_RSS_bytes']<=plan['maximum_host_RSS_GiB']*2**30
            or not 0<row['observed_peak_GPU_MiB']<=plan['maximum_GPU_MiB']
            or rt['status']!='task_returned_under_checked_cutensor_runtime' or rt['pid']!=row['owner']['pid']
            or rt['program']!=plan['execution_driver'] or rt['versions']!=profile['versions']
            or rt['profile_sha256']!=api.digest(Path(plan['runtime_profile'])) or rt['library_sha256']!=profile['library_sha256']):
        raise ValueError('Local resource, exclusive residency or contraction runtime gate failed')
    proof={**plan['producer_sha256'],**plan['input_sha256'],**row['output_sha256']}
    for p in (Path(__file__).resolve(),campaign/'plan.json',campaign/'launch.json',path,guard/'owner.json',guard/'state.json',host_path,runtime_path):proof[str(p)]=api.digest(p)
    result=dict(status='frozen_local_electronic_state_and_resource_pilot_independently_checked',checked=time.time(),
        campaign=str(campaign),task=plan['task'],settings=plan['settings'],energies={k:v['E_eV'] for k,v in records.items()},
        SCF_checkpoint_checks=checkpoints,observed_peak_GPU_MiB=row['observed_peak_GPU_MiB'],
        peak_host_RSS_bytes=host['peak_RSS_bytes'],actual_device=row['actual_device'],source_sha256=proof,
        geometry_optimized=False,thermochemistry_computed=False,primary_candidate_admission=False,
        interaction_energy_assembled=False,basis_superposition_correction_computed=False,
        large_polymer_or_Hessian_resource_qualification=False,
        scope='Numerical fixed-geometry local monomer/pair energy and measured resource point; no interaction/free-energy/coupling or large-system admission')
    out.parent.mkdir(parents=True,exist_ok=True);atomic_write_json(out,result);return result


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--campaign',type=Path,required=True);parser.add_argument('--out',type=Path,required=True);a=parser.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    import os
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    row=review(a.campaign.resolve(),a.out.resolve());print(json.dumps({k:row[k] for k in ['status','observed_peak_GPU_MiB']}))


if __name__=='__main__':main()
