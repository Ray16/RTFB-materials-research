#!/usr/bin/env python3
"""Independently reconstruct one closed-shell charged polymer component G."""
import argparse
import importlib.util
import json
import math
import os
from pathlib import Path
import re

import h5py
import numpy as np
from ase.io import read
from pyscf import gto

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify


def load(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def component_sums(component,nuclear,energy_Ha):
    if component['mult']!=1 or component['charge']<=0 or component['spin_partition_eV']!=0. or component['primary_model_admission'] or component['material_admission']:
        raise ValueError('Explicit alternate-model closed-shell ionic component required')
    energy=energy_Ha*27.211386245988
    if not all(math.isfinite(v) for v in (energy,component['E_geometry_level_eV'],component['E_geometry_level_Ha'])) or abs(component['E_geometry_level_eV']-energy)>1e-8 or abs(component['E_geometry_level_Ha']-energy_Ha)>1e-10:
        raise ValueError('Geometry-level electronic parent differs')
    totals={}
    for key in ('25','50','100'):
        value=energy+nuclear['nuclear_by_cutoff_eV'][key]
        if (not all(math.isfinite(v) for v in (value,component['nuclear_by_cutoff_eV'][key],component['G_by_cutoff_eV'][key])) or abs(component['nuclear_by_cutoff_eV'][key]-nuclear['nuclear_by_cutoff_eV'][key])>1e-10
                or abs(component['G_by_cutoff_eV'][key]-value)>1e-8):raise ValueError('Exact electronic/nuclear/spin component sum does not close')
        totals[key]=value
    return totals


def review(namespace):
    spec=verify(namespace);case=Path(spec['case']);sources=Path(spec['sources']);used=[namespace]
    def read_json(path):used.append(path);return json.loads(path.read_text())
    plan=read_json(case/'plan.json');stage=read_json(case/'status.json');placements=list(case.parent.glob('placement*/launch.json'))
    if len(placements)!=1:raise ValueError('One exact charged native placement required')
    launch=read_json(placements[0]);completion=read_json(placements[0].parent/'completion.json')
    owner=read_json(case/'residency_guard/owner.json');residency=read_json(case/'residency_guard/state.json')
    host=read_json(case/'host_guard/status.json');runtime=read_json(case/'tensor_runtime/status.json')
    owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','charged_native_closed_owners')
    if (not plan.get('charged_polymer_native_component') or plan['role']!='charged_polymer_native_component'
            or plan['mult']!=1 or plan['charge']<=0 or stage['status']!='charged_polymer_native_component_finished_review_pending'
            or stage['native_G_qualified'] or owners.alive(stage['owner']) is not False
            or completion['returncode']!=0 or completion['budget_stop'] is not None or not owner.get('guard_stopped')
            or not residency['own_residency_observed'] or residency['foreign_resident_PIDs']
            or host['status']!='effective_memory_guard_finished_profile_review_pending'
            or runtime['status']!='task_returned_under_checked_cutensor_runtime'
            or host['maximum_RSS_GiB']!=80 or host['minimum_free_GiB']!=32 or host['free_RAM_bytes']<32*2**30
            or max(host['peak_tree_RSS_bytes'],host['self_VmHWM_bytes'])>80*2**30
            or not 0<stage['observed_peak_GPU_MiB']<=31744
            or not stage['owner']['pid']==owner['pid']==host['pid']==runtime['pid']
            or not stage['owner']['process_start']==owner['process_start']==host['process_start']
            or not stage['owner']['host']==owner['host']==host['host']==launch['host']
            or not owner['CPU_affinity']==host['CPU_affinity']==runtime['CPU_affinity']==launch['CPU_affinity']):
        raise ValueError('Normal closed native component and exact resource/ownership chain required')
    source=read_json(Path(plan['source_stationary_review']))
    if (source['status']!='charged_polymer_solution_stationary_geometry_independently_checked'
            or source['graph_key']!=plan['prepared_case']['key'] or source['candidate_alias']!=plan['exact_candidate_alias']
            or source['charge']!=plan['charge'] or source['mult']!=plan['mult'] or source['solvent_model']!=plan['solvent_model']):
        raise ValueError('Stationary source identity/model differs')
    geometry_path=sources/'input.xyz';atoms=read(geometry_path);used.append(geometry_path);count=len(atoms)
    matrix_path=case/'Hessian.npy';matrix=np.load(matrix_path,allow_pickle=False);used.append(matrix_path)
    center=read_json(case/'points/0000/evaluation.json');coordinates=np.asarray(center['coordinates_Bohr']);gradient=np.asarray(center['gradient_Ha_Bohr'])
    if (atoms.get_chemical_symbols()!=plan['symbols'] or matrix.dtype!=np.float64 or matrix.shape!=(3*count,3*count) or not np.isfinite(matrix).all()
            or abs(matrix-matrix.T).max()>1e-6 or gradient.shape!=(count,3) or not np.isfinite(gradient).all()
            or abs(gradient).max()>.00045 or not np.allclose(coordinates,atoms.positions/.529177210903,rtol=0,atol=1e-10)
            or abs(center['E_Ha']-source['E_Ha'])>1e-8):raise ValueError('Native matrix, stationarity or source geometry/energy fails')
    hetero=next((i for i,s in enumerate(plan['symbols']) if s=='O'),next((i for i,s in enumerate(plan['symbols']) if s=='Fe'),count-1))
    required={(0,1),(hetero,2),*map(tuple,plan['extra_FD_columns'])}
    required.update(tuple(map(int,np.unravel_index(i,gradient.shape))) for i in np.argsort(abs(gradient).ravel())[-2:])
    columns=read_json(case/'columns.json');derivatives=read_json(case/'derivatives.json')
    if not derivatives['complete'] or set(map(tuple,columns['columns']))!=required:raise ValueError('Full declared sampled derivative coverage required')
    checker=load(ROOT/'scripts/validation/review_screening_polymer_hessian.py','charged_native_independent_columns')
    errors=checker.column_checks(matrix,center,required,derivatives['checks'],read_json)
    points=sorted((case/'points').glob('*/evaluation.json'))
    if len(points)!=1+4*len(required):raise ValueError('Native SCF point inventory differs')
    for point_path in points:
        point=read_json(point_path);checkpoint=Path(point['SCF_checkpoint']);settings=read_json(checkpoint.with_name('settings.json'));solved=read_json(checkpoint.with_name('SCF_result.json'));used.append(checkpoint)
        if (not point['SCF_converged'] or abs(point['spin'].get('spin_contam',math.inf))>.1
                or any(settings[key]!=value for key,value in dict(charge=plan['charge'],mult=1,basis='def2-svp',xc='r2scan',disp='d4',conv_tol=1e-12,conv_tol_grad=1e-8,max_cycle=200,grid_level=3).items())
                or not np.allclose(settings['coordinates_Bohr'],point['coordinates_Bohr'],rtol=0,atol=1e-12)
                or not solved['converged'] or abs(solved['energy_Ha']-point['E_Ha'])>1e-10):raise ValueError('Actual charged native SCF/state differs')
        with h5py.File(checkpoint,'r') as stream:
            mol=gto.loads(stream['mol'][()])
            if (mol.charge!=plan['charge'] or mol.spin!=0 or mol.basis!='def2-svp'
                    or not np.allclose(mol.atom_coords(),point['coordinates_Bohr'],rtol=0,atol=1e-12)
                    or abs(float(stream['scf/e_tot'][()])-point['E_Ha'])>1e-10):raise ValueError('Native SCF checkpoint differs')
    solver=read_json(case/'execution/response_solver.json')
    if (solver['status']!='response_equations_solved_full_physical_QC_separate' or solver['equations_changed'] or solver['level_shift_added']
            or sum(row['shape'][0] for row in solver['calls'])!=3*count
            or any(row['solver_info']!=0 or row['precision']!='float64' or row['level_shift_added'] or not 0<=row['true_residual_norm']<=row['solver_absolute_target'] for row in solver['calls'])):
        raise ValueError('Full original-equation charged response required')
    nuclear=read_json(case/'nuclear.json');symmetry=read_json(case/'symmetry.json');log_path=case/'symmetry_thermo.log';used.append(log_path)
    match=re.search(r'rotational number\s+(\d+)',log_path.read_text())
    if match is None or int(match[1])!=symmetry['symmetry_number']:raise ValueError('Exact-geometry symmetry evidence differs')
    modes=load(ROOT/'scripts/validation/replay_archived_polymer_native.py','charged_native_independent_modes')
    equations=load(ROOT/'scripts/validation/audit/native_partition.py','charged_native_independent_nuclear')
    from redox.qm import qrrho
    masses=np.asarray(nuclear['atomic_masses_amu']);frequencies,rank=modes.projected_modes(matrix,atoms.positions,masses)
    if (not np.array_equal(masses,qrrho.atomic_masses(atoms.numbers)) or rank!=6 or len(frequencies)!=3*count-6 or min(frequencies)<-5
            or max(abs(frequencies-np.asarray(nuclear['frequencies_cm'])))>1e-5):raise ValueError('Independent native mode/minimum reconstruction failed')
    thermal={str(k):float(equations.nuclear_equations(frequencies,atoms.positions,masses,symmetry['symmetry_number'],k)) for k in (25,50,100)}
    if max(abs(v-nuclear['nuclear_by_cutoff_eV'][k]) for k,v in thermal.items())>1e-8:raise ValueError('Independent native nuclear equations differ')
    component=read_json(case/'thermodynamic_component.json');totals=component_sums(component,nuclear,center['E_Ha'])
    if (component['model']!=plan['solvent_model'] or component['geometry_sha256']!=sha256_file(geometry_path)
            or any(component[k]!=plan[v] for k,v in [('charge','charge'),('mult','mult'),('candidate_alias','exact_candidate_alias'),('quantum_settings','numerical_settings')])
            or component['graph_key']!=plan['prepared_case']['key']):raise ValueError('Native component identity/model/geometry binding differs')
    verify(namespace)
    return dict(status='charged_polymer_native_component_independently_reconstructed',namespace=str(namespace),
        candidate_alias=plan['exact_candidate_alias'],graph_key=plan['prepared_case']['key'],charge=plan['charge'],mult=1,
        model=plan['solvent_model'],quantum_settings=plan['numerical_settings'],geometry_sha256=sha256_file(geometry_path),
        atoms=count,maximum_column_error=max(errors),minimum_frequency_cm=float(min(frequencies)),nuclear_by_cutoff_eV=thermal,
        E_geometry_level_Ha=center['E_Ha'],G_by_cutoff_eV=totals,native_G_qualified=True,primary_model_admission=False,
        conformer_ensemble_qualified=False,material_admission=False,peak_GPU_MiB=stage['observed_peak_GPU_MiB'],
        source_sha256={str(p):sha256_file(p) for p in used},producer_sha256=sha256_file(__file__))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--spec',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out.exists():raise FileExistsError('Preserve prior native reviews')
    result=review(a.spec.resolve());atomic_write_json(a.out,result);print(json.dumps({k:result[k] for k in ('status','graph_key','maximum_column_error','minimum_frequency_cm')}))
