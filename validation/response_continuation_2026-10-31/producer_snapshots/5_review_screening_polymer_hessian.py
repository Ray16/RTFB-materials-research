#!/usr/bin/env python3
"""Independently review a complete corrected polymer Hessian; assemble no G."""
import argparse
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path

import h5py
import numpy as np
from ase.io import read as read_atoms

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def column_checks(matrix,center,columns,checks,read_json):
    """Reconstruct displaced coordinates and every declared derivative norm."""
    coordinates=np.asarray(center['coordinates_Bohr'],dtype=float)
    required={(atom,axis,step) for atom,axis in columns for step in (.0005,.001)}
    observed=set();errors=[]
    for check in checks:
        key=(check['atom'],check['axis'],check['step_Bohr'])
        if key not in required or key in observed:raise ValueError('Duplicate or undeclared column/step')
        observed.add(key);atom,axis,step=key;vectors=[]
        for sign,label in [(-1,'minus'),(1,'plus')]:
            point=read_json(Path(check[label]));expected=coordinates.copy();expected[atom,axis]+=sign*step
            gradient=np.asarray(point['gradient_Ha_Bohr'],dtype=float)
            if (not point['SCF_converged'] or gradient.shape!=coordinates.shape or not np.isfinite(gradient).all()
                    or not np.allclose(point['coordinates_Bohr'],expected,atol=1e-12,rtol=0)
                    or abs(point['spin'].get('spin_contam',math.inf))>.1):
                raise ValueError('Actual displaced converged state or coordinates differ')
            vectors.append(gradient)
        finite=(vectors[1]-vectors[0]).ravel()/(2*step)
        error=float(np.linalg.norm(finite-matrix[:,3*atom+axis])/max(np.linalg.norm(finite),1e-14))
        if not math.isfinite(error) or error>.01 or abs(error-check['relative_error'])>1e-10:
            raise ValueError('Reconstructed original derivative gate failed')
        errors.append(error)
    if observed!=required:raise ValueError('Complete original two-step column panel required')
    return errors


def review(namespace,placement):
    spec=verify(namespace);case=Path(spec['case']);used=[namespace]
    def read(path):
        used.append(path);return json.loads(path.read_text())
    plan=read(case/'plan.json');stage=read(case/'status.json')
    owner=read(case/'residency_guard/owner.json');residency=read(case/'residency_guard/state.json')
    host=read(case/'host_guard/status.json');runtime=read(case/'tensor_runtime/status.json')
    launch=read(placement/'launch.json');completion=read(placement/'completion.json')
    owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','full_hessian_closed_owners')
    if (plan['role']!='screening_polymer_Hessian_control' or plan['screening_variant']!='corrected'
            or stage['status']!='screening_polymer_Hessian_control_finished_review_pending'
            or stage['native_G_qualified'] or owners.alive(stage['owner']) is not False
            or completion['returncode']!=0 or completion['budget_stop'] is not None
            or not owner['guard_stopped'] or residency['foreign_resident_PIDs'] or not residency['own_residency_observed']
            or host['status']!='effective_memory_guard_finished_profile_review_pending'
            or runtime['status']!='task_returned_under_checked_cutensor_runtime'
            or host['maximum_RSS_GiB']!=80 or host['minimum_free_GiB']!=32 or host['free_RAM_bytes']<32*2**30
            or max(host['peak_tree_RSS_bytes'],host['self_VmHWM_bytes'])>80*2**30
            or not 0<stage['observed_peak_GPU_MiB']<=31744
            or not stage['owner']['pid']==owner['pid']==host['pid']==runtime['pid']
            or not stage['owner']['process_start']==owner['process_start']==host['process_start']
            or not owner['CPU_affinity']==host['CPU_affinity']==runtime['CPU_affinity']==launch['CPU_affinity']
            or len(owner['CPU_affinity'])!=2 or runtime['program']!=str(Path(spec['sources'])/'worker.py')):
        raise ValueError('Closed exact full-Hessian control and normal resource chain required')
    reference=plan['checked_native_reference'];atoms=read_atoms(reference['geometry']);count=len(atoms)
    used.extend(Path(reference[k]) for k in ('geometry','Hessian','partition'))
    matrix_path=case/'Hessian.npy';matrix=np.load(matrix_path,allow_pickle=False);used.append(matrix_path)
    center=read(case/'points/0000/evaluation.json');coordinates=np.asarray(center['coordinates_Bohr'])
    gradient=np.asarray(center['gradient_Ha_Bohr'])
    hetero=next((i for i,s in enumerate(plan['symbols']) if s=='O'),next((i for i,s in enumerate(plan['symbols']) if s=='Fe'),count-1))
    required={(0,1),(hetero,2)}
    required.update(tuple(map(int,np.unravel_index(i,gradient.shape))) for i in np.argsort(abs(gradient).ravel())[-2:])
    required.add((plan['failed_column']['atom'],plan['failed_column']['axis']))
    columns=read(case/'columns.json');derivatives=read(case/'derivatives.json')
    if (matrix.dtype!=np.float64 or matrix.shape!=(3*count,3*count) or not np.isfinite(matrix).all()
            or abs(matrix-matrix.T).max()>1e-6 or gradient.shape!=(count,3) or not np.isfinite(gradient).all()
            or abs(gradient).max()>.00045 or not center['SCF_converged']
            or not np.allclose(coordinates,atoms.positions/.529177210903,atol=1e-10,rtol=0)
            or {tuple(v) for v in columns['columns']}!=required or not derivatives['complete']):
        raise ValueError('Retained geometry, stationarity, full matrix or column coverage fails')
    errors=column_checks(matrix,center,required,derivatives['checks'],read)
    if len(list((case/'points').glob('*/evaluation.json')))!=1+4*len(required):raise ValueError('SCF point inventory differs')
    for path in sorted((case/'points').glob('*/evaluation.json')):
        point=read(path)
        if point.get('SCF_state_restored'):
            if path!=case/'points/0000/evaluation.json' or not plan.get('resume_components'):
                raise ValueError('Only the declared center may be restored')
            frame=point['source_frame']
            if frame!=plan['resume_components']['frame']:raise ValueError('Restored orbital frame differs')
            if any(sha256_file(p)!=h for p,h in frame['source_sha256'].items()):raise ValueError('Restored source changed')
            energy=read(case/'restart_energy.json')
            if (not energy['original_grid_arrays_available'] or abs(energy['difference_Ha'])>1e-8
                    or abs(energy['recreated_energy_Ha']-point['E_Ha'])>1e-8):raise ValueError('Exact-grid source energy does not reproduce')
            used.extend(Path(p) for p in frame['source_sha256'])
        else:
            checkpoint=Path(point['SCF_checkpoint']);settings=read(checkpoint.with_name('settings.json'))
            solved=read(checkpoint.with_name('SCF_result.json'));used.append(checkpoint)
            if (any(settings[k]!=v for k,v in dict(charge=0,mult=1,basis='def2-svp',xc='r2scan',disp='d4',conv_tol=1e-12,conv_tol_grad=1e-8,max_cycle=200,grid_level=3).items())
                    or not np.allclose(settings['coordinates_Bohr'],point['coordinates_Bohr'],atol=1e-12,rtol=0)
                    or not solved['converged'] or abs(solved['energy_Ha']-point['E_Ha'])>1e-10):
                raise ValueError('Actual SCF settings, geometry or result differ')
            with h5py.File(checkpoint,'r') as handle:
                if abs(float(handle['scf/e_tot'][()])-point['E_Ha'])>1e-10:raise ValueError('Actual SCF checkpoint differs')
    if plan.get('resume_components'):
        reuse=read(case/'component_reuse/reuse.json');capture=read(case/'solution_H1/identity.json')
        if (not reuse['PCM_H1_recomputed'] or reuse['stage_calls']!={'partial_hess_elec':1,'make_h1':1}
                or capture['scope']!='full SMD H1 including PCM contribution'
                or sha256_file(capture['path'])!=capture['file_sha256']):
            raise ValueError('Complete solvent H1 must survive component reuse')
        used.append(Path(capture['path']))
    solver=read(case/'execution/response_solver.json')
    if (solver['status']!='response_equations_solved_full_physical_QC_separate' or solver['equations_changed'] or solver['level_shift_added']
            or sum(row['shape'][0] for row in solver['calls'])!=3*count
            or any(row['solver_info']!=0 or row['precision']!='float64' or row['level_shift_added']
                   or not 0<=row['true_residual_norm']<=row['solver_absolute_target'] for row in solver['calls'])):
        raise ValueError('Full source-bound original-equation response coverage required')
    response_reused=0
    if plan.get('resume_response'):
        continuation=read(case/'response_reuse/status.json');capture=read(case/'solution_H1/identity.json')
        if (capture.get('equilibrium_solvation') is not True or continuation['equations_changed']
                or continuation['original_solver_sha256']!=plan['resume_response']['original_solver_source_sha256']
                or len(continuation['calls'])!=len(solver['calls'])):
            raise ValueError('A source-bound equilibrium-response continuation is required')
        for index,(candidate,solved) in enumerate(zip(continuation['calls'],solver['calls'])):
            if candidate['index']!=index or candidate['reused']!=solved.get('checkpoint_reused',False):
                raise ValueError('Response continuation decisions differ from the actual solver records')
            if candidate['reused']:
                if (not 0<=candidate['fresh_true_residual_norm']<=candidate['current_absolute_target']
                        or candidate['fresh_true_residual_norm']!=solved['true_residual_norm']
                        or candidate['current_absolute_target']!=solved['solver_absolute_target']):
                    raise ValueError('Reused block lacks its fresh unchanged-equation residual check')
                response_reused+=1
        for path,digest in plan['resume_response']['source_sha256'].items():
            if sha256_file(path)!=digest:raise ValueError('Response continuation backing changed')
            used.append(Path(path))
    nuclear=read(case/'nuclear.json');parent=read(Path(reference['partition']))
    masses=np.asarray(nuclear['atomic_masses_amu'])
    modes=load(ROOT/'scripts/validation/replay_archived_polymer_native.py','independent_full_hessian_modes')
    equations=load(ROOT/'scripts/validation/audit/native_partition.py','independent_full_hessian_nuclear')
    if masses.shape!=(count,) or not np.array_equal(masses,parent['nuclear']['atomic_masses_amu']):raise ValueError('Exact geometric atom masses differ')
    frequencies,rank=modes.projected_modes(matrix,atoms.positions,masses)
    frequency_error=float(abs(frequencies-np.asarray(nuclear['frequencies_cm'])).max())
    if rank!=6 or len(frequencies)!=3*count-6 or min(frequencies)<-5 or frequency_error>1e-5:
        raise ValueError('Independent internal projection/minimum check failed')
    thermal={str(k):float(equations.nuclear_equations(frequencies,atoms.positions,masses,parent['symmetry_number'],k)) for k in (25,50,100)}
    thermal_errors={k:abs(v-nuclear['nuclear_by_cutoff_eV'][k]) for k,v in thermal.items()}
    if max(thermal_errors.values())>1e-8:raise ValueError('Independent nuclear equations do not reconstruct')
    original=np.load(reference['Hessian'],allow_pickle=False)
    verify(namespace)
    return dict(status='corrected_full_polymer_Hessian_independently_reconstructed',namespace=str(namespace),
        atoms=count,columns=[list(v) for v in sorted(required)],maximum_column_error=max(errors),minimum_frequency_cm=float(min(frequencies)),
        frequency_error_cm=frequency_error,nuclear_errors_eV=thermal_errors,nuclear_by_cutoff_eV=thermal,
        maximum_Hessian_difference_from_reference_Ha_Bohr2=float(abs(matrix-original).max()),
        nuclear_difference_from_reference_eV={k:v-parent['nuclear']['nuclear_by_cutoff_eV'][k] for k,v in thermal.items()},
        exact_component_reuse_exercised=bool(plan.get('resume_components')),equilibrium_solvent_H1_retained=True,
        response_checkpoint_blocks_reused_after_fresh_operator_check=response_reused,
        residuals_independently_reevaluated=False,native_G_assembled=False,native_G_qualified=False,complete_polymer_ensemble_qualified=False,
        source_sha256={str(p):sha256_file(p) for p in used},producer_sha256=sha256_file(__file__))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('spec','placement','out'):parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if args.out.exists():raise FileExistsError('Preserve earlier physical reviews')
    result=review(args.spec.resolve(),args.placement.resolve());atomic_write_json(args.out,result)
    print(json.dumps({k:v for k,v in result.items() if k in ('status','maximum_column_error','minimum_frequency_cm','nuclear_difference_from_reference_eV','exact_component_reuse_exercised')}))
