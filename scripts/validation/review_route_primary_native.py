#!/usr/bin/env python3
"""Independently reconstruct a first primary native component, including sulfur."""
import argparse
import json
import math
import os
from pathlib import Path

import numpy as np
from ase.io import read as read_atoms

from redox.core.common import ROOT,HARTREE_EV
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify,load
from redox.workflows.route_primary_native import validate_plan
from redox.validation.execution_review import review_closure,freeze_source
from redox.validation.primary_checkpoints import checkpoint_series
from redox.validation.protonation import hydrogen_topology_status
from redox.workflows.ionic_stationary_cycles import screening_signature


def review(namespace):
    spec=verify(namespace);case=Path(spec['case']);sources=Path(spec['sources']);used=[namespace]
    def read(p):used.append(Path(p));return json.loads(Path(p).read_text())
    plan=read(case/'plan.json');validate_plan(plan)
    owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','first_native_closed_owners')
    stage=review_closure(case,read,owners.alive,expected_status='route_primary_native_component_finished_review_pending',
        host_GiB=40,GPU_MiB=31744,worker=sources/'worker.py')
    used.extend([freeze_source(case,Path(__file__)),freeze_source(case,Path(validate_plan.__code__.co_filename)),
        freeze_source(case,Path(review_closure.__code__.co_filename)),freeze_source(case,Path(checkpoint_series.__code__.co_filename))])
    implementation=case/'screening_implementation/implementation.json';compiled=implementation.with_name('screening_kernel.cu')
    signature=screening_signature(read(implementation),sources,compiled);used.append(compiled)
    energy_review=read(sources/'electronic_review.json')
    if (energy_review['status']!='benchmark_screening_energy_control_independently_checked'
            or not energy_review['actual_checkpoint_frames_verified'] or not energy_review['declared_agreement_passed']
            or any(sha256_file(p)!=h for p,h in energy_review['source_sha256'].items())):
        raise ValueError('Exact independently reviewed electronic backing required')
    result=read(case/'result.json');native_path=Path(result['native_record_path'])
    if not native_path.resolve().is_relative_to((case/'state').resolve()):raise ValueError('Own native record required')
    native=read(native_path);matrix_path=native_path.with_suffix('')/'hessian_cartesian.npy';used.append(matrix_path)
    matrix=np.load(matrix_path,allow_pickle=False);geometry=case/'state/opt.xyz';used.append(geometry);atoms=read_atoms(geometry)
    if (sha256_file(geometry)!=plan['geometry_sha256'] or hydrogen_topology_status(plan['smiles'],geometry)!='ok'
            or (native['identity']['charge'],native['identity']['mult'])!=(plan['charge'],plan['mult'])
            or native['status']!='derivative_checked_solution_thermal' or not native['stationary_minimum']
            or not native['grid_response'] or not 0<=native['gradient_max']<=.00045
            or matrix.dtype!=np.float64 or matrix.shape!=(117,117) or not np.isfinite(matrix).all()
            or abs(matrix-matrix.T).max()>1e-6):raise ValueError('Own exact native geometry/state/gradient/matrix required')
    extra=read(case/'extra_derivatives.json')
    if extra!=result['extra_derivative_checks']:raise ValueError('Extra sulfur checks changed')
    fd=load(ROOT/'scripts/validation/review_molecular_native.py','first_native_independent_FD')
    used.append(freeze_source(case,Path(fd.__file__)))
    full=dict(native,derivative_checks=[*native['derivative_checks'],*extra])
    columns=[(0,1),(38,2),*map(tuple,plan['extra_columns'])]
    checks,paths=fd.finite_difference_review(case,geometry,full,matrix,expected_columns=columns);used.extend(paths)
    series=checkpoint_series(case/'SCF_evaluations',geometry,
        dict(charge=plan['charge'],mult=plan['mult'],basis='def2-svp',xc='r2scan',disp='d4',
            conv_tol=1e-12,conv_tol_grad=1e-8,max_cycle=200,grid_level=3,solvent_class='SMD'),
        read,used,expected_count=21,native_diagnostics=native['smd_spin_diagnostics'])
    equations=load(ROOT/'scripts/validation/audit/native_partition.py','first_native_independent_partition')
    used.append(freeze_source(case,Path(equations.__file__)))
    partition=equations.read_partition(native_path,geometry)
    if partition is None:raise ValueError('Independent masses, projected modes, minimum, symmetry and nuclear reconstruction failed')
    projection=read(Path(partition['path']))
    used.extend(Path(p) for p in partition['evidence_sha256'])
    energy=read(sources/'corrected_electronic.json');spin=-8.617333262145e-5*298.15*math.log(plan['mult'])
    if abs(spin-partition['g_spin_partition_eV'])>1e-8:raise ValueError('Ground-state spin partition differs')
    values={c:energy['e_smd_eV']+partition['nuclear_by_cutoff_eV'][c]+spin for c in ('25','50','100')}
    if (energy['energy_protocol']!=plan['electronic_protocol'] or not energy['converged_smd']
            or abs(energy['e_smd_Ha']*HARTREE_EV-energy['e_smd_eV'])>1e-8
            or any(not math.isfinite(v) or abs(v-result['G_by_cutoff_eV'][c])>1e-8
                or abs(partition['nuclear_by_cutoff_eV'][c]-result['nuclear_by_cutoff_eV'][c])>1e-8 for c,v in values.items())):
        raise ValueError('Own corrected electronic/nuclear/spin sum does not close')
    verify(namespace)
    return dict(status='route_primary_native_component_independently_reconstructed',namespace=str(namespace),
        candidate_alias='phenothiazine',key=plan['source_id'],smiles=plan['smiles'],charge=plan['charge'],mult=plan['mult'],
        atoms=39,native_nao=397,electronic_protocol=plan['electronic_protocol'],native_protocol=plan['native_protocol'],thermal_protocol=plan['thermal_protocol'],
        geometry_sha256=sha256_file(geometry),G_by_cutoff_eV=values,corrected_electronic_eV=energy['e_smd_eV'],
        nuclear_by_cutoff_eV=partition['nuclear_by_cutoff_eV'],spin_partition_eV=spin,
        minimum_frequency_cm=min(projection['frequencies_cm']),maximum_column_error=max(r['relative_error'] for r in checks),
        checked_derivative_columns=columns,actual_SCF_checkpoint_frames_checked=series['actual_SCF_frames_checked'],
        peak_GPU_MiB=stage['observed_peak_GPU_MiB'],numerical_implementation=signature,
        native_G_qualified=True,redox_potential_qualified=False,geometry_ensemble_qualified=False,material_admission=False,
        source_sha256={str(p):sha256_file(p) for p in used})


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--spec',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out.exists():raise FileExistsError('Preserve earlier independent native reviews')
    result=review(a.spec.resolve());atomic_write_json(a.out,result)
    print(json.dumps({k:result[k] for k in ('status','key','minimum_frequency_cm','maximum_column_error')}))
