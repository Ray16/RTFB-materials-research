#!/usr/bin/env python3
"""Independently reconstruct a closed primary molecular native control."""
import argparse
import json
import math
import os
import re
import hashlib

from pathlib import Path

REVIEWER_SOURCE_BYTES=Path(__file__).read_bytes()

import numpy as np
from ase.io import read as read_atoms

from redox.core.common import ROOT,HARTREE_EV
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify,load
from redox.validation.execution_review import review_closure,freeze_source,SOURCE_BYTES
from redox.validation.scf_frame import freeze_reviewer,check_spin_contamination
from redox.workflows.ionic_stationary_cycles import screening_signature


def review(namespace):
    spec=verify(namespace);case=Path(spec['case']);sources=Path(spec['sources']);used=[namespace]
    def read(path):used.append(path);return json.loads(path.read_text())
    plan=read(case/'plan.json')
    if plan['role']!='primary_molecular_Hessian_control' or plan['geometry_repair_permitted']:
        raise ValueError('Fixed primary-model control required')
    owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','primary_native_closed_owner')
    stage=review_closure(case,read,owners.alive,
        expected_status='primary_molecular_Hessian_control_finished_review_pending',host_GiB=40,
        GPU_MiB=31744,worker=sources/'worker.py')
    implementation=case/'screening_implementation/implementation.json'
    compiled=implementation.with_name('screening_kernel.cu')
    signature=screening_signature(read(implementation),sources,compiled);used.append(compiled)
    used.extend([freeze_reviewer(case),freeze_source(case,Path(review_closure.__code__.co_filename),SOURCE_BYTES),freeze_source(case,Path(__file__),REVIEWER_SOURCE_BYTES)])
    result=read(case/'result.json');native_path=Path(result['native_record_path']);native=read(native_path)
    matrix_path=native_path.with_suffix('')/'hessian_cartesian.npy';used.append(matrix_path)
    matrix=np.load(matrix_path,allow_pickle=False);geometry=case/'state/opt.xyz';atoms=read_atoms(geometry);used.append(geometry)
    if (sha256_file(geometry)!=sha256_file(sources/'reference_geometry.xyz')
            or native['identity']['charge']!=plan['charge'] or native['identity']['mult']!=plan['mult']
            or native['status']!='derivative_checked_solution_thermal' or not native['stationary_minimum']
            or not native['grid_response'] or not 0<=native['gradient_max']<=.00045
            or matrix.shape!=(3*len(atoms),3*len(atoms)) or matrix.dtype!=np.float64
            or not np.isfinite(matrix).all() or abs(matrix-matrix.T).max()>1e-6):
        raise ValueError('Source geometry, native state, gradient or Hessian checks fail')
    check_spin_contamination(native['smd_spin_diagnostics'].get('spin_contam'))
    fd=load(ROOT/'scripts/validation/review_molecular_native.py','primary_native_independent_FD')
    used.append(freeze_source(case,Path(fd.__file__)))
    checks,paths=fd.finite_difference_review(case,geometry,native,matrix);used.extend(paths)
    from redox.validation import primary_checkpoints
    used.append(freeze_source(case,Path(primary_checkpoints.__file__),primary_checkpoints.SOURCE_BYTES))
    series=primary_checkpoints.checkpoint_series(case/'SCF_evaluations',geometry,
        dict(charge=plan['charge'],mult=plan['mult'],basis='def2-svp',xc='r2scan',disp='d4',
             conv_tol=1e-12,conv_tol_grad=1e-8,max_cycle=200,grid_level=3,solvent_class='SMD'),
        read,used,expected_count=9,native_diagnostics=native['smd_spin_diagnostics'])
    frame_errors=series['frame_errors'];spin_rows=series['spins']
    projections=[]
    for p in (case/'state/projected_thermal').glob('*.json'):
        x=read(p)
        if x.get('identity',{}).get('raw_native_sha256')==sha256_file(native_path):projections.append((p,x))
    if len(projections)!=1:raise ValueError('One source-bound projected partition required')
    projection_path,partition=projections[0]
    for artifact,digest in partition['artifacts_sha256'].items():
        if sha256_file(artifact)!=digest:raise ValueError('Projected partition artifact changed')
        used.append(Path(artifact))
    if (partition['status']!='projected_internal_qrrho_checked'
            or partition['identity']['geometry_sha256']!=sha256_file(geometry)
            or partition['identity']['Cartesian_Hessian_sha256']!=sha256_file(matrix_path)):
        raise ValueError('Projected native partition binding differs')
    for cutoff in (25,50,100):
        log=projection_path.with_suffix('')/f'thermo_CLI_{cutoff}.log';used.append(log)
        match=re.search(r'rotational number\s+(\d+)',log.read_text())
        if match is None or int(match[1])!=partition['symmetry_number']:
            raise ValueError('Exact-geometry CLI symmetry evidence differs')
    modes=load(ROOT/'scripts/validation/replay_archived_polymer_native.py','primary_native_modes')
    equations=load(ROOT/'scripts/validation/audit/native_partition.py','primary_native_nuclear_equations')
    used.extend([freeze_source(case,Path(modes.__file__)),freeze_source(case,Path(equations.__file__))])
    masses=np.asarray(partition['atomic_masses_amu'])
    reference=read(sources/'reference_partition.json')
    if not np.array_equal(masses,reference['atomic_masses_amu']):raise ValueError('Exact source mass convention differs')
    frequencies,rank=modes.projected_modes(matrix,atoms.positions,masses)
    if (rank!=6 or not np.isfinite(frequencies).all() or min(frequencies)<-5
            or np.max(abs(frequencies-np.asarray(partition['frequencies_cm'])))>1e-5):
        raise ValueError('Independent projected minimum check failed')
    nuclear={str(c):float(equations.nuclear_equations(frequencies,atoms.positions,masses,partition['symmetry_number'],c)) for c in (25,50,100)}
    spin=-8.617333262145e-5*298.15*math.log(plan['mult'])
    if abs(spin-partition['g_spin_partition_eV'])>1e-8:raise ValueError('Ground-multiplet spin term differs')
    energy=read(sources/'corrected_electronic.json');values={c:energy['e_smd_eV']+v+spin for c,v in nuclear.items()}
    if (any(abs(nuclear[c]-partition['nuclear_by_cutoff_eV'][c])>1e-8 or abs(values[c]-result['G_by_cutoff_eV'][c])>1e-8 for c in values)
            or abs(energy['e_smd_eV']-energy['e_smd_Ha']*HARTREE_EV)>1e-8):
        raise ValueError('Independent composite E/nuclear/spin sum differs')
    reference_matrix=np.load(sources/'reference_hessian.npy',allow_pickle=False);used.append(sources/'reference_hessian.npy')
    if reference_matrix.shape!=matrix.shape:raise ValueError('Reference matrix dimensions differ')
    verify(namespace)
    return dict(status='primary_molecular_native_control_independently_reconstructed',namespace=str(namespace),
        key=plan['source_id'],charge=plan['charge'],mult=plan['mult'],atoms=len(atoms),
        G_by_cutoff_eV=values,nuclear_by_cutoff_eV=nuclear,spin_partition_eV=spin,
        reference_G_difference_eV={c:values[c]-plan['reference_checked_G_by_cutoff_eV'][c] for c in values},
        reference_nuclear_difference_eV={c:nuclear[c]-reference['nuclear_by_cutoff_eV'][c] for c in values},
        maximum_Hessian_difference_Ha_Bohr2=float(abs(matrix-reference_matrix).max()),
        minimum_frequency_cm=float(min(frequencies)),maximum_column_error=max(r['relative_error'] for r in checks),
        actual_SCF_checkpoint_frames_checked=len(frame_errors),maximum_checkpoint_coordinate_error_Bohr=max(frame_errors),
        actual_SCF_wavefunction_spins_checked=len(spin_rows),
        maximum_reconstructed_spin_contamination=max(abs(r['spin_contamination']) for r in spin_rows),
        producer_spin_diagnostic_decimal_places=4,
        peak_GPU_MiB=stage['observed_peak_GPU_MiB'],native_component_requalified=True,
        numerical_implementation=signature,
        new_complete_free_energy_cycles_qualified=False,conformer_ensemble_qualified=False,material_admission=False,
        source_sha256={str(p):sha256_file(p) for p in used},producer_sha256=hashlib.sha256(REVIEWER_SOURCE_BYTES).hexdigest())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--spec',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out.exists():raise FileExistsError('Preserve earlier primary native reviews')
    r=review(a.spec.resolve());atomic_write_json(a.out,r);print(json.dumps({k:r[k] for k in ('status','key','minimum_frequency_cm','maximum_column_error')}))
