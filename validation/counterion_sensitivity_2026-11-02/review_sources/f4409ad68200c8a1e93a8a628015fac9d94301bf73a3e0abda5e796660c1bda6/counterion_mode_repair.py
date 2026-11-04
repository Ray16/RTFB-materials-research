"""Diagnose a closed failed minimum and prepare both signed-mode relaxations.

The original Hessian, SCFs and failed attempt remain untouched. No frequency is
made positive and no G is assembled. Each branch needs fresh optimization,
electronic energy, native Hessian and independent physical review.
"""
import argparse
import copy
import json
import os
from pathlib import Path
import shutil

import numpy as np
from ase.io import read, write
from pyscf import gto
from pyscf.hessian import thermo

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation import primary_checkpoints
from redox.validation.execution_review import freeze_source
from redox.workflows.molecular_native_namespace import verify, load

SOURCE_BYTES=Path(__file__).read_bytes()


def diagnose(namespace, output):
    if output.exists():raise FileExistsError('Preserve previous failed-minimum diagnoses')
    spec=verify(namespace);case=Path(spec['case']);sources=Path(spec['sources']);used=[namespace]
    def read_json(path):used.append(Path(path));return json.loads(Path(path).read_text())
    plan=read_json(case/'plan.json');stage=read_json(case/'status.json')
    if (plan['role']!='primary_counterion_thermochemistry' or (plan['charge'],plan['mult'])!=(1,1)
            or stage['status']!='counterion_thermochemistry_failed_preserved'
            or stage.get('error')!="RuntimeError('Source DFT geometry is not a minimum of its own solution surface')"):
        raise ValueError('An actual diagnosed native minimum failure is required')
    owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','counterion_failed_mode_owners')
    placements=list(case.parent.glob('placement*/launch.json'))
    if len(placements)!=1:raise ValueError('One actual scientific placement required')
    launch=read_json(placements[0]);completion=read_json(placements[0].with_name('completion.json'))
    guard=read_json(case/'residency_guard/owner.json');residency=read_json(case/'residency_guard/state.json')
    host=read_json(case/'host_guard/status.json')
    identities=[stage['owner'],dict(host=launch['host'],pid=launch['pid'],process_start=launch['process_start']),launch['budget_observer_owner']]
    if (any(owners.alive(owner) is not False for owner in identities) or completion['returncode'] in (None,0)
            or completion['budget_stop'] is not None or not guard['guard_stopped']
            or residency['foreign_resident_PIDs'] or not residency['own_residency_observed']
            or (case/'residency_guard/external_interruption.json').exists()
            or max(host['peak_tree_RSS_bytes'],host['self_VmHWM_bytes'])>40*2**30
            or host['free_RAM_bytes']<32*2**30):
        raise ValueError('Closed physical failure without resource or residency interference required')
    geometry=case/'state/opt.xyz';atoms=read(geometry);used.append(geometry)
    paths=list((case/'state/native_thermal').glob('*/hessian_cartesian.npy'))
    if len(paths)!=1:raise ValueError('One complete preserved native Hessian required')
    matrix=np.load(paths[0],allow_pickle=False);used.append(paths[0]);n=len(atoms)
    if matrix.shape!=(3*n,3*n) or not np.isfinite(matrix).all() or abs(matrix-matrix.T).max()>1e-6:
        raise ValueError('Complete finite symmetric Cartesian Hessian required')
    series=primary_checkpoints.checkpoint_series(case/'SCF_evaluations',geometry,
        dict(charge=1,mult=1,basis='def2-svp',xc='r2scan',disp='d4',conv_tol=1e-12,
             conv_tol_grad=1e-8,max_cycle=200,grid_level=3,solvent_class='SMD'),read_json,used,expected_count=9)
    evaluations=[]
    for p in sorted((case/'SCF_evaluations').glob('*/*/settings.json')):
        value=read_json(p);gpath=p.with_name('gradient.npy');used.append(gpath)
        evaluations.append((np.asarray(value['coordinates_Bohr']),np.load(gpath,allow_pickle=False)))
    center=evaluations[0][0]
    def gradient_at(target):
        found=[gradient for coords,gradient in evaluations if np.max(abs(coords-target))<1e-10]
        if len(found)!=1:raise ValueError('Each declared mode-audit gradient must have one exact SCF frame')
        return found[0]
    checks=[]
    for atom,axis in [(0,1),(n-1,2)]:
        for step in (.0005,.001):
            plus=center.copy();minus=center.copy();plus[atom,axis]+=step;minus[atom,axis]-=step
            numerical=(gradient_at(plus)-gradient_at(minus)).ravel()/(2*step)
            error=float(np.linalg.norm(numerical-matrix[3*atom+axis])/np.linalg.norm(numerical))
            if not np.isfinite(error) or error>.01:raise ValueError('Saved finite-difference column does not qualify the failed Hessian')
            checks.append(dict(atom=atom,axis=axis,step_Bohr=step,relative_error=error))
    qrrho=load(Path(spec['archive'])/'src/redox/qm/qrrho.py','failed_counterion_qualified_mass_source')
    mass=qrrho.atomic_masses(atoms.numbers)
    modes=load(ROOT/'scripts/validation/replay_archived_polymer_native.py','failed_counterion_independent_modes')
    frequencies,rank=modes.projected_modes(matrix,atoms.positions,mass)
    mol=gto.M(atom=list(zip(atoms.get_chemical_symbols(),atoms.positions.tolist())),basis='def2-svp',charge=1,spin=0,verbose=0)
    harmonic=thermo.harmonic_analysis(mol,matrix.reshape(n,3,n,3).transpose(0,2,1,3),mass=mass,imaginary_freq=False)
    if rank!=6 or min(frequencies)>=-5 or np.max(abs(np.asarray(harmonic['freq_wavenumber'])-frequencies))>1e-5:
        raise ValueError('The independent projected unstable mode must agree with the molecular analysis')
    index=int(np.argmin(frequencies));mode=np.asarray(harmonic['norm_mode'][index]);mode/=np.linalg.norm(mode,axis=1).max()
    mode_path=output.with_suffix('.mode.npy');np.save(mode_path,mode,allow_pickle=False);used.append(mode_path)
    used.extend([freeze_source(output.parent,Path(__file__),SOURCE_BYTES),
        freeze_source(output.parent,Path(primary_checkpoints.__file__),primary_checkpoints.SOURCE_BYTES),
        freeze_source(output.parent,Path(modes.__file__)),Path(qrrho.__file__)])
    result=dict(status='counterion_unstable_native_mode_independently_reconstructed',namespace=str(namespace),
        geometry=str(geometry),geometry_sha256=sha256_file(geometry),Hessian=str(paths[0]),Hessian_sha256=sha256_file(paths[0]),
        minimum_frequency_cm=float(min(frequencies)),frequencies_cm=frequencies.tolist(),rigid_rank=rank,
        negative_internal_modes_below_gate=int(sum(frequencies<-5)),actual_native_SCF_frames_checked=series['actual_SCF_frames_checked'],
        sampled_derivative_checks=checks,maximum_column_error=max(r['relative_error'] for r in checks),
        signed_mode=str(mode_path),signed_mode_sha256=sha256_file(mode_path),mode_normalization='maximum Cartesian atom displacement norm = 1',
        proposed_maximum_atom_displacement_A=.1,native_G_qualified=False,association_free_energy_computed=False,
        source_sha256={str(p):sha256_file(p) for p in used})
    atomic_write_json(output,result);return result


def prepare(diagnosis,output):
    from redox.validation.protonation import hydrogen_topology_status
    report=json.loads(diagnosis.read_text())
    if (report['status']!='counterion_unstable_native_mode_independently_reconstructed'
            or report['native_G_qualified'] or any(sha256_file(p)!=h for p,h in report['source_sha256'].items())):
        raise ValueError('The complete independently reconstructed failed mode is required')
    failed=verify(Path(report['namespace']));failed_sources=Path(failed['sources'])
    geometry_review=json.loads((failed_sources/'geometry_review.json').read_text());original=verify(Path(geometry_review['namespace']))
    old_plan=json.loads((Path(original['case'])/'plan.json').read_text())
    if output.exists() or any(output.is_relative_to(Path(original[k])) for k in ('prior_root','archive','case','sources')):
        raise FileExistsError('Fresh repair branches outside every scientific source required')
    atoms=read(report['geometry']);mode=np.load(report['signed_mode'],allow_pickle=False)
    output.mkdir(parents=True);rows=[]
    for sign in (-1,1):
        root=output/('minus' if sign<0 else 'plus');case=root/'case';sources=root/'sources';case.mkdir(parents=True)
        shutil.copytree(original['sources'],sources,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        shutil.copyfile(ROOT/'scripts/validation/run_counterion_geometry.py',sources/'worker.py')
        shutil.copyfile(diagnosis,sources/'failed_mode_diagnosis.json');shutil.copyfile(report['signed_mode'],sources/'failed_mode.npy')
        (sources/'mode_repair_preparer.py').write_bytes(SOURCE_BYTES)
        seed=atoms.copy();seed.positions+=sign*.1*mode;write(case/'input.xyz',seed,format='xyz')
        if hydrogen_topology_status(old_plan['smiles'],case/'input.xyz')!='ok':raise ValueError('The signed-mode seed changed the exact graph')
        plan=copy.deepcopy(old_plan)
        plan.update(seed=dict(source='verified negative native mode',sign=sign,maximum_atom_displacement_A=.1,
            source_geometry_sha256=report['geometry_sha256'],source_Hessian_sha256=report['Hessian_sha256'],
            source_minimum_frequency_cm=report['minimum_frequency_cm']),
            optimizer_convergence=dict(convergence_energy=1e-9,convergence_grms=5e-6,convergence_gmax=1e-5,
                convergence_drms=1e-4,convergence_dmax=2e-4),
            mode_repair_of=report['namespace'],original_geometry_qualification_transferred=False,
            source_sha256={**old_plan['source_sha256'],**{str(p):sha256_file(p) for p in [*sources.glob('*'),case/'input.xyz'] if p.is_file()}})
        atomic_write_json(case/'plan.json',plan);spec=dict(original,case=str(case),sources=str(sources),plan_sha256=sha256_file(case/'plan.json'))
        atomic_write_json(root/'namespace.json',spec);verify(root/'namespace.json')
        rows.append(dict(key='N0_above_mode_'+('minus' if sign<0 else 'plus'),namespace=str(root/'namespace.json'),namespace_sha256=sha256_file(root/'namespace.json'),sign=sign))
    result=dict(status='two_signed_counterion_mode_repair_geometries_prepared',states=rows,source_diagnosis=str(diagnosis),
        source_diagnosis_sha256=sha256_file(diagnosis),native_G_qualified=False,quantum_jobs_submitted=0,
        selection='Both branches retained; fresh geometry and native minimum checks required before any thermodynamic selection')
    atomic_write_json(output/'manifest.json',result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('diagnose','prepare'))
    for key in ('spec','diagnosis','out'):p.add_argument('--'+key,type=Path)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out is None:p.error('A fresh output path is required')
    if a.action=='diagnose':
        if a.spec is None:p.error('Diagnosis requires the closed failed namespace')
        result=diagnose(a.spec.resolve(),a.out.resolve());print(json.dumps({k:result[k] for k in ('status','minimum_frequency_cm','maximum_column_error')}))
    else:
        if a.diagnosis is None:p.error('Preparation requires the verified failed-mode diagnosis')
        result=prepare(a.diagnosis.resolve(),a.out.resolve());print(json.dumps(dict(status=result['status'],branches=len(result['states']))))
