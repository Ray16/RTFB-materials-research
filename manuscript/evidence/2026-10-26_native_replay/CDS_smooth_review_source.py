#!/usr/bin/env python
"""Reconstruct full smooth-CDS control checks from raw evaluations and Hessians."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def checked_controls(root):
    path=Path(root)/'independent_review.json';row=json.loads(path.read_text())
    if row['status']!='experimental_smooth_quantum_controls_independently_checked' or row['primary_candidate_admission']:
        raise ValueError('Independent two-spin quantum model controls required')
    api.verify_sources(row['source_sha256'])
    return {str(path):api.digest(path),**row['source_sha256']}


def projected_modes(matrix,xyz,mass):
    centered=xyz-np.average(xyz,axis=0,weights=mass)
    vectors=[(np.sqrt(mass)[:,None]*axis).ravel() for axis in np.eye(3)]
    vectors += [(np.sqrt(mass)[:,None]*np.cross(centered,axis)).ravel() for axis in np.eye(3)]
    rigid=[]
    for vector in vectors:
        residual=vector.copy()
        for _ in range(2):
            for unit in rigid:residual-=unit*np.dot(unit,residual)
        if np.linalg.norm(residual)>1e-10*max(1.,np.linalg.norm(vector)):
            rigid.append(residual/np.linalg.norm(residual))
    basis=np.array(rigid).T
    values,vectors=np.linalg.eigh(np.eye(matrix.shape[0])-basis@basis.T)
    internal=vectors[:,values>.5]
    roots=np.repeat(np.sqrt(mass*(1.660539040e-27/9.10938356e-31)),3)
    weighted=(matrix+matrix.T)/(2*np.outer(roots,roots))
    roots=np.linalg.eigvalsh(internal.T@weighted@internal)
    return np.sign(roots)*np.sqrt(abs(roots))*219474.63067,len(rigid)


def review(root):
    from ase.io import read
    model=api.load_module('smooth_control_qualification_model',ROOT/'scripts/pipeline/smooth_cds_control.py')
    equations=api.load_module('smooth_control_independent_nuclear',ROOT/'scripts/validation/audit/native_partition.py')
    owner=api.load_module('smooth_control_stopped_owners',ROOT/'scripts/pipeline/run_polymer_spacing.py')
    profile_path=ROOT/'environments/cutensor22/runtime.json';runtime_profile=json.loads(profile_path.read_text())
    proof=model.component_evidence(ROOT/'results/validation/smooth_CDS_2026-10-25/pilot')
    proof.update({str(p):api.digest(p) for p in (Path(__file__).resolve(),ROOT/'scripts/pipeline/smooth_cds_control.py',
        ROOT/'scripts/validation/audit/native_partition.py',ROOT/'src/redox/qm/qrrho.py',
        ROOT/'src/redox/qm/xtb_kernel.py',profile_path)})
    controls=[]
    for mult in (1,2):
        path=root/'controls'/str(mult)/'status.json';record=json.loads(path.read_text())
        if record['status']!='two_model_quantum_force_Hessian_and_nuclear_controls_passed' or owner.owner_alive(record['owner']):
            raise ValueError('Both full quantum controls must finish normally')
        api.verify_sources(record['source_sha256']);api.verify_sources(record['input_sha256']);api.verify_sources(record['output_sha256'])
        rt_path=root/'runtime'/str(mult)/'status.json';rt=json.loads(rt_path.read_text())
        if (rt['status']!='task_returned_under_checked_cutensor_runtime' or rt['pid']!=record['owner']['pid']
                or rt['versions']!=runtime_profile['versions'] or rt['profile_sha256']!=api.digest(profile_path)
                or rt['source_sha256']!=api.digest(ROOT/'scripts/fleet/tensor_runtime_worker.py')):
            raise ValueError('Exact completed contraction runtime must match the quantum control')
        guard=root/'guards'/str(mult)
        residency=json.loads((guard/'state.json').read_text());identity=json.loads((guard/'owner.json').read_text())
        if not identity['guard_stopped'] or not residency['own_residency_observed'] or residency['foreign_resident_PIDs'] or (guard/'interruption.json').exists():
            raise ValueError('Quantum control did not retain exclusive residency')
        if not 0<record['observed_peak_GPU_MiB']<=14336:raise ValueError('Control resource ceiling failed')
        results={}
        for arm in ('legacy','smooth'):
            directory=Path(record['raw_root'])/arm
            baseline=json.loads((directory/'SCF_evaluations/0000/evaluation.json').read_text())
            xyz=np.array(baseline['coordinates_Bohr']);gradient=np.load(directory/'gradient.npy')
            matrix=np.load(directory/'hessian.npy');force_errors=[];column_errors=[];index=1
            if matrix.shape!=(12,12) or abs(matrix-matrix.T).max()>1e-6 or not np.isfinite(matrix).all():
                raise ValueError('Invalid four-atom control Hessian')
            for atom in range(4):
                for axis in range(3):
                    for step in (.0005,.001):
                        a,b=[json.loads((directory/f'SCF_evaluations/{i:04d}/evaluation.json').read_text()) for i in (index,index+1)];index+=2
                        expected=np.zeros_like(xyz);expected[atom,axis]=step
                        if (not np.allclose(np.array(a['coordinates_Bohr'])-xyz,expected,rtol=0,atol=1e-12)
                                or not np.allclose(np.array(b['coordinates_Bohr'])-xyz,-expected,rtol=0,atol=1e-12)):
                            raise ValueError('Finite-difference SCF geometry binding differs')
                        if any(not row['SCF_converged'] or abs(row['spin'].get('spin_contam',float('inf')))>.1 for row in (a,b)):
                            raise ValueError('Perturbed SCF/spin check failed')
                        force_errors.append(abs((a['E_Ha']-b['E_Ha'])/(2*step)-gradient[atom,axis]))
                        if (atom,axis) in {(0,1),(1,2)}:
                            finite=(np.array(a['gradient_Ha_Bohr'])-np.array(b['gradient_Ha_Bohr'])).ravel()/(2*step)
                            col=matrix[:,3*atom+axis]
                            column_errors.append(float(np.linalg.norm(finite-col)/max(np.linalg.norm(col),1e-14)))
            if max(force_errors)>1e-6 or max(column_errors)>.01 or abs(gradient).max()>.00045:
                raise ValueError('Independent force/Hessian/stationarity reconstruction failed')
            nuclear=json.loads((directory/'nuclear.json').read_text());atoms=read(directory/'opt.xyz')
            frequencies,rank=projected_modes(matrix,atoms.positions,np.array(nuclear['atomic_masses_amu']))
            error=float(np.max(abs(frequencies-np.array(nuclear['frequencies_cm']))))
            if rank!=6 or error>1e-5 or min(frequencies)<-5:raise ValueError('Independent projected-mode check failed')
            corrections={str(k):equations.nuclear_equations(frequencies,atoms.positions,
                np.array(nuclear['atomic_masses_amu']),1,k) for k in (25,50,100)}
            if max(abs(corrections[k]-nuclear['nuclear_by_cutoff_eV'][k]) for k in corrections)>1e-8:
                raise ValueError('Independent nuclear equations differ')
            results[arm]=dict(maximum_force_error_Ha_Bohr=float(max(force_errors)),
                maximum_Hessian_column_relative_error=max(column_errors),tested_force_components=12,
                tested_Hessian_columns=2,steps_Bohr=[.0005,.001],mode_reconstruction_error_cm=error,
                minimum_internal_frequency_cm=float(min(frequencies)),nuclear_by_cutoff_eV=corrections)
            if arm=='smooth':
                state=json.loads((directory/'model/model.json').read_text())
                if state['status']!='experimental_smooth_CDS_context_finished_QC_separate' or any(state['profile'][k]!=v for k,v in model.MODEL.items()):
                    raise ValueError('Actual smooth model differs from its declared definition')
        for p in (path,rt_path,guard/'owner.json',guard/'state.json'):proof[str(p)]=api.digest(p)
        proof.update(record['source_sha256']);proof.update(record['input_sha256']);proof.update(record['output_sha256'])
        controls.append(dict(mult=mult,arms=results,peak_GPU_MiB=record['observed_peak_GPU_MiB']))
    method_id=hashlib.sha256(json.dumps(model.MODEL,sort_keys=True).encode()).hexdigest()
    result=dict(status='experimental_smooth_quantum_controls_independently_checked',checked=time.time(),
        method=model.MODEL,portable_method_id=method_id,controls=controls,source_sha256=proof,
        primary_candidate_admission=False,large_geometry_or_native_G_qualified=False,
        scope='Two four-atom full quantum method controls; imposed comparison symmetry number one. All large states and their model-specific thermochemistry remain separate.')
    atomic_write_json(root/'independent_review.json',result)
    print(json.dumps(dict(status=result['status'],controls=len(controls),method_id=method_id)))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--root',type=Path,required=True)
    args=parser.parse_args();review(args.root.resolve())


if __name__=='__main__':main()
