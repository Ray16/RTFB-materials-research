#!/usr/bin/env python
"""Independently reconstruct alternate-model native state records."""
import argparse
import json
import math
from pathlib import Path
import re
import time

import numpy as np

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def review(campaign,out):
    plan=api.read_plan(campaign);api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    path=campaign/'experimental_native/status.json';row=json.loads(path.read_text())
    owners=api.load_module('alternate_native_review_stopped_owner',ROOT/'scripts/pipeline/run_polymer_spacing.py')
    if row['status']!='experimental_native_finished_independent_review_pending' or owners.owner_alive(row['owner']):
        raise ValueError('A completed native stage with stopped ownership is required')
    api.verify_sources(row['output_sha256']);raw=Path(row['raw_root']);xyz=raw/'geometry.xyz'
    if api.digest(xyz)!=plan['scored_geometry_sha256']:raise ValueError('Native and electronic geometry differ')
    native=json.loads((raw/'experimental_partition.json').read_text());parent=json.loads(Path(plan['electronic_review']).read_text())
    api.verify_sources(parent['source_sha256'])
    driver_path=ROOT/'scripts/pipeline/smooth_polymer_native.py'
    driver=api.load_module('alternate_native_review_scope_gates',driver_path)
    driver.validate_partition(native,parent['E_smd_eV'],plan['experimental_method_id'])
    if native['native_settings']!=driver.NATIVE_SETTINGS or plan['experimental_native_settings']!=driver.NATIVE_SETTINGS:
        raise ValueError('Actual native numerical settings differ')
    hessian=np.load(raw/'hessian.npy',allow_pickle=False);gradient=np.load(raw/'gradient.npy',allow_pickle=False)
    from ase.io import read
    atoms=read(xyz);count=len(atoms);symbols=atoms.get_chemical_symbols();coords=atoms.positions/.529177210903
    if (hessian.shape!=(3*count,3*count) or gradient.shape!=(count,3)
            or not np.isfinite(hessian).all() or not np.isfinite(gradient).all()
            or abs(hessian-hessian.T).max()>1e-6 or abs(gradient).max()>.00045):
        raise ValueError('Native Hessian/gradient shape, symmetry or stationarity failed')
    hetero=next((i for i,s in enumerate(symbols) if s=='O'),next((i for i,s in enumerate(symbols) if s=='Fe'),count-1))
    required={(0,1),(hetero,2)}
    required.update(tuple(map(int,np.unravel_index(i,gradient.shape))) for i in np.argsort(abs(gradient).ravel())[-2:])
    if {tuple(v) for v in native['required_columns']}!=required:
        raise ValueError('Native protocol and largest-force column coverage differs')
    baseline=json.loads((raw/'SCF_evaluations/0000/evaluation.json').read_text())
    if (not np.allclose(baseline['coordinates_Bohr'],coords,rtol=0,atol=1e-10)
            or not np.allclose(baseline['gradient_Ha_Bohr'],gradient,rtol=0,atol=1e-13)):
        raise ValueError('Baseline derivative/geometry binding differs')
    files=sorted((raw/'SCF_evaluations').glob('*/evaluation.json'))
    if len(files)!=1+4*len(required):raise ValueError('Raw displaced-gradient coverage is incomplete')
    index=1;checks=[]
    for atom,axis in sorted(required):
        for step in (.0005,.001):
            points=[json.loads((raw/f'SCF_evaluations/{i:04d}/evaluation.json').read_text()) for i in (index,index+1)];index+=2
            vectors=[]
            for sign,point in zip((-1,1),points):
                expected=coords.copy();expected[atom,axis]+=sign*step
                if (not np.allclose(point['coordinates_Bohr'],expected,rtol=0,atol=1e-12)
                        or not point['SCF_converged'] or abs(point['spin'].get('spin_contam',math.inf))>.1
                        or point['experimental_method_id']!=plan['experimental_method_id']):
                    raise ValueError('Displaced SCF/gradient/model binding failed')
                vectors.append(np.asarray(point['gradient_Ha_Bohr']))
            finite=(vectors[1]-vectors[0]).ravel()/(2*step);analytic=hessian[:,3*atom+axis]
            error=float(np.linalg.norm(finite-analytic)/max(np.linalg.norm(finite),1e-14))
            if not math.isfinite(error) or error>.01:raise ValueError('Independent native Hessian-column test failed')
            stored=next(v for v in native['derivative_checks'] if (v['atom'],v['axis'],v['step_Bohr'])==(atom,axis,step))
            if abs(error-stored['relative_error'])>1e-10:raise ValueError('Stored derivative error does not reconstruct')
            checks.append(error)
    controls=api.load_module('alternate_native_independent_mode_projection',ROOT/'scripts/validation/review_smooth_quantum.py')
    equations=api.load_module('alternate_native_independent_nuclear_equations',ROOT/'scripts/validation/audit/native_partition.py')
    nuclear=native['nuclear'];frequencies,rank=controls.projected_modes(hessian,atoms.positions,np.array(nuclear['atomic_masses_amu']))
    error=float(np.max(abs(frequencies-np.array(nuclear['frequencies_cm']))))
    if rank!=6 or len(frequencies)!=3*count-6 or error>1e-5 or min(frequencies)<-5:
        raise ValueError('Independent internal modes/minimum check failed')
    sigma=int(re.search(r'rotational number\s+(\d+)',(raw/'symmetry_thermo.log').read_text())[1])
    if sigma<1 or sigma!=native['symmetry_number']:raise ValueError('Geometry symmetry number differs')
    errors={str(k):abs(equations.nuclear_equations(frequencies,atoms.positions,np.array(nuclear['atomic_masses_amu']),sigma,k)
        -nuclear['nuclear_by_cutoff_eV'][str(k)]) for k in (25,50,100)}
    if not all(math.isfinite(v) and v<1e-8 for v in errors.values()):raise ValueError('Independent native nuclear sums differ')
    model_path=campaign/'model/model.json';model=json.loads(model_path.read_text())
    if (model['status']!='experimental_smooth_CDS_context_finished_QC_separate'
            or model['energy_calls']<1 or model['gradient_calls']<1
            or any(model['profile'][k]!=v for k,v in plan['experimental_model'].items())):
        raise ValueError('Actual solvent-surface model differs')
    runtime_path=campaign/'tensor_runtime/status.json';rt=json.loads(runtime_path.read_text());profile_path=Path(plan['runtime_profile']);profile=json.loads(profile_path.read_text())
    if (rt['status']!='task_returned_under_checked_cutensor_runtime' or rt['pid']!=row['owner']['pid']
            or rt['program']!=plan['execution_driver'] or rt['profile_sha256']!=api.digest(profile_path)
            or rt['versions']!=profile['versions'] or rt['library_sha256']!=profile['library_sha256']):
        raise ValueError('Actual native contraction runtime differs')
    guard=campaign/'residency_guard';identity=json.loads((guard/'owner.json').read_text());residency=json.loads((guard/'state.json').read_text())
    host_path=campaign/'host_guard/status.json';host=json.loads(host_path.read_text())
    if (not identity['guard_stopped'] or residency['foreign_resident_PIDs'] or not residency['own_residency_observed']
            or (guard/'interruption.json').exists() or host['status']!='program_returned_review_required'
            or not 0<host['peak_RSS_bytes']<=80*2**30 or not 0<row['observed_peak_GPU_MiB']<=31744):
        raise ValueError('Native ownership/residency or resource gate failed')
    proof={**plan['producer_sha256'],**plan['input_sha256'],**row['output_sha256']}
    for p in (Path(__file__).resolve(),campaign/'plan.json',path,model_path,runtime_path,guard/'owner.json',guard/'state.json',host_path,raw/'symmetry_thermo.log'):
        proof[str(p)]=api.digest(p)
    result=dict(status='alternate_native_state_independently_reconstructed',checked=time.time(),campaign=str(campaign),
        task=plan['tasks'][0],experimental_model=plan['experimental_model'],experimental_method_id=plan['experimental_method_id'],
        G_eV=native['G_by_cutoff_eV']['50'],G_by_cutoff_eV=native['G_by_cutoff_eV'],
        minimum_frequency_cm=float(min(frequencies)),frequency_reconstruction_error_cm=error,
        maximum_sampled_column_error=max(checks),tested_Hessian_columns=len(required),
        SCF_gradient_evaluations=len(files),nuclear_reconstruction_errors_eV=errors,
        observed_peak_GPU_MiB=row['observed_peak_GPU_MiB'],peak_host_RSS_bytes=host['peak_RSS_bytes'],
        source_sha256=proof,primary_candidate_admission=False,conformer_convergence_established=False,
        scope='Numerically and physically checked native state of the distinct smooth model; no legacy cache substitution or material eligibility')
    out.parent.mkdir(parents=True,exist_ok=True);atomic_write_json(out,result)
    print(json.dumps({k:result[k] for k in ('status','G_eV','minimum_frequency_cm','tested_Hessian_columns')}))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--campaign',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    review(a.campaign.resolve(),a.out.resolve())


if __name__=='__main__':main()
