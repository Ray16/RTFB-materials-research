#!/usr/bin/env python
"""Reconstruct alternate-model geometry/force gates without primary admission."""
import argparse
import hashlib
import json
from pathlib import Path
import time

import numpy as np

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def review(campaign,out):
    model=api.load_module('smooth_geometry_review_model',ROOT/'scripts/pipeline/smooth_cds_control.py')
    controls=api.load_module('smooth_geometry_review_controls',ROOT/'scripts/validation/review_smooth_model_panel.py')
    owner=api.load_module('smooth_geometry_review_owner',ROOT/'scripts/pipeline/run_polymer_spacing.py')
    plan=api.read_plan(campaign);api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
    path=campaign/'geometry/status.json';row=json.loads(path.read_text())
    if (row['status']!='experimental_smooth_geometry_checked_energy_and_native_QC_pending'
            or owner.owner_alive(row['owner']) or row['primary_candidate_admission']):
        raise ValueError('A finished alternate-model geometry with no primary admission is required')
    api.verify_sources(row['output_sha256'])
    proof=controls.checked_controls(Path(plan['full_quantum_controls']))
    sd=Path(row['state_dir'])
    if (sd/'result.json').exists():raise ValueError('Alternate model must not create a legacy primary geometry record')
    result_path=sd/'experimental_geometry_result.json';result=json.loads(result_path.read_text())
    if (not result['geometry_optimization_converged'] or not result['optimization_SCF_converged']
            or result['primary_candidate_admission'] or result['experimental_model']!=model.MODEL):
        raise ValueError('Alternate geometry optimizer/model metadata failed')
    model_path=campaign/'model/model.json';actual=json.loads(model_path.read_text())
    if (actual['status']!='experimental_smooth_CDS_context_finished_QC_separate'
            or any(actual['profile'][k]!=v for k,v in model.MODEL.items())
            or actual['energy_calls']<1 or actual['gradient_calls']<1):
        raise ValueError('Declared alternate energy and analytic gradient must be exercised')
    api.verify_sources(actual['profile']['implementation_sha256'])
    portable=hashlib.sha256(json.dumps(model.MODEL,sort_keys=True).encode()).hexdigest()
    if portable!=plan['experimental_method_id']:raise ValueError('Alternate method identity changed')
    forces=sd/'total_energy_force_check';baseline=json.loads((forces/'baseline.json').read_text())
    coords=np.array(baseline['coordinates_Bohr']);gradient=np.array(baseline['gradient_Ha_Bohr'])
    report=json.loads((forces/'summary.json').read_text());errors=[];tested=set()
    for item in report['checks']:
        atom,axis,step=item['atom'],item['axis'],item['step_Bohr'];tested.add((atom,axis))
        values=[]
        for sign in (-1,1):
            p=forces/f'atom{atom}_axis{axis}_h{step}_sign{sign}.json';point=json.loads(p.read_text())
            expected=coords.copy();expected[atom,axis]+=sign*step
            if not point['SCF_converged'] or not np.allclose(point['coordinates_Bohr'],expected,rtol=0,atol=1e-12):
                raise ValueError('Force finite-difference geometry/SCF binding differs')
            values.append(point['E_Ha']);proof[str(p)]=api.digest(p)
        numeric=(values[1]-values[0])/(2*step);error=abs(numeric-gradient[atom,axis])
        if (not np.isfinite(error) or error>1e-6 or abs(error-item['absolute_error'])>1e-12
                or step not in (.0005,.001)):
            raise ValueError('Independent alternate-model energy/force reconstruction failed')
        errors.append(float(error))
    largest={tuple(map(int,np.unravel_index(i,gradient.shape))) for i in np.argsort(abs(gradient).ravel())[-2:]}
    required=largest|{tuple(v) for v in plan['extra_force_axes']}
    if not required.issubset(tested) or len(report['checks'])!=2*len(tested):
        raise ValueError('Largest and declared critical force coordinates must all be checked at both steps')
    if not np.isfinite(gradient).all() or abs(gradient).max()>.00045:
        raise ValueError('Independent stationarity threshold failed')
    audit_path=ROOT/'scripts/validation/audit_polymer_geometries.py'
    audit=api.load_module('smooth_final_exact_graph_audit',audit_path).audit(plan['tasks'][0],sd/'opt.xyz')
    if audit['issues']:raise ValueError('Final graph/stereochemistry failed')
    rt_path=campaign/'tensor_runtime/status.json';rt=json.loads(rt_path.read_text());profile=Path(plan['runtime_profile'])
    expected=json.loads(profile.read_text())
    if (rt['status']!='task_returned_under_checked_cutensor_runtime' or rt['pid']!=row['owner']['pid']
            or rt['profile_sha256']!=api.digest(profile) or rt['versions']!=expected['versions']
            or rt['library_sha256']!=expected['library_sha256'] or rt['total_GPU_bytes']<31*2**30):
        raise ValueError('Actual completed 32-GB-class runtime differs')
    guard=campaign/'residency_guard';identity=json.loads((guard/'owner.json').read_text());residency=json.loads((guard/'state.json').read_text())
    host_path=campaign/'host_guard/status.json';host=json.loads(host_path.read_text())
    if (not identity['guard_stopped'] or residency['foreign_resident_PIDs'] or not residency['own_residency_observed']
            or (guard/'interruption.json').exists() or host['status']!='program_returned_review_required'
            or not 0<host['peak_RSS_bytes']<=32*2**30 or not 0<row['observed_peak_GPU_MiB']<=31744):
        raise ValueError('Pilot runtime ownership or memory validation failed')
    for p in (Path(__file__).resolve(),campaign/'plan.json',path,result_path,model_path,rt_path,host_path,
              forces/'baseline.json',forces/'summary.json',guard/'owner.json',guard/'state.json',sd/'opt.xyz',audit_path):
        proof[str(p)]=api.digest(p)
    proof.update(plan['producer_sha256']);proof.update(plan['input_sha256']);proof.update(row['output_sha256'])
    result=dict(status='alternate_smooth_geometry_and_resource_pilot_independently_checked',checked=time.time(),
        campaign=str(campaign),state_dir=str(sd),method=model.MODEL,portable_method_id=portable,
        gradient_max_Ha_Bohr=float(abs(gradient).max()),maximum_force_error_Ha_Bohr=max(errors),
        tested_force_axes=[list(v) for v in sorted(tested)],steps_Bohr=[.0005,.001],
        observed_peak_GPU_MiB=row['observed_peak_GPU_MiB'],peak_host_RSS_bytes=host['peak_RSS_bytes'],
        graph_review=audit,source_sha256=proof,primary_candidate_admission=False,
        active_electronic_energy_qualified=False,native_Hessian_or_G_qualified=False,
        scope='New-model stationary geometry and measured resource point only; source legacy failure preserved, no primary record or free energy substituted')
    out.parent.mkdir(parents=True,exist_ok=True);atomic_write_json(out,result)
    print(json.dumps({k:result[k] for k in ('status','gradient_max_Ha_Bohr','maximum_force_error_Ha_Bohr','observed_peak_GPU_MiB')}))


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--campaign',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True);args=parser.parse_args();review(args.campaign.resolve(),args.out.resolve())


if __name__=='__main__':main()
