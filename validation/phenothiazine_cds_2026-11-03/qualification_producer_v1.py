"""Qualify a targeted CDS numerical step for fresh full native controls.

This admits no G: exact saved neutral failures identify the component and all
fresh DFT/Hessian/derivative/minimum/nuclear checks remain mandatory.
"""
import argparse
import json
import os
from pathlib import Path

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json


def check_components(paths):
    import numpy as np
    expected={'phenothiazine_neu_mode_minus','phenothiazine_neu_mode_plus'}
    reports=[json.loads(p.read_text()) for p in paths]
    if len(reports)!=2 or {r['key'] for r in reports}!=expected:
        raise ValueError('Both exact neutral branch component controls required')
    sources={str(p):sha256_file(p) for p in paths};summaries=[]
    for path,r in zip(paths,reports):
        if (r['status']!='legacy_CDS_Hessian_step_component_study_not_full_QC'
                or r['native_G_qualified'] or r['CUDA_initialized']
                or r['CDS_energy_and_gradient_functions_changed'] or r['electrostatic_model_changed']
                or any(sha256_file(p)!=h for p,h in r['source_sha256'].items())):
            raise ValueError('Unchanged-model CPU component backing required')
        arrays=path.with_suffix('.npz')
        if str(arrays) not in r['source_sha256']:raise ValueError('Exact component arrays required')
        with np.load(arrays,allow_pickle=False) as data:
            baseline=data['saved_full'];residual=data['residual_without_default_CDS']
            if (baseline.shape!=(117,117) or not np.isfinite(baseline).all()
                    or not np.isfinite(residual).all()
                    or not np.allclose(baseline-data['CDS_0.0001'],residual,rtol=0,atol=1e-14)):
                raise ValueError('Saved full/default-CDS decomposition differs')
            if (abs(baseline-baseline.T).max()<=1e-6 or abs(residual-residual.T).max()>1e-6
                    or abs(data['CDS_1e-06']-data['CDS_1e-07']).max()>1e-8):
                raise ValueError('Default failure, residual symmetry and fine-step convergence required')
            candidate=residual+data['CDS_1e-06']
            if abs(candidate-candidate.T).max()>1e-6:
                raise ValueError('Finer component does not resolve the existing symmetry threshold')
        if (len(r['energy_gradient_checks'])!=2 or any(c['tested_components']!=117
                or c['maximum_energy_gradient_error_Ha_Bohr']>1e-6 for c in r['energy_gradient_checks'])):
            raise ValueError('Every unchanged CDS gradient component must agree with energy differences')
        sources.update(r['source_sha256'])
        summaries.append(dict(key=r['key'],geometry_sha256=r['geometry_sha256'],component_report=str(path),
            fine_pair_maximum_change_Ha_Bohr2=r['fine_pair_maximum_change_Ha_Bohr2'],
            residual_maximum_asymmetry_Ha_Bohr2=r['residual_without_default_CDS_maximum_asymmetry_Ha_Bohr2']))
    return dict(status='two_exact_neutral_CDS_step_controls_qualified_for_fresh_native_tests',
        step_Bohr=1e-6,library_default_step_Bohr=1e-4,controls=summaries,
        full_native_gates_qualified=False,native_G_qualified=False,redox_potential_qualified=False,
        scientific_model_changed=False,source_sha256=sources,producer_sha256=sha256_file(__file__))


def checked_qualification(path):
    row=json.loads(path.read_text());paths=[Path(c['component_report']) for c in row['controls']]
    repeated=check_components(paths)
    if row!=repeated:raise ValueError('CDS component qualification or producer changed')
    return row


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--component-review',type=Path,action='append',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out.exists():raise FileExistsError('Preserve prior component qualifications')
    result=check_components([p.resolve() for p in a.component_review]);atomic_write_json(a.out,result)
    print(json.dumps({k:result[k] for k in ('status','step_Bohr','native_G_qualified')}))
