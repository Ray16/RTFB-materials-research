#!/usr/bin/env python3
"""Compare a preserved derivative failure with its fully checked native repeat.

Acceptance by the declared gate and diagnosis of a numerical discrepancy are
separate conclusions. The original failed state never gains thermochemistry.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.namespace_sources import NamespaceSources,read_checked_review
from redox.workflows.archived_native_namespace import verify_spec


def compare(original_path,repeat_path):
    original=verify_spec(original_path);repeat=verify_spec(repeat_path);checked,binding=read_checked_review(repeat_path,'native')
    cases=[Path(s['attempt']) for s in (original,repeat)];plans=[json.loads((p/'plan.json').read_text()) for p in cases]
    fields=('tasks','quantum_settings','experimental_model','experimental_method_id','experimental_native_settings','scored_geometry_sha256','maximum_host_RSS_GiB','maximum_peak_GPU_MiB','execution_driver')
    if any(plans[0][key]!=plans[1][key] for key in fields):raise ValueError('Repeat changed a chemical, numerical or geometry field')
    stage=json.loads((cases[0]/'experimental_native/status.json').read_text())
    if stage['status']!='experimental_native_failed_preserved' or stage.get('error')!="ValueError('Alternate native sampled derivative gate failed')":raise ValueError('A preserved sampled-derivative failure is required')
    readers=[NamespaceSources(s) for s in (original,repeat)]
    for reader,plan in zip(readers,plans):reader.verify_all(plan['producer_sha256']);reader.verify_all(plan['input_sha256'])
    roots=[Path(s['raw_parent'])/s['key'] for s in (original,repeat)];used=[original_path,repeat_path,*[p/'plan.json' for p in cases],cases[0]/'experimental_native/status.json',Path(binding['review'])]
    matrices=[np.load(p/'hessian.npy',allow_pickle=False) for p in roots]
    if matrices[0].shape!=matrices[1].shape or not all(np.isfinite(m).all() for m in matrices):raise ValueError('Finite matched Hessian dimensions required')
    points=[];checks=[]
    for root in roots:
        rows={}
        for path in sorted((root/'SCF_evaluations').glob('*/evaluation.json')):
            row=json.loads(path.read_text());coords=np.asarray(row['coordinates_Bohr'],dtype=np.float64);gradient=np.asarray(row['gradient_Ha_Bohr'])
            if row.get('SCF_converged') is not True or coords.shape!=gradient.shape or not np.isfinite([row['E_Ha']]).all() or not np.isfinite(gradient).all():raise ValueError('Finite converged matched SCF point required')
            key=hashlib.sha256(coords.tobytes()).hexdigest()
            if key in rows:raise ValueError('Repeated point in one native calculation')
            rows[key]=dict(row,path=str(path));used.append(path)
        points.append(rows);checks.append(json.loads((root/'derivative_checks.json').read_text()));used.extend([root/'hessian.npy',root/'gradient.npy',root/'derivative_checks.json'])
    common=sorted(set(points[0])&set(points[1]));differences=[]
    for key in common:
        a,b=points[0][key],points[1][key]
        differences.append(dict(coordinates_sha256=key,original=a['path'],repeat=b['path'],energy_difference_Ha=b['E_Ha']-a['E_Ha'],maximum_gradient_difference_Ha_Bohr=float(abs(np.asarray(b['gradient_Ha_Bohr'])-np.asarray(a['gradient_Ha_Bohr'])).max())))
    column_rows=[]
    for check in checks[1]['checks']:
        key=(check['atom'],check['axis'],check['step_Bohr']);old=next((c['relative_error'] for c in checks[0]['checks'] if (c['atom'],c['axis'],c['step_Bohr'])==key),None)
        column_rows.append(dict(atom=key[0],axis=key[1],step_Bohr=key[2],original_relative_error=old,repeat_relative_error=check['relative_error']))
    error_limit=plans[1]['experimental_native_settings']['FD_relative_limit']
    if not checks[1]['complete'] or max(c['relative_error'] for c in checks[1]['checks'])>error_limit:raise ValueError('The repeat does not pass its unchanged declared derivative gate')
    delta=matrices[1]-matrices[0];hashes={str(p):sha256_file(p) for p in used}
    for reader in readers:reader.recheck()
    if any(sha256_file(p)!=h for p,h in hashes.items()):raise ValueError('Comparison source changed')
    return dict(status='preserved_native_derivative_failure_and_checked_repeat_compared',scientific_fields_unchanged=list(fields),matched_SCF_points=len(common),SCF_point_differences=differences,sampled_columns=column_rows,
        maximum_Hessian_difference_Ha_Bohr2=float(abs(delta).max()),relative_Hessian_difference=float(np.linalg.norm(delta)/np.linalg.norm(matrices[0])),
        unchanged_declared_relative_error_limit=error_limit,repeat_maximum_relative_error=max(c['relative_error'] for c in checks[1]['checks']),
        repeat_passes_declared_native_gate=True,repeat_review=binding,repeat_G_by_cutoff_eV=checked['G_by_cutoff_eV'],
        original_failed_state_qualified=False,gradient_reproducibility_cause_identified=False,finite_difference_step_convergence_established=False,
        source_sha256=hashes,producer_sha256=sha256_file(__file__))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('original','repeat','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out.exists():raise FileExistsError('Preserve previous comparisons')
    r=compare(a.original.resolve(),a.repeat.resolve());atomic_write_json(a.out,r)
    print(json.dumps({k:r[k] for k in ('status','matched_SCF_points','maximum_Hessian_difference_Ha_Bohr2','repeat_maximum_relative_error')}))
