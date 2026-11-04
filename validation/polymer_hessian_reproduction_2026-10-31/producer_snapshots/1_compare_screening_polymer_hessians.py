#!/usr/bin/env python3
"""Compare separately reviewed fresh/restarted full molecular Hessian controls.

Every underlying control must pass its ordinary physical review first. This
comparison measures numerical reproduction, not chemical accuracy or G.
"""
import argparse
from itertools import combinations
import json
import math
import os
from pathlib import Path

import numpy as np

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify


def pair(left,right,limits):
    identity=('original_namespace','original_namespace_sha256','numerical_settings','symbols')
    if any(left['plan'][key]!=right['plan'][key] for key in identity):
        raise ValueError('One retained geometry, chemical/model identity and protocol required')
    a,b=np.asarray(left['center']['coordinates_Bohr']),np.asarray(right['center']['coordinates_Bohr'])
    n=len(left['plan']['symbols'])
    if (a.shape!=(n,3) or b.shape!=(n,3) or not np.isfinite(a).all() or not np.isfinite(b).all()
            or not np.allclose(a,b,atol=1e-12,rtol=0)):
        raise ValueError('The retained Cartesian geometries differ')
    matrices=[np.asarray(row['Hessian']) for row in (left,right)]
    gradients=[np.asarray(row['center']['gradient_Ha_Bohr']) for row in (left,right)]
    frequencies=[np.asarray(row['nuclear']['frequencies_cm']) for row in (left,right)]
    for matrix,gradient,frequency,row in zip(matrices,gradients,frequencies,(left,right)):
        if (matrix.shape!=(3*n,3*n) or matrix.dtype!=np.float64 or not np.isfinite(matrix).all()
                or gradient.shape!=(n,3) or not np.isfinite(gradient).all()
                or frequency.shape!=(3*n-6,) or not np.isfinite(frequency).all()
                or not math.isfinite(row['center']['E_Ha'])):
            raise ValueError('Complete finite matched Hessian/force/mode/energy arrays required')
    difference=matrices[0]-matrices[1];denominator=max(np.linalg.norm(m) for m in matrices)
    if not denominator>0:raise ValueError('Nonzero Hessian norms required')
    thermal={key:abs(left['nuclear']['nuclear_by_cutoff_eV'][key]-right['nuclear']['nuclear_by_cutoff_eV'][key]) for key in ('25','50','100')}
    if not all(math.isfinite(v) for v in thermal.values()):raise ValueError('Finite nuclear terms required')
    values=dict(maximum_center_energy_difference_Ha=abs(left['center']['E_Ha']-right['center']['E_Ha']),
        maximum_center_gradient_difference_Ha_Bohr=float(abs(gradients[0]-gradients[1]).max()),
        maximum_absolute_Hessian_difference_Ha_Bohr2=float(abs(difference).max()),
        maximum_relative_Frobenius_Hessian_difference=float(np.linalg.norm(difference)/denominator),
        maximum_internal_frequency_difference_cm=float(abs(frequencies[0]-frequencies[1]).max()),
        maximum_nuclear_difference_eV=max(thermal.values()))
    if any(not math.isfinite(limits[key]) or limits[key]<=0 for key in values):raise ValueError('Positive finite declared tolerances required')
    passed={key:value<=limits[key] for key,value in values.items()}
    return dict(left=left['namespace'],right=right['namespace'],differences=values,nuclear_differences_eV=thermal,
        passed=passed,all_within_declared_tolerances=all(passed.values()))


def compare(reviews,config):
    hashes={};rows=[]
    def bind(path,digest=None):
        path=Path(path);actual=hashes.get(str(path))
        if actual is None:actual=sha256_file(path);hashes[str(path)]=actual
        if digest is not None and actual!=digest:raise ValueError('Changed independent-review backing: '+str(path))
    bind(config);bind(Path(__file__));limits=json.loads(config.read_text())
    for path in reviews:
        bind(path);record=json.loads(path.read_text())
        if (record['status']!='corrected_full_polymer_Hessian_independently_reconstructed'
                or record['native_G_qualified']):raise ValueError('Completed independently reviewed corrected controls required')
        for source,digest in record['source_sha256'].items():bind(source,digest)
        namespace=Path(record['namespace']);spec=verify(namespace);case=Path(spec['case'])
        payload={name:json.loads((case/filename).read_text()) for name,filename in [('plan','plan.json'),('center','points/0000/evaluation.json'),('nuclear','nuclear.json')]}
        for filename in ['plan.json','points/0000/evaluation.json','nuclear.json','Hessian.npy']:
            expected=record['source_sha256'].get(str(case/filename))
            if expected is None:raise ValueError('Independent review must bind every numerical comparison input')
            bind(case/filename,expected)
        if payload['plan']['screening_variant']!='corrected':raise ValueError('No mixed screening implementation comparison')
        rows.append(dict(namespace=str(namespace),case=str(case.resolve()),Hessian=np.load(case/'Hessian.npy',allow_pickle=False),**payload))
    if len(rows)<2 or len({row['namespace'] for row in rows})!=len(rows) or len({row['case'] for row in rows})!=len(rows):raise ValueError('Distinct independently reviewed calculations required')
    pairs=[pair(a,b,limits) for a,b in combinations(rows,2)]
    if any(sha256_file(path)!=digest for path,digest in hashes.items()):raise ValueError('Comparison backing changed during evaluation')
    return dict(status='checked_corrected_full_Hessian_controls_compared',controls=len(rows),comparisons=pairs,
        all_pairs_within_declared_tolerances=all(row['all_within_declared_tolerances'] for row in pairs),
        declared_tolerances=limits,native_G_qualified=False,chemical_accuracy_established=False,
        full_polymer_ensemble_qualified=False,source_sha256=hashes,producer_sha256=sha256_file(__file__))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--review',type=Path,action='append',required=True)
    p.add_argument('--config',type=Path,default=ROOT/'config/polymer_hessian_reproduction.json');p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out.exists():raise FileExistsError('Preserve previous comparison records')
    result=compare([path.resolve() for path in a.review],a.config.resolve());atomic_write_json(a.out,result)
    print(json.dumps({k:result[k] for k in ('status','controls','all_pairs_within_declared_tolerances')}))
    if not result['all_pairs_within_declared_tolerances']:raise SystemExit(2)
