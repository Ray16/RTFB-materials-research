#!/usr/bin/env python
"""Reconstruct the full fixed-orbital derivative comparison in bounded RAM."""
import argparse
import json
import math
from pathlib import Path
import time

from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def compare(plan_path,out):
    import numpy as np
    specification=json.loads(plan_path.read_text())
    roots=[Path(specification[name]) for name in ('reference','test')]
    arrays=[];sources={str(plan_path):api.digest(plan_path),str(Path(__file__).resolve()):api.digest(__file__)};identities=[]
    for root in roots:
        plan=api.read_plan(root);status_path=root/'status.json';status=json.loads(status_path.read_text())
        if status['status']!='isolated_DF_derivative_resource_probe_complete_comparison_pending':
            raise ValueError('Both isolated derivative probes must finish within their resource limits')
        api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256']);api.verify_sources(status['output_sha256'])
        state_path=Path(plan['raw_root'])/'fixed_state.json';state=json.loads(state_path.read_text());identities.append(state)
        array_path=Path(plan['raw_root'])/'DF_derivative_MO_occ.npy';arrays.append(np.load(array_path,mmap_mode='r'))
        for p in (root/'plan.json',status_path,state_path,array_path):sources[str(p)]=api.digest(p)
    if (identities[0]['source_checkpoint']!=identities[1]['source_checkpoint']
            or identities[0]['coordinates_Bohr']!=identities[1]['coordinates_Bohr']
            or identities[0]['basis']!=identities[1]['basis'] or arrays[0].shape!=arrays[1].shape
            or identities[0]['SCF_performed'] or identities[1]['SCF_performed']):
        raise ValueError('The full derivative arrays use different fixed states')
    if any(array.dtype!=np.float64 for array in arrays):
        raise ValueError('Full-precision FP64 derivative arrays are required')
    maximum=0.;norm=0.;residual=0.
    for i in range(arrays[0].shape[0]):
        reference=np.asarray(arrays[0][i]);test=np.asarray(arrays[1][i])
        if not np.isfinite(reference).all() or not np.isfinite(test).all():raise ValueError('Nonfinite derivative array')
        delta=test-reference;maximum=max(maximum,float(np.max(abs(delta))));norm+=float(np.sum(reference*reference));residual+=float(np.sum(delta*delta))
    relative=math.sqrt(residual/max(norm,1e-30))
    if maximum>specification['maximum_absolute_difference'] or relative>specification['maximum_relative_Frobenius_difference']:
        raise ValueError('Cross-runtime fixed-state derivative comparison exceeds its predeclared bounds')
    result=dict(status='full_fixed_orbital_DF_derivative_arrays_match',checked=time.time(),shape=list(arrays[0].shape),
        maximum_absolute_difference=maximum,relative_Frobenius_difference=relative,source_sha256=sources,
        SCF_recomputed=False,full_native_Hessian_qualified=False,G_qualified=False)
    atomic_write_json(out,result);print(json.dumps({k:v for k,v in result.items() if k!='source_sha256'}))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--plan',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--wait-hours',type=float,default=0)
    args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    import os
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    if not math.isfinite(args.wait_hours) or not 0<=args.wait_hours<=48:
        parser.error('Use a bounded wait up to 48 hours')
    deadline=time.time()+args.wait_hours*3600
    spec=json.loads(args.plan.read_text())
    while args.wait_hours:
        states=[]
        for role in ('reference','test'):
            p=Path(spec[role])/'status.json';states.append(json.loads(p.read_text()).get('status','') if p.exists() else '')
        if any('failed' in state for state in states):raise ValueError('A derivative probe failed; comparison not admitted')
        if all(state=='isolated_DF_derivative_resource_probe_complete_comparison_pending' for state in states):break
        if time.time()>=deadline:raise TimeoutError('Derivative comparison wait ended; partial records preserved')
        time.sleep(15)
    compare(args.plan.resolve(),args.out.resolve())


if __name__=='__main__':main()
