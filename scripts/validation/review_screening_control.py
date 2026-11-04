#!/usr/bin/env python3
"""Independently check closed mask/integration controls; withhold SCF/Hessian use."""
import argparse
import importlib.util
import json
import math
import os
from pathlib import Path

import numpy as np

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows import molecular_native_namespace as execution


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def review(namespace):
    spec=execution.verify(namespace);case=Path(spec['case']);plan=json.loads((case/'plan.json').read_text())
    paths={name:case/name for name in ['status.json','host_guard/status.json','tensor_runtime/status.json','residency_guard/owner.json','residency_guard/state.json','independent_mask.npy','water_comparison.json']}
    data={name:json.loads(path.read_text()) for name,path in paths.items() if path.suffix=='.json'}
    stage=data['status.json'];host=data['host_guard/status.json'];runtime=data['tensor_runtime/status.json'];owner=data['residency_guard/owner.json'];residency=data['residency_guard/state.json']
    placement=case.parent/'placement';launch=json.loads((placement/'launch.json').read_text());completion=json.loads((placement/'completion.json').read_text())
    owners=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','screening_control_closed_owner')
    if (stage['status']!='screening_barrier_masks_and_requested_integration_controls_checked' or stage['native_G_qualified'] or stage['implementation_qualified']
            or completion['returncode']!=0 or completion['budget_stop'] is not None or owners.alive(stage['owner']) is not False
            or not owner.get('guard_stopped') or residency['foreign_resident_PIDs'] or not residency['own_residency_observed']
            or host['status']!='effective_memory_guard_finished_profile_review_pending' or host['maximum_RSS_GiB']!=plan['maximum_host_RSS_GiB']
            or max(host['peak_tree_RSS_bytes'],host['self_VmHWM_bytes'])>plan['maximum_host_RSS_GiB']*2**30 or host['free_RAM_bytes']<32*2**30
            or runtime['status']!='task_returned_under_checked_cutensor_runtime' or runtime['program']!=str(Path(spec['sources'])/'worker.py')
            or not stage['owner']['pid']==owner['pid']==host['pid']==runtime['pid']
            or not stage['owner']['process_start']==owner['process_start']==host['process_start']
            or not owner['CPU_affinity']==host['CPU_affinity']==runtime['CPU_affinity']==launch['CPU_affinity'] or len(host['CPU_affinity'])!=2):
        raise ValueError('Closed exact control and normal ownership/resource chain required')
    reference=np.load(paths['independent_mask.npy'],allow_pickle=False)
    if reference.dtype!=np.int8 or reference.shape!=(2,512):raise ValueError('Declared synthetic independent mask required')
    import hashlib
    digest=hashlib.sha256(reference.tobytes()).hexdigest();rows=stage['comparisons']
    if len(rows)!=60 or {(r['variant'],r['repeat']) for r in rows}!={(v,i) for v in ['installed','corrected'] for i in range(30)}:raise ValueError('Complete two-variant repeated mask comparison required')
    for row in rows:
        if type(row['mismatched_shell_blocks']) is not int or row['mismatched_shell_blocks']<0:raise ValueError('Invalid mismatch count')
        if row['variant']=='corrected' and (row['mismatched_shell_blocks'] or row['mask_sha256']!=digest):raise ValueError('Corrected mask does not match independent CPU data')
    water=json.loads((Path(spec['sources'])/'water.json').read_text());results=data['water_comparison.json']['rows'];limit=water['sorted_unsorted_absolute_limit']
    if len(results)!=2 or [r['variant'] for r in results]!=['installed','corrected'] or data['water_comparison.json']['native_G_qualified']:raise ValueError('Both declared integration variants required')
    differences={}
    for row in results:
        values=row['values']
        if len(values)!=2 or [v['sorted'] for v in values]!=[True,False] or values[0]['grid_points']!=values[1]['grid_points']:raise ValueError('Matched sorted/unsorted grid counts required')
        number=values[0]['numeric_density_integral']-values[1]['numeric_density_integral'];energy=values[0]['XC_energy_Ha']-values[1]['XC_energy_Ha']
        if not all(math.isfinite(v) for v in [number,energy]) or abs(number-row['density_integral_difference'])>1e-15 or abs(energy-row['XC_energy_difference_Ha'])>1e-15:raise ValueError('Integration differences do not reconstruct')
        differences[row['variant']]=dict(density_integral_difference=number,XC_energy_difference_Ha=energy)
    if max(abs(v) for v in differences['corrected'].values())>limit:raise ValueError('Corrected control fails the unchanged upstream tolerance')
    used=[namespace,case/'plan.json',*paths.values(),placement/'launch.json',placement/'completion.json',Path(spec['sources'])/'water.json',case/'water_implementation/implementation.json',case/'water_implementation/screening_kernel.cu',case/'water_implementation/corrected_host.py']
    hashes={str(p):sha256_file(p) for p in used}
    if any(sha256_file(path)!=value for path,value in plan['source_sha256'].items()):raise ValueError('Control source changed')
    return dict(status='screening_masks_and_fixed_density_integration_independently_checked',namespace=str(namespace),device=launch['device'],mask_comparisons=60,corrected_mask_errors=0,
        integration_differences=differences,upstream_absolute_limit=limit,installed_control_exceeds_limit=max(abs(v) for v in differences['installed'].values())>limit,
        peak_host_RSS_bytes=max(host['peak_tree_RSS_bytes'],host['self_VmHWM_bytes']),SCF_gradient_Hessian_qualification_complete=False,production_default_changed=False,native_G_qualified=False,
        proceed_to_SCF_force_qualification=True,source_sha256=hashes,producer_sha256=sha256_file(__file__))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--spec',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=2);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out.exists():raise FileExistsError('Preserve earlier qualification reports')
    r=review(a.spec.resolve());atomic_write_json(a.out,r);print(json.dumps({k:r[k] for k in ['status','device','installed_control_exceeds_limit','integration_differences','proceed_to_SCF_force_qualification']}))
