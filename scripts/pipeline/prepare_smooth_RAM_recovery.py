#!/usr/bin/env python
"""Prepare bounded RAM experiments for preserved alternate-energy attempts.

The larger ceiling is prospective resource headroom, not a new scientific
protocol. Exact coordinates, gas/SMD equations and all energy/QC gates remain.
"""
import argparse
import json
import math
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def prepare(failed, reference, resource_review, out, bound):
    if not math.isfinite(bound) or not 80<bound<=128:raise ValueError('Use a finite bounded RAM experiment above 80 and at most 128 GiB')
    owner=api.load_module('RAM_recovery_stopped_owner',ROOT/'scripts/pipeline/run_polymer_spacing.py')
    original=api.load_module('RAM_recovery_exact_geometry_preparation',ROOT/'scripts/pipeline/smooth_polymer_energy.py')
    point=api.read_plan(reference);resource=json.loads(resource_review.read_text())
    api.verify_sources(resource['source_sha256'])
    out.mkdir(parents=True,exist_ok=False);entries=[]
    for source in failed:
        old=api.read_plan(source);api.verify_sources(old['producer_sha256']);api.verify_sources(old['input_sha256'])
        state_path=source/'experimental_energy/status.json';state=json.loads(state_path.read_text())
        guard_path=source/'host_guard/status.json';guard=json.loads(guard_path.read_text())
        if (state['status']!='experimental_hybrid_energy_failed_preserved' or owner.owner_alive(state['owner'])
                or guard['maximum_RSS_GiB']!=80 or guard['peak_RSS_bytes']<=80*2**30
                or guard['free_RAM_bytes']<32*2**30 or (source/'residency_guard/interruption.json').exists()):
            raise ValueError('A stopped, preserved 80-GiB RSS exceedance with free-RAM headroom is required')
        geometry,review=Path(old['source_geometry_campaign']),Path(old['geometry_review'])
        key=old['tasks'][0]['key'];case=out/key
        original.prepare(geometry,review,case);plan=api.read_plan(case)
        if api.digest(case/'input.xyz')!=api.digest(source/'input.xyz'):
            raise ValueError('Resource recovery changed the exact scored coordinates')
        for p in (source/'plan.json',state_path,guard_path,resource_review,reference/'plan.json'):
            plan['input_sha256'][str(p)]=api.digest(p)
        plan['input_sha256'].update(resource['source_sha256'])
        plan['producer_sha256'][str(Path(__file__).resolve())]=api.digest(__file__)
        plan.update(execution_driver=str(ROOT/'scripts/pipeline/smooth_polymer_energy_bounded.py'),
            maximum_host_RSS_GiB=bound,prospective_resource_pilot=True,
            measured_alternate_model_electronic_reference=str(reference),
            preserved_failed_attempt=str(source),electronic_RAM_class_review=str(resource_review),
            resource_scope='Prospective larger RAM ceiling after observed RSS exceedance; exact geometry/equations/QC retained')
        atomic_write_json(case/'plan.json',plan)
        entries.append(dict(key=key,campaign=str(case),geometry_campaign=str(geometry),
            geometry_plan_sha256=api.digest(geometry/'plan.json'),geometry_review=str(review),
            existing_checked_reference=False,plan_sha256=None))
    scripts=[Path(__file__).resolve(),ROOT/'scripts/pipeline/prepare_smooth_electronic_batch.py',
        ROOT/'scripts/pipeline/coordinate_smooth_electronic_batch.py',ROOT/'scripts/pipeline/smooth_polymer_energy_bounded.py',
        ROOT/'scripts/pipeline/smooth_polymer_energy.py',ROOT/'scripts/pipeline/dispatch_smooth_energy_pilot.py',
        ROOT/'scripts/pipeline/coordinate_geometry_batch.py',ROOT/'environments/cutensor22/runtime.json',resource_review]
    atomic_write_json(out/'batch.json',dict(status='prospective_RAM_recovery_experiments_prepared',stage='electronic',
        entries=entries,worker=str(ROOT/'scripts/pipeline/smooth_polymer_energy_bounded.py'),
        preparer=str(ROOT/'scripts/pipeline/prepare_smooth_electronic_batch.py'),reference=str(reference),
        resource_review=str(resource_review),maximum_RSS_GiB=bound,minimum_free_RAM_GiB=32,minimum_GPU_MiB=32768,
        experimental_method_id=point['experimental_method_id'],primary_candidate_admission=False,native_G_computed=False,
        source_sha256={str(p):api.digest(p) for p in scripts},existing_reference_entries=0,created=time.time(),
        quantum_values_transferred_to_other_geometries=False))
    print(json.dumps(dict(prepared=len(entries),prospective_host_RSS_ceiling_GiB=bound,old_attempts_preserved=True)))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--failed',type=Path,action='append',required=True)
    for name in ('reference','resource-review','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--bound',type=float,default=128);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    prepare([v.resolve() for v in a.failed],a.reference.resolve(),a.resource_review.resolve(),a.out.resolve(),a.bound)


if __name__=='__main__':main()
