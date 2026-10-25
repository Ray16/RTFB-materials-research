#!/usr/bin/env python
"""Prepare remaining conformers after each occupied graph's reviewed pilot.

The completed input-1 pilot remains an external exact-state reference. Only
the other nineteen declared inputs are optimized; search hit counts and pilot
energies never become conformer weights or values for a different geometry.
"""
import argparse
import json
import os
from pathlib import Path
import shutil

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def prepare(reviews, sources, out):
    builder=api.load_module('reviewed_occupied_conformer_resource_expansion',ROOT/'scripts/pipeline/prepare_smooth_conformer_batch.py')
    choices=[]
    for path in sorted(reviews.glob('*.json')):
        report=json.loads(path.read_text())
        if report.get('status')!='alternate_smooth_geometry_and_resource_pilot_independently_checked':continue
        reference=Path(report['campaign']);point=api.read_plan(reference);task=point['tasks'][0]
        if not task['occupied_sites']:continue
        source=sources/task['model_key'];original=api.read_plan(source)
        if len(original['tasks'])!=20 or original['tasks'][0]['key']!=task['key']:
            raise ValueError('Pilot must be the exact first input of the declared twenty-start graph')
        api.verify_sources(report['source_sha256'])
        choices.append((path,reference,task,source,original))
    if not choices:raise ValueError('No reviewed occupied resource class is ready')
    if len({t['model_key'] for _,_,t,_,_ in choices})!=len(choices):raise ValueError('Duplicate occupied graph resource references')
    out.mkdir(parents=True,exist_ok=False);published=[]
    for review,reference,task,source,old in choices:
        key=task['model_key'];filtered=out/'filtered_inputs'/key;filtered.mkdir(parents=True)
        remaining=[t for t in old['tasks'] if t['key']!=task['key']]
        if len(remaining)!=19:raise ValueError('Exactly nineteen uncomputed input identities are required')
        for item in remaining:shutil.copytree(source/'live_inputs'/item['key'],filtered/'live_inputs'/item['key'])
        inputs=dict(old['input_sha256']);inputs[str(source/'plan.json')]=api.digest(source/'plan.json')
        inputs[str(review)]=api.digest(review)
        producers=dict(old['producer_sha256']);producers[str(Path(__file__).resolve())]=api.digest(__file__)
        plan=dict(old,tasks=remaining,input_sha256=inputs,producer_sha256=producers,
            requested_initial_starts=19,original_twenty_start_source=str(source),
            existing_checked_reference=str(reference),cohort_size_including_reference=20,
            purpose='Remaining nineteen geometry inputs; exact completed pilot retained separately')
        atomic_write_json(filtered/'plan.json',plan)
        batch=out/'batches'/(out.name+'__'+key)
        builder.prepare([filtered],reference,review,batch)
        published.append(dict(model_key=key,new_geometry_tasks=19,existing_checked_reference=str(reference),
            resource_review=str(review),batch=str(batch),batch_sha256=api.digest(batch/'batch.json')))
    atomic_write_json(out/'manifest.json',dict(status='reviewed_occupied_classes_expanded',graphs=published,
        new_geometry_tasks=19*len(published),existing_exact_references=len(published),
        conformer_search_hit_weights_used=False,conformer_convergence_established=False,
        quantum_results_transferred=False,source_sha256={str(p):api.digest(p) for p in (
            Path(__file__).resolve(),ROOT/'scripts/pipeline/prepare_smooth_conformer_batch.py')}))
    print(json.dumps(dict(graphs=len(published),new_geometry_tasks=19*len(published),references_not_recomputed=len(published))))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('reviews','sources','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    prepare(a.reviews.resolve(),a.sources.resolve(),a.out.resolve())


if __name__=='__main__':main()
