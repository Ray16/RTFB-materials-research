#!/usr/bin/env python
"""Connect checked occupied-graph resource points to all matching geometries."""
import argparse
import json
import math
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def prepare(batches,pilot_batch,pilot_reviews,reference,resource_review,out,bound):
    helper=api.load_module('occupied_electronic_existing_dependency_contract',ROOT/'scripts/pipeline/prepare_smooth_electronic_batch.py')
    point=api.read_plan(reference);witness=json.loads(resource_review.read_text());api.verify_sources(witness['source_sha256'])
    helper.resource_compatible(point,point,witness)
    peak=witness['peak_host_RSS_bytes']/2**30
    if not math.isfinite(bound) or not 1.25*peak<=bound<=point['maximum_host_RSS_GiB']:
        raise ValueError('Successor RAM guard must retain measured headroom and not exceed the qualified pilot ceiling')
    out.mkdir(parents=True,exist_ok=False);declared=[]
    for batch in batches:
        manifest=json.loads((batch/'batch.json').read_text());api.verify_sources(manifest['source_sha256'])
        for entry in manifest['entries']:helper.resource_compatible(api.read_plan(Path(entry['campaign'])),point,witness)
        destination=out/(out.name+'__'+batch.name)
        helper.declare(batch,ROOT/'results/reviews/2026-10-25/smooth_occupied_expansion_review/checked',reference,resource_review,destination)
        declared.append(destination)
    # Retain the independently checked first pose for each configuration as
    # an explicit geometry subset, without launching another geometry job.
    old_path=pilot_batch/'batch.json';old=json.loads(old_path.read_text());matching=[];reviews={}
    for entry in old['entries']:
        plan=api.read_plan(Path(entry['campaign']))
        try:helper.resource_compatible(plan,point,witness)
        except ValueError:continue
        path=pilot_reviews/(pilot_batch.name+'__'+entry['key']+'.json')
        if not path.exists():raise ValueError('A matching pilot geometry lacks its independent review')
        matching.append(entry);reviews[entry['key']]=str(path)
    if not matching:raise ValueError('No matching pilot geometries in the declared source batch')
    subset=out/'matching_geometry_pilots';subset.mkdir()
    atomic_write_json(subset/'batch.json',dict(old,entries=matching,source_sha256={**old['source_sha256'],str(old_path):api.digest(old_path)}))
    destination=out/(out.name+'__pilot_electronic_states');helper.declare(subset,pilot_reviews,reference,resource_review,destination);declared.append(destination)
    totals=0;existing=0
    for batch in declared:
        path=batch/'batch.json';manifest=json.loads(path.read_text());manifest.update(maximum_RSS_GiB=bound,
            resource_budget_derivation=dict(measured_peak_GiB=peak,guard_GiB=bound,headroom_ratio=bound/peak,
                numerical_protocol_changed=False),parent_followup=str(out))
        manifest['source_sha256'][str(Path(__file__).resolve())]=api.digest(__file__)
        if batch.name.endswith('__pilot_electronic_states'):
            for entry in manifest['entries']:entry['geometry_review']=reviews[entry['key']]
            for review in reviews.values():manifest['source_sha256'][review]=api.digest(review)
        atomic_write_json(path,manifest);totals+=len(manifest['entries']);existing+=manifest['existing_reference_entries']
    proposed=set()
    for batch in declared:
        for entry in json.loads((batch/'batch.json').read_text())['entries']:
            if entry['existing_checked_reference']:continue
            case=Path(entry['campaign']);raw=ROOT/'calcs/polymer_models'/case.parent.name/case.name
            source_raw=Path(api.read_plan(Path(entry['geometry_campaign']))['raw_root'])
            if raw==source_raw or str(raw) in proposed or raw.exists():
                raise ValueError('Electronic namespace collides with an existing/source calculation: '+str(raw))
            proposed.add(str(raw))
    report=dict(status='occupied_same_constitution_electronic_followup_declared',created=time.time(),
        batches=[str(p) for p in declared],declared_electronic_states=totals,identical_checked_references=existing,
        fresh_electronic_states=totals-existing,maximum_host_RSS_GiB=bound,measured_reference_peak_GiB=peak,
        source_sha256={str(p):api.digest(p) for p in [Path(__file__).resolve(),reference/'plan.json',resource_review,*[b/'batch.json' for b in declared]]},
        primary_candidate_admission=False,native_G_qualification_inferred=False)
    atomic_write_json(out/'manifest.json',report);return report


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--geometry-batch',type=Path,action='append',required=True)
    for name in ('pilot-batch','pilot-reviews','reference','resource-review','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--host-bound',type=float,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    result=prepare([p.resolve() for p in a.geometry_batch],a.pilot_batch.resolve(),a.pilot_reviews.resolve(),a.reference.resolve(),a.resource_review.resolve(),a.out.resolve(),a.host_bound)
    print(json.dumps({k:result[k] for k in ['status','declared_electronic_states','fresh_electronic_states','maximum_host_RSS_GiB']}))


if __name__=='__main__':main()
