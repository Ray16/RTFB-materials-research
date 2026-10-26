#!/usr/bin/env python
"""Prepare remaining frozen local states using measured graph/size classes."""
import argparse
import json
import math
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def prepare(construction,witness_path,exclude_reviews,out):
    witness=json.loads(witness_path.read_text());api.verify_sources(witness['source_sha256'])
    expected_atoms=witness['task']['atoms'];peak=witness['peak_host_RSS_bytes']/2**30
    # Resource guards only: a minimum 8 GiB and >4x observed RSS headroom.
    bound=max(8,8*math.ceil((4*peak+2)/8))
    exclude=set()
    for path in exclude_reviews:
        row=json.loads(path.read_text());api.verify_sources(row['source_sha256'])
        plan=api.read_plan(Path(row['campaign']));exclude.add(plan['geometry_sha256'])
    manifest=json.loads((construction/'manifest.json').read_text());selections=[]
    for case in manifest['entries']:
        for mol in case['molecules']:
            if mol['atoms']==expected_atoms and mol['geometry_sha256'] not in exclude:
                selections.append(case['key']+':'+mol['id'])
    if not selections:raise ValueError('No uncomputed local geometries in the measured class')
    old=api.load_module('measured_local_fresh_state_preparer',ROOT/'scripts/pipeline/prepare_local_motif_followup.py')
    old.prepare(construction,witness_path,selections,out)
    campaign=json.loads((out/'campaign.json').read_text())
    for entry in campaign['entries']:
        case=Path(entry['campaign']);plan=api.read_plan(case)
        if any((case/name).exists() for name in ('launch.json','electronic','residency_guard')):
            raise ValueError('Cannot change resource declarations after a scientific attempt')
        plan.update(maximum_host_RSS_GiB=bound,resource_bound_experiment=True,
            resource_bound_derivation=dict(measured_peak_GiB=peak,reference=str(witness_path),
                rule='At least 8 GiB, rounded up to 8 GiB above 4x measured peak plus 2 GiB'),
            reference_energy_transferred=False)
        for p in (Path(__file__).resolve(),ROOT/'scripts/pipeline/dispatch_local_motif_energy.py'):
            plan['producer_sha256'][str(p)]=api.digest(p)
        atomic_write_json(case/'plan.json',plan);entry['plan_sha256']=api.digest(case/'plan.json')
    campaign.update(status='remaining_same_class_local_energies_prepared',maximum_RSS_GiB=bound,
        minimum_free_RAM_GiB=32,minimum_GPU_MiB=15360,created=time.time(),
        worker=str(ROOT/'scripts/pipeline/local_motif_energy.py'),dispatcher=str(ROOT/'scripts/pipeline/dispatch_local_motif_energy.py'),
        source_sha256={str(p):api.digest(p) for p in (Path(__file__).resolve(),construction/'manifest.json',witness_path,
            ROOT/'scripts/pipeline/prepare_local_motif_followup.py',ROOT/'scripts/pipeline/dispatch_local_motif_energy.py')})
    # Test a far-pose pair before expanding the remaining 38-atom class.
    if expected_atoms==38:
        sentinel=next(e['key'] for e in campaign['entries'] if e['key']=='occ_2_4__pose0__cap109__AB')
        campaign['far_pair_resource_sentinel']=sentinel
    atomic_write_json(out/'campaign.json',campaign)
    return campaign


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('construction','resource-review','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--exclude-review',type=Path,action='append',default=[]);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    row=prepare(a.construction.resolve(),a.resource_review.resolve(),[p.resolve() for p in a.exclude_review],a.out.resolve())
    print(json.dumps(dict(prepared=len(row['entries']),host_bound_GiB=row['maximum_RSS_GiB'])))


if __name__=='__main__':main()
