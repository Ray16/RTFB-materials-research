#!/usr/bin/env python
"""Assemble explicit counterpoise attempts only after independent AB validation."""
import argparse
import json
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

COLLECTOR=ROOT/'scripts/validation/collect_counterpoise_manifest.py'


def classify(status,owner_live):
    state=status.get('status')
    if state=='refined_AB_independent_stationarity_passed':return 'ready'
    if state in ('refined_AB_independent_stationarity_failed_preserved','refinement_ended_without_reviewable_result','refinement_observer_failed_preserved'):
        return 'upstream_failed_preserved'
    if owner_live:return 'waiting_for_independent_AB_check'
    return 'observer_stopped_without_qualified_result'


def watch(manifest,observer,out,hours):
    if out.exists():raise FileExistsError('Preserve previous collection follow-ups')
    sources={str(p):api.digest(p) for p in (manifest,COLLECTOR,Path(__file__).resolve())}
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    out.mkdir(parents=True);record=dict(status='waiting_for_independent_AB_check',started=time.time(),deadline=time.time()+hours*3600,
        owner=dict(host=os.uname().nodename,pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]),source_sha256=sources,full_counterpoise_computed=False)
    try:
        while time.time()<record['deadline']:
            api.verify_sources(sources)
            path=observer/'status.json';upstream=json.loads(path.read_text());owner=upstream['owner']
            if owner['host']!=os.uname().nodename:raise ValueError('Cannot establish remote observer ownership')
            proc=Path('/proc')/str(owner['pid'])/'stat'
            fields=proc.read_text().rsplit(')',1)[1].split() if proc.exists() else None
            live=fields is not None and fields[19]==owner['process_start'] and fields[0] not in ('Z','X')
            state=classify(upstream,live);record.update(status=state,checked=time.time(),upstream_owner=owner)
            if state=='ready':
                result_path=Path(upstream['result'])
                if api.digest(result_path)!=upstream['result_sha256']:raise ValueError('Independent AB result changed')
                result=api.load_module('five_term_counterpoise_checked_collection',COLLECTOR).collect(manifest,out/'interaction.json')
                record.update(status='five_term_counterpoise_checked_collection_complete',result=str(out/'interaction.json'),result_sha256=api.digest(out/'interaction.json'),full_counterpoise_computed=result['full_counterpoise_computed']);return
            if state!='waiting_for_independent_AB_check':return
            atomic_write_json(out/'status.json',record);time.sleep(15)
        record['status']='finite_collection_observation_ended_sources_preserved'
    except BaseException as error:
        record.update(status='checked_collection_failed_preserved',error=repr(error));raise
    finally:
        record['finished']=time.time();atomic_write_json(out/'status.json',record)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('manifest','observer','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--max-hours',type=float,default=24);a=p.parse_args()
    if not 0<a.max_hours<=48:raise ValueError('Finite observation window required')
    watch(a.manifest.resolve(),a.observer.resolve(),a.out.resolve(),a.max_hours)
