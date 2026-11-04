#!/usr/bin/env python3
"""Count distinct physically replayed polymer states without inferring ensembles."""
import argparse
import json
import math
import os
from pathlib import Path

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json


def census(replays,configurations,site_cases,count):
    rows=[];sources={};keys=set();methods=set();graphs={}
    if len(configurations)!=len(set(configurations)) or len(site_cases)!=len(set(site_cases)) or count<1:
        raise ValueError('Distinct declared graph/configuration cases and positive input count required')
    for path in replays:
        data=json.loads(path.read_text())
        if (data['status'] not in ('archived_polymer_native_physical_records_replayed','namespace_native_physical_records_replayed')
                or data['states']!=len(data['records']) or data['quantum_calculations_run'] or data['global_ensemble_convergence_proven']):
            raise ValueError('Complete, scoped physical replay required')
        for producer,digest in data['producer_sha256'].items():
            if sha256_file(producer)!=digest:raise ValueError('Replay producer changed')
            sources[producer]=digest
        sources[str(path.resolve())]=sha256_file(path)
        for row in data['records']:
            group=(row['configuration'],tuple(row['occupied_sites']))
            if (row['key'] in keys or group[0] not in configurations or group[1] not in site_cases
                    or row['status']!='relocated_native_physical_replay_passed'
                    or set(row['G_by_cutoff_eV'])!={'25','50','100'}
                    or any(type(v) not in (int,float) or not math.isfinite(v) for v in row['G_by_cutoff_eV'].values())
                    or not math.isfinite(row['maximum_column_error']) or not 0<=row['maximum_column_error']<=.01
                    or not math.isfinite(row['minimum_frequency_cm']) or row['minimum_frequency_cm']<-5):
                raise ValueError('Duplicate, undeclared, incomplete or unqualified native state')
            label='_'.join(map(str,group[1])) if group[1] else 'none'
            expected_prefix=f'{group[0]}__occ_{label}__conf'
            if not row['key'].startswith(expected_prefix):raise ValueError('State key differs from its graph/configuration')
            index=int(row['key'][len(expected_prefix):])
            if not 1<=index<=count:raise ValueError('Conformer index is outside the declared campaign')
            if graphs.setdefault(group,row['graph_sha256'])!=row['graph_sha256']:raise ValueError('Different exact graphs share a census cell')
            methods.add(row['method_id']);keys.add(row['key']);rows.append(row)
    if len(methods)!=1:raise ValueError('One explicit solvent model per census required')
    groups=[]
    for configuration in configurations:
        for sites in site_cases:
            selected=[r['key'] for r in rows if r['configuration']==configuration and tuple(r['occupied_sites'])==sites]
            label='_'.join(map(str,sites)) if sites else 'none';expected=[f'{configuration}__occ_{label}__conf{i:03d}' for i in range(1,count+1)]
            groups.append(dict(configuration=configuration,occupied_sites=list(sites),checked=len(selected),expected=count,
                checked_keys=sorted(selected),missing_keys=sorted(set(expected)-set(selected)),
                complete_declared_batch=len(selected)==count,ensemble_convergence_established=False))
    if any(sha256_file(p)!=h for p,h in sources.items()):raise ValueError('Replay evidence changed during census')
    return dict(status='distinct_source_replayed_native_state_census',checked_states=len(rows),expected_states=len(groups)*count,
        groups=groups,method_id=next(iter(methods)),source_replays=[dict(path=str(p.resolve()),sha256=sha256_file(p)) for p in replays],
        source_sha256=sources,global_ensemble_convergence_proven=False,balanced_occupied_ensemble_cycle_established=False,
        experimental_occupancy_distribution_inferred=False,primary_candidate_admission=False)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--replay',type=Path,action='append',required=True)
    p.add_argument('--configuration',action='append');p.add_argument('--sites',action='append');p.add_argument('--count',type=int,default=20)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out.exists():raise FileExistsError('Preserve previous census snapshots')
    sites=[tuple() if s=='none' else tuple(map(int,s.split(','))) for s in (a.sites or ['none','2','3','4','2,3','2,4'])]
    result=census([p.resolve() for p in a.replay],a.configuration or ['all_same','alternating'],sites,a.count)
    result['producer_sha256']=sha256_file(__file__);a.out.parent.mkdir(parents=True,exist_ok=True);atomic_write_json(a.out,result)
    print(json.dumps(dict(status=result['status'],checked_states=result['checked_states'],expected_states=result['expected_states'])))
