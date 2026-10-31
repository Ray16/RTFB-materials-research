#!/usr/bin/env python3
"""Measure retained native outputs for capacity planning, never GPU admission.

Observed completed sizes are neither peak bounds nor storage reservations.
The estimate excludes new electronic stages and other users' future writes.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import time

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json


def measure(census_path,replay_path,volume):
    census=json.loads(census_path.read_text());replay=json.loads(replay_path.read_text())
    if census['status']!='distinct_source_replayed_native_state_census' or replay['status']!='namespace_native_physical_records_replayed':
        raise ValueError('A checked census and completed namespace replay are required')
    if not any(Path(row['path']).resolve()==replay_path and row['sha256']==sha256_file(replay_path) for row in census['source_replays']):
        raise ValueError('The sample replay is not bound by this exact census')
    keys={key for group in census['groups'] for key in group['checked_keys']}
    bindings={}
    sources={str(path):sha256_file(path) for path in (census_path,replay_path,Path(__file__).resolve())}
    for binding in replay['bindings']:
        for kind in ('namespace','review'):
            path=Path(binding[kind]);digest=binding[kind+'_sha256']
            if sha256_file(path)!=digest:raise ValueError('Changed sampled namespace/review')
            sources[str(path)]=digest
        spec=json.loads(Path(binding['namespace']).read_text())
        if spec['key'] in bindings:raise ValueError('Repeated sample namespace')
        bindings[spec['key']]=spec
    measured=[];classes={}
    for row in replay['records']:
        if row['key'] not in keys or row['method_id']!=census['method_id'] or row['status']!='relocated_native_physical_replay_passed':
            raise ValueError('The sample is not a checked state of the declared census')
        spec=bindings[row['key']];raw=Path(spec['raw_parent'])/spec['key'];inventory={}
        for path in raw.rglob('*'):
            if path.is_symlink():raise ValueError('Do not extrapolate storage from linked output trees')
            if path.is_file():
                stat=path.stat();inventory[str(path.relative_to(raw))]=dict(bytes=stat.st_size,allocated_bytes=stat.st_blocks*512,mtime_ns=stat.st_mtime_ns)
        if not inventory:raise ValueError('The completed raw output tree is empty')
        allocated=sum(item['allocated_bytes'] for item in inventory.values());grafts=len(row['occupied_sites'])
        classes.setdefault(grafts,[]).append(allocated)
        measured.append(dict(key=row['key'],grafts=grafts,raw=str(raw),files=len(inventory),allocated_bytes=allocated,logical_bytes=sum(item['bytes'] for item in inventory.values()),inventory=inventory))
    missing={}
    for group in census['groups']:
        if group['checked']!=len(group['checked_keys']) or not 0<=group['checked']<=group['expected']:raise ValueError('Invalid census group counts')
        grafts=len(group['occupied_sites']);missing[grafts]=missing.get(grafts,0)+group['expected']-group['checked']
    if sum(missing.values())!=census['expected_states']-census['checked_states']:raise ValueError('Census totals differ')
    unmeasured=[g for g,n in missing.items() if n and g not in classes]
    estimate=None if unmeasured else sum(n*max(classes[g]) for g,n in missing.items() if n)
    for path,digest in sources.items():
        if sha256_file(path)!=digest:raise ValueError('Storage measurement source changed')
    for sample in measured:
        raw=Path(sample['raw']);current={str(p.relative_to(raw)) for p in raw.rglob('*') if p.is_file()}
        if current!=set(sample['inventory']):raise ValueError('Sample output inventory changed')
        for name,old in sample['inventory'].items():
            stat=(raw/name).stat()
            if (stat.st_size,stat.st_blocks*512,stat.st_mtime_ns)!=(old['bytes'],old['allocated_bytes'],old['mtime_ns']):raise ValueError('Sample size changed during measurement')
    return dict(status='source_bound_retained_native_storage_measurement_not_peak_or_admission_qualification',UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),volume=str(volume),free_bytes=shutil.disk_usage(volume).free,
        measured=measured,missing_native_inputs_by_graft_count=missing,unmeasured_pending_graft_counts=unmeasured,
        observed_retained_allocated_bytes_by_graft_count={g:dict(samples=len(v),minimum=min(v),maximum=max(v)) for g,v in classes.items()},
        remaining_retained_native_bytes_at_observed_class_maximum=estimate,source_sha256=sources,
        temporary_or_in_flight_growth_fully_accounted=False,upstream_energy_outputs_included=False,measured_final_sizes_are_peak_bounds=False,
        other_users_future_writes_controlled=False,admission_policy_changed=False,quantum_jobs_launched=False)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('census','sample-replay','volume','out'):p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if args.out.exists():raise FileExistsError('Preserve earlier storage observations')
    report=measure(args.census.resolve(),args.sample_replay.resolve(),args.volume.resolve());atomic_write_json(args.out,report)
    print(json.dumps(dict(status=report['status'],samples=len(report['measured']),remaining_native_GiB=None if report['remaining_retained_native_bytes_at_observed_class_maximum'] is None else report['remaining_retained_native_bytes_at_observed_class_maximum']/2**30)))
