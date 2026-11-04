#!/usr/bin/env python3
"""Extend an immutable declared-source inventory with explicitly hashed seeds."""
import argparse
import json
import os
from pathlib import Path

from redox.core.hashing import sha256_file
from redox.validation.source_inventory import inventory


def extend(root, baseline_summary, seed_manifest, output, summary):
    if output.exists() or summary.exists():
        raise FileExistsError('Preserve earlier inventory observations')
    baseline_record=json.loads(baseline_summary.read_text())
    baseline_path=Path(baseline_record['inventory'])
    if sha256_file(baseline_path)!=baseline_record['inventory_sha256']:
        raise ValueError('Baseline inventory bytes changed')
    baseline=json.loads(baseline_path.read_text())
    seeds=json.loads(seed_manifest.read_text())
    # Known snapshot paths are candidates only: inventory hashes every new use.
    known=[dict(path=str(root/p),sha256=h) for p,h in seeds.get('known_snapshot_sha256',{}).items()]
    result=inventory(root,root.parent/'archive/0_redox',seeds['seeds_sha256'],known,baseline=baseline)
    result.update(baseline_inventory=str(baseline_path),baseline_inventory_sha256=sha256_file(baseline_path),
        addendum_seed_manifest=str(seed_manifest),source_sha256={str(p):sha256_file(p) for p in
        (baseline_summary,seed_manifest,Path(__file__),Path(inventory.__code__.co_filename))})
    output.write_text(json.dumps(result,indent=2)+'\n')
    record=dict(status=result['status'],inventory=str(output),inventory_sha256=sha256_file(output),
        verified_files=len(result['files']),unresolved_bindings=len(result['unresolved']),
        unique_content_bytes=result['unique_content_bytes'],new_files_hash_checked=result['new_files_hash_checked'],
        inherited_files=result['inherited_files'],baseline_raw_files_rehashed=False,
        current_full_raw_archive_created=False,public_deposition_complete=False,
        seed_manifest=str(seed_manifest),seed_manifest_sha256=sha256_file(seed_manifest))
    summary.write_text(json.dumps(record,indent=2)+'\n')
    return record


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('root','baseline-summary','seeds','out','summary'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    print(json.dumps(extend(a.root.resolve(),a.baseline_summary.resolve(),a.seeds.resolve(),a.out.resolve(),a.summary.resolve())))
