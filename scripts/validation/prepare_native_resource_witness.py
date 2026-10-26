#!/usr/bin/env python
"""Project a checked native state into a small successor-job resource witness.

Only temporary first-order response arrays in the execution directory are
excluded from successor dependencies. The full original review and arrays stay
in place. Final Hessians, gradients, partitions, FD records and code remain
strictly hash-checked. No verification cache or numerical substitution is used.
"""
import argparse
import json
import os
from pathlib import Path
import re
import time

from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def project_sources(mapping,raw):
    required={};intermediate={};execution=(raw/'execution').resolve()
    for name,value in mapping.items():
        path=Path(name).resolve()
        if path.parent==execution and re.fullmatch(r'(?:rks|uks)_make_h1_\d+\.npy',path.name):
            intermediate[name]=value
        else:required[name]=value
    return required,intermediate


def prepare(source,out):
    original=json.loads(source.read_text())
    if original['status']!='alternate_native_state_independently_reconstructed':raise ValueError('Independently checked native state required')
    api.verify_sources(original['source_sha256'])
    campaign=Path(original['campaign']);plan=api.read_plan(campaign);raw=Path(plan['raw_root'])
    required,intermediate=project_sources(original['source_sha256'],raw)
    # These final-state objects must never disappear from the strict closure.
    for path in (raw/'hessian.npy',raw/'gradient.npy',raw/'geometry.xyz',raw/'experimental_partition.json'):
        if str(path) not in required:raise ValueError('Final-state closure is incomplete: '+str(path))
    for path in (source,Path(__file__).resolve()):required[str(path)]=api.digest(path)
    result=dict(original,source_sha256=required,
        artifact_kind='derived_resource_witness_with_checked_final_state',resource_only_projection=True,
        derivation=dict(prepared=time.time(),original_review=str(source),original_review_sha256=api.digest(source),
            original_full_source_manifest=original['source_sha256'],
            temporary_response_dependencies_omitted=intermediate,
            omitted_files_preserved=True,original_full_source_manifest_verified_at_projection=True,
            numerical_fields_changed=False,final_state_dependencies_all_retained=True),
        scope='Derived successor-job resource witness and unchanged checked final state. Temporary response arrays remain in the full original audit but are not numerical inputs to successor calculations; no new QC or thermochemistry claimed.')
    if out.exists():raise ValueError('Preserve an existing resource witness')
    out.parent.mkdir(parents=True,exist_ok=True);atomic_write_json(out,result);return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--review',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    row=prepare(a.review.resolve(),a.out.resolve());print(json.dumps(dict(status=row['artifact_kind'],
        excluded_intermediate_arrays=len(row['derivation']['temporary_response_dependencies_omitted']),
        strict_final_state_sources=len(row['source_sha256']))))


if __name__=='__main__':main()
