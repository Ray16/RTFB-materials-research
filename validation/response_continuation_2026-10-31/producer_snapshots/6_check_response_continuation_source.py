#!/usr/bin/env python3
"""Read-only check of a full-Hessian response prefix, including live ownership.

This never starts, stops or qualifies a quantum calculation. A live owner makes
the source ineligible for continuation even when its completed arrays bind.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.qm.hessian_restart import read_continuation_frame
from redox.qm.response_resume import inspect_source
from redox.workflows.molecular_native_namespace import verify
from redox.workflows.screening_control import response_continuation_closure


def check(namespace):
    producers=[Path(__file__),Path(read_continuation_frame.__code__.co_filename),Path(inspect_source.__code__.co_filename),Path(response_continuation_closure.__code__.co_filename)]
    hashes={str(p):sha256_file(p) for p in producers}
    spec=verify(namespace);case=Path(spec['case']);stage=json.loads((case/'status.json').read_text())
    definition=importlib.util.spec_from_file_location('response_prefix_owners',ROOT/'scripts/pipeline/observe_archived_native_namespace.py')
    owners=importlib.util.module_from_spec(definition);definition.loader.exec_module(owners)
    live=owners.alive(stage['owner']);frame,_,source=read_continuation_frame(case);bound=inspect_source(case,frame)
    placements=list(case.parent.glob('placement*/completion.json'));reason='live_or_unreachable_owner' if live is not False else 'no_closed_completion'
    eligible=False
    if len(placements)==1:
        try:
            reason=response_continuation_closure(stage,json.loads((case/'residency_guard/owner.json').read_text()),
                json.loads(placements[0].read_text()),lambda owner:live)
            eligible=True
        except ValueError as error:reason=str(error)
    if any(sha256_file(path)!=digest for path,digest in hashes.items()):raise ValueError('Source-prefix auditor changed during the check')
    return dict(status='completed_response_prefix_source_binding_checked',checked=time.time(),namespace=str(namespace),
        namespace_sha256=sha256_file(namespace),captured_stage=stage,owner_live=live,closed_source_eligible=eligible,
        closure=reason,original_component_case=str(source),frame=frame,response_prefix=bound,
        GPU_job_started=False,full_Hessian_or_G_qualified=False,source_sha256=hashes,producer_sha256=sha256_file(__file__))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,action='append',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out.exists():raise FileExistsError('Preserve earlier source-prefix audits')
    records=[check(path.resolve()) for path in a.source]
    atomic_write_json(a.out,dict(status='read_only_response_prefix_preflight',records=records,GPU_jobs_submitted=0))
    print(json.dumps([dict(source=r['namespace'],live=r['owner_live'],closed_source_eligible=r['closed_source_eligible'],
        completed_blocks=r['response_prefix']['completed_blocks'],completed_RHS=r['response_prefix']['completed_RHS']) for r in records]))
