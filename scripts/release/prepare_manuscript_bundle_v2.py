#!/usr/bin/env python
"""Prepare the reviewer bundle after checking its projected-thermal closure."""
import argparse
import hashlib
import json
from pathlib import Path,PurePosixPath
import subprocess

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api

CONFIG=ROOT/'config/manuscript_reviewer_bundle_v2.json'
PRODUCER=ROOT/'scripts/release/prepare_manuscript_bundle.py'


def included(path,declared):
    candidate=PurePosixPath(path)
    if candidate.is_absolute() or '..' in candidate.parts:return False
    return any(candidate==PurePosixPath(base) or PurePosixPath(base) in candidate.parents for base in declared)


def check_git_closure(root,config,head):
    declared=config['source_paths']
    if not included('library/manifest.csv',declared):raise ValueError('Candidate library manifest is not included')
    for needed in ['src/redox/qm/xtb_runtime.py','scripts/pipeline/run_portable_thermal.py']:
        if not included(needed,declared):raise ValueError('Qualified xTB-prefix code is not included')
        subprocess.check_output(['git','cat-file','-e',head+':'+needed],cwd=root)
    data=subprocess.check_output(['git','show',head+':config/projected_thermal_validation.json'],cwd=root)
    certificate=json.loads(data);expected=certificate.pop('certificate_content_sha256')
    if hashlib.sha256(json.dumps(certificate,sort_keys=True,separators=(',',':')).encode()).hexdigest()!=expected:
        raise ValueError('Projected thermal certificate content differs')
    for path,digest in certificate['evidence_sha256'].items():
        if not included(path,declared):raise ValueError('Missing projected certificate dependency: '+path)
        content=subprocess.check_output(['git','show',head+':'+path],cwd=root)
        if hashlib.sha256(content).hexdigest()!=digest:raise ValueError('Projected certificate dependency changed: '+path)
    return len(certificate['evidence_sha256'])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    import os
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    config=json.loads(CONFIG.read_text());head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    count=check_git_closure(ROOT,config,head)
    producer=api.load_module('immutable_previous_bundle_producer',PRODUCER);result=producer.prepare(ROOT,CONFIG,a.out)
    if result['source_commit']!=head:raise ValueError('Source commit changed after closure validation')
    print(json.dumps({'status':result['status'],'source_commit':head,'projected_certificate_dependencies':count}))


if __name__=='__main__':main()
