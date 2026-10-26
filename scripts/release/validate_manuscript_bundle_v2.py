#!/usr/bin/env python
"""Replay a consolidated bundle and verify its relocated certificate closure."""
import argparse
import hashlib
import json
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_closure(project):
    manifest=project/'library/manifest.csv'
    if not manifest.is_file():raise ValueError('Relocated candidate manifest is missing')
    for p in ['src/redox/qm/xtb_runtime.py','scripts/pipeline/run_portable_thermal.py']:
        if not (project/p).is_file():raise ValueError('Relocated xTB-prefix entry point is missing')
    path=project/'config/projected_thermal_validation.json';certificate=json.loads(path.read_text());expected=certificate.pop('certificate_content_sha256')
    if hashlib.sha256(json.dumps(certificate,sort_keys=True,separators=(',',':')).encode()).hexdigest()!=expected:
        raise ValueError('Relocated projected-thermal certificate changed')
    for relative,value in certificate['evidence_sha256'].items():
        target=(project/relative).resolve()
        if not target.is_relative_to(project.resolve()) or not target.is_file() or sha(target)!=value:
            raise ValueError('Relocated projected certificate dependency differs: '+relative)
    return dict(candidate_manifest_sha256=sha(manifest),projected_certificate_sha256=sha(path),
                projected_certificate_dependencies=len(certificate['evidence_sha256']))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--archive',type=Path,required=True);p.add_argument('--bundle-name',required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--python',type=Path,default=Path('/usr/bin/python3'));a=p.parse_args()
    base=api.load_module('immutable_previous_bundle_validator',ROOT/'scripts/release/validate_manuscript_bundle.py')
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    import os
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    report=base.validate(a.archive,a.bundle_name,a.out,a.python)
    closure=verify_closure(a.out.resolve()/'project');closure.update(status='relocated_candidate_manifest_and_projected_certificate_closure_checked',
        checked=time.time(),bundle_sha256=report['archive_sha256'],source_commit=report['source_commit'],
        native_GPU_or_projected_entrypoint_requalified=False,experimental_accuracy_inferred=False,validator_sha256=sha(__file__))
    (a.out.resolve()/'closure.json').write_text(json.dumps(closure,indent=2)+'\n');print(json.dumps(closure,indent=2))


if __name__=='__main__':main()
