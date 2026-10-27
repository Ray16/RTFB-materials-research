#!/usr/bin/env python
"""Prepare the current local-motif manuscript bundle with six raw-state packages."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity


CONFIG=ROOT/'config/manuscript_reviewer_bundle_v4.json'


def check_epoch(root,config,head,closure):
    declared=config['source_paths']
    pointer=json.loads(subprocess.check_output(['git','show',head+':manuscript/evidence/current.json'],cwd=root))
    if not closure.included(pointer['manifest'],declared):raise ValueError('Current manuscript evidence manifest is excluded')
    manifest=json.loads(subprocess.check_output(['git','show',head+':'+pointer['manifest']],cwd=root))
    for path,expected in manifest['snapshot_sha256'].items():
        if not closure.included(path,declared):raise ValueError('Current manuscript snapshot file is excluded: '+path)
        content=subprocess.check_output(['git','show',head+':'+path],cwd=root)
        if hashlib.sha256(content).hexdigest()!=expected:raise ValueError('Committed manuscript snapshot differs: '+path)
    return len(manifest['snapshot_sha256'])


def prepare(out):
    config=json.loads(CONFIG.read_text());head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
    closure=api.load_module('bundle_v3_git_closure',ROOT/'scripts/release/prepare_manuscript_bundle_v2.py')
    dependencies=closure.check_git_closure(ROOT,config,head);evidence_files=check_epoch(ROOT,config,head,closure)
    if {x['id'] for x in config['packages']}!={'primary_review','first_site','numerical_controls','neighbor_cycle','tempo_addendum','local_motif'}:
        raise ValueError('Expected six distinct reviewer packages')
    producer=api.load_module('bundle_v3_existing_producer',ROOT/'scripts/release/prepare_manuscript_bundle.py')
    result=producer.prepare(ROOT,CONFIG,out)
    if result['source_commit']!=head:raise ValueError('Source commit changed after closure validation')
    (out/'README.md').write_text('''# Manuscript reviewer bundle (local preparation)

The committed source snapshot contains the current manuscript evidence epoch.
Six overlapping raw-state packages cover the primary panel, first-site and
neighbor cycles, numerical controls, TEMPO and local-motif electronic records. Counts
are not additive. Extract each archive and use its README for CPU replay.

Cached replay does not repeat original physical QC or establish experimental
accuracy, broad GPU portability or conformer convergence. Third-party source
articles retain their own licenses. Data license, author metadata, public
deposition and DOI are pending.
''')
    result['files_sha256']['README.md']=api.digest(out/'README.md')
    result.update(v4_preparer_sha256=api.digest(__file__),current_evidence_files=evidence_files,
                  projected_certificate_dependencies=dependencies,raw_package_count=6)
    atomic_write_json(out/'manifest.json',result)
    return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    r=prepare(a.out.resolve());print(json.dumps({'status':r['status'],'source_commit':r['source_commit'],
                                                 'packages':len(r['packages']),'evidence_files':r['current_evidence_files']}))
