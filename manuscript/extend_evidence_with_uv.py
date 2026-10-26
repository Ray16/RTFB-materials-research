#!/usr/bin/env python
"""Extend the frozen NMR epoch with independently reviewed DIPEA UV evidence."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time

ROOT=Path(__file__).resolve().parents[1]


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def extend(name):
    if not name or Path(name).name!=name:raise ValueError('One snapshot name required')
    pointer=json.loads((ROOT/'manuscript/evidence/current.json').read_text())
    parent=ROOT/pointer['manifest'];old=json.loads(parent.read_text());old_evidence=parent.parent;old_figures=ROOT/pointer['figures']
    evidence=ROOT/'manuscript/evidence'/name;figures=ROOT/'results/figures/manuscript'/name
    if evidence.exists() or figures.exists():raise FileExistsError('Never overwrite an evidence epoch')
    mapping={}
    for path,expected in old['snapshot_sha256'].items():
        source=ROOT/path
        if sha(source)!=expected:raise ValueError('Parent snapshot changed: '+path)
        if source.parent==old_evidence:dest=evidence/source.name
        elif source.parent==old_figures:dest=figures/source.name
        else:raise ValueError('Unexpected parent snapshot location: '+path)
        mapping[dest]=source
    extra_root=ROOT/'results/validation/dipea_uv_anchor_2026-10-26'
    extras={'dipea_uv_anchor.json':extra_root/'summary.json','dipea_uv_review.json':extra_root/'independent_review.json',
        'dipea_uv_tests.xml':extra_root/'tests.xml','dipea_uv_source.py':ROOT/'scripts/validation/assess_dipea_uv_anchor.py',
        'dipea_uv_review_source.py':ROOT/'scripts/validation/review_dipea_uv_anchor.py'}
    assessment=json.loads(extras['dipea_uv_anchor.json'].read_text());review=json.loads(extras['dipea_uv_review.json'].read_text())
    if (assessment['status']!='primary_DIPEA_UV_protonation_constant_source_checked_conditional_comparison'
            or review['status']!='independent_UV_source_and_comparison_review_checked'):
        raise ValueError('The UV assessment and independent review must pass')
    for record in (assessment,review):
        for path,expected in record['source_sha256'].items():
            if sha(path)!=expected:raise ValueError('Experimental source changed: '+path)
    mapping.update({evidence/key:source for key,source in extras.items()})
    identities={str(p.relative_to(ROOT)):sha(p) for p in set(mapping.values())}
    evidence.mkdir();figures.mkdir()
    for destination,source in mapping.items():shutil.copyfile(source,destination)
    if any(sha(ROOT/path)!=expected for path,expected in identities.items()):
        raise RuntimeError('Source changed during extension; do not publish this epoch')
    manifest=dict(kind='immutable_draft_evidence_snapshot',captured_UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        parent_manifest=pointer['manifest'],parent_manifest_sha256=sha(parent),source_sha256=identities,
        snapshot_sha256={str(p.relative_to(ROOT)):sha(p) for p in mapping},
        scope='Previous frozen evidence plus independent 25 C DIPEA UV measurement; no primary reaction changed')
    (evidence/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return dict(manifest=str((evidence/'manifest.json').relative_to(ROOT)),figures=str(figures.relative_to(ROOT)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--name',required=True);a=p.parse_args()
    print(json.dumps(extend(a.name),indent=2))
