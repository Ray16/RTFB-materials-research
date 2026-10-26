#!/usr/bin/env python
"""Extend an immutable manuscript epoch without re-reading live result feeds."""
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
    old_pointer=json.loads((ROOT/'manuscript/evidence/current.json').read_text())
    old_manifest=ROOT/old_pointer['manifest'];prior=json.loads(old_manifest.read_text())
    old_evidence=old_manifest.parent;old_figures=ROOT/old_pointer['figures']
    evidence=ROOT/'manuscript/evidence'/name;figures=ROOT/'results/figures/manuscript'/name
    if evidence.exists() or figures.exists():raise FileExistsError('Never overwrite an evidence epoch')
    mapping={}
    for path,expected in prior['snapshot_sha256'].items():
        source=ROOT/path
        if sha(source)!=expected:raise ValueError('Parent snapshot changed: '+path)
        if source.parent==old_evidence:dest=evidence/source.name
        elif source.parent==old_figures:dest=figures/source.name
        else:raise ValueError('Unexpected parent snapshot location: '+path)
        mapping[dest]=source
    extras={'dipea_experimental_anchor.json':ROOT/'results/validation/dipea_experimental_anchor_2026-10-26/summary.json',
        'dipea_experimental_anchor_source.py':ROOT/'scripts/validation/assess_dipea_experimental_anchor.py',
        'dipea_experimental_anchor_tests.xml':ROOT/'results/validation/dipea_experimental_anchor_2026-10-26/tests.xml'}
    mapping.update({evidence/name:source for name,source in extras.items()})
    identities={str(p.relative_to(ROOT)):sha(p) for p in set(mapping.values())}
    anchor=json.loads(extras['dipea_experimental_anchor.json'].read_text())
    for path,expected in anchor['source_sha256'].items():
        if sha(path)!=expected:raise ValueError('DIPEA assessment source changed: '+path)
    evidence.mkdir();figures.mkdir()
    for destination,source in mapping.items():shutil.copyfile(source,destination)
    if any(sha(ROOT/path)!=expected for path,expected in identities.items()):
        raise RuntimeError('Source changed during extension; do not publish this epoch')
    manifest=dict(kind='immutable_draft_evidence_snapshot',captured_UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        parent_manifest=old_pointer['manifest'],parent_manifest_sha256=sha(old_manifest),
        source_sha256=identities,snapshot_sha256={str(p.relative_to(ROOT)):sha(p) for p in mapping},
        scope='Previous frozen numerical evidence plus measured DIPEA salt NMR estimate; no primary reaction changed')
    (evidence/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return dict(manifest=str((evidence/'manifest.json').relative_to(ROOT)),figures=str(figures.relative_to(ROOT)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--name',required=True);a=p.parse_args()
    print(json.dumps(extend(a.name),indent=2))
