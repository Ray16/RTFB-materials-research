#!/usr/bin/env python
"""Freeze a checked TEMPO route and local potential beside the prior epoch."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def extend(name):
    if not name or Path(name).name != name:
        raise ValueError('One snapshot name required')
    pointer=json.loads((ROOT/'manuscript/evidence/current.json').read_text())
    parent=ROOT/pointer['manifest'];previous=json.loads(parent.read_text())
    old_evidence,old_figures=parent.parent,ROOT/pointer['figures']
    evidence,figures=ROOT/'manuscript/evidence'/name,ROOT/'results/figures/manuscript'/name
    if evidence.exists() or figures.exists():raise FileExistsError('Never overwrite an evidence epoch')
    mapping={}
    for path,expected in previous['snapshot_sha256'].items():
        source=ROOT/path
        if sha(source)!=expected:raise ValueError('Parent snapshot changed: '+path)
        if source.parent==old_evidence:dest=evidence/source.name
        elif source.parent==old_figures:dest=figures/source.name
        else:raise ValueError('Unexpected parent snapshot location: '+path)
        mapping[dest]=source
    base=ROOT/'results/grafting/tempo_exact_support_2026-10-26'
    extras={'tempo_route.json':base/'route_review.json','tempo_route_screen.json':base/'route_screen.json',
            'tempo_native_potential.json':base/'native_potential_review.json',
            'tempo_decision.json':base/'all_axis_decision.json','tempo_summary.json':base/'summary.json',
            'tempo_registry_reconciliation.json':base/'registry_reconciliation.json'}
    route,screen,potential,decision,summary=[json.loads(extras[k].read_text()) for k in
        ('tempo_route.json','tempo_route_screen.json','tempo_native_potential.json','tempo_decision.json','tempo_summary.json')]
    if (route['candidate']['id']!='tempo' or screen['candidates'][0]['route_status']!='proceed'
            or potential['status']!='TEMPO_local_oxidation_native_partition_independently_reconstructed'
            or decision['decision']['status']!='pending' or decision['decision']['all_required_axes_pass']
            or summary['actual_route_deltaG_computed'] or summary['quantum_jobs_submitted']):
        raise ValueError('TEMPO route, native potential or material decision is not reviewable')
    for record in (route,potential,decision):
        for path,expected in record['source_sha256'].items():
            if sha(path)!=expected:raise ValueError('TEMPO source changed: '+path)
    for key,source in extras.items():mapping[evidence/key]=source
    figure=ROOT/'results/figures/grafting/tempo_exact_route.pdf'
    if sha(figure)!=summary['figure_pdf_sha256']:raise ValueError('TEMPO figure changed')
    mapping[figures/'tempo_exact_route.pdf']=figure
    identities={str(p.relative_to(ROOT)):sha(p) for p in set(mapping.values())}
    evidence.mkdir();figures.mkdir()
    for destination,source in mapping.items():shutil.copyfile(source,destination)
    if any(sha(ROOT/path)!=expected for path,expected in identities.items()):
        raise RuntimeError('Source changed during extension; do not publish this epoch')
    manifest=dict(kind='immutable_draft_evidence_snapshot',captured_UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
        source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        parent_manifest=pointer['manifest'],parent_manifest_sha256=sha(parent),source_sha256=identities,
        snapshot_sha256={str(p.relative_to(ROOT)):sha(p) for p in mapping},
        scope='Prior evidence plus exact-support TEMPO route, local oxidation model and pending material decision')
    (evidence/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return dict(manifest=str((evidence/'manifest.json').relative_to(ROOT)),figures=str(figures.relative_to(ROOT)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--name',required=True);a=p.parse_args()
    print(json.dumps(extend(a.name),indent=2))
