#!/usr/bin/env python
"""Extend the TEMPO epoch with the independently replayed raw-state addendum."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time
from xml.etree import ElementTree

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
    release=ROOT/'results/releases/tempo_addendum';current=json.loads((release/'current.json').read_text())
    relocation=json.loads((release/'relocation.json').read_text())
    package=ROOT/current['prepared_package'];archive=ROOT/current['compressed_archive']
    if (current['status']!='TEMPO_raw_archive_and_two_relocated_replays_checked'
            or relocation['status']!='TEMPO_archive_and_two_relocated_replays_checked'
            or current['raw_states']!=4 or current['Hessians_reconstructed']!=4
            or current['public_archive_deposited'] or current['public_DOI_created']
            or sha(archive)!=current['archive_sha256'] or sha(package/'manifest.json')!=current['package_manifest_sha256']
            or relocation['archive_sha256']!=current['archive_sha256']
            or relocation['native_result']['native_reconstruction']['Hessians_reconstructed']!=4
            or relocation['stdlib_isolation']['original_project_open_attempts']
            or relocation['native_isolation']['original_project_open_attempts']):
        raise ValueError('TEMPO raw archive or isolated replay is not qualified')
    tests=ElementTree.parse(release/'negative_tests.xml').getroot()
    suites=[tests] if tests.tag=='testsuite' else list(tests.iter('testsuite'))
    if sum(int(s.get('tests','0')) for s in suites)!=6 or any(int(s.get('failures','0')) or int(s.get('errors','0')) for s in suites):
        raise ValueError('TEMPO addendum adverse-input checks did not pass')
    extras={'tempo_archive.json':release/'current.json','tempo_archive_relocation.json':release/'relocation.json',
            'tempo_archive_manifest.json':package/'manifest.json','tempo_archive_tests.xml':release/'negative_tests.xml'}
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
        scope='Prior TEMPO evidence plus isolated raw-state arithmetic and Hessian replay; no new quantum results')
    (evidence/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return dict(manifest=str((evidence/'manifest.json').relative_to(ROOT)),figures=str(figures.relative_to(ROOT)))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--name',required=True);a=p.parse_args()
    print(json.dumps(extend(a.name),indent=2))
