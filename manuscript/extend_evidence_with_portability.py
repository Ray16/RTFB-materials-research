#!/usr/bin/env python
"""Extend the frozen paper evidence with checked CPU partition relocation."""
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
    if not name or Path(name).name!=name:raise ValueError('One epoch name required')
    pointer=json.loads((ROOT/'manuscript/evidence/current.json').read_text());parent=ROOT/pointer['manifest']
    old=json.loads(parent.read_text());old_figures=ROOT/pointer['figures']
    evidence=ROOT/'manuscript/evidence'/name;figures=ROOT/'results/figures/manuscript'/name
    if evidence.exists() or figures.exists():raise FileExistsError('Preserve frozen epochs')
    mapping={}
    for relative,expected in old['snapshot_sha256'].items():
        source=ROOT/relative
        if sha(source)!=expected:raise ValueError('Parent snapshot changed')
        destination=(evidence/source.name if source.parent==parent.parent else
                     figures/source.name if source.parent==old_figures else None)
        if destination is None:raise ValueError('Unexpected parent snapshot location')
        mapping[destination]=source
    review=ROOT/'results/validation/relocated_projected_2026-10-26';data=json.loads((review/'summary_verified.json').read_text())
    if (data['status']!='relocated_projected_partitions_and_physical_guards_independently_checked'
            or len(data['controls'])!=2 or data['negative_controls']!=10 or data['original_checkout_file_accesses']!=0
            or data['xtb_CLI_executions']!=6 or data['quantum_calculations_run']
            or data['native_Hessian_execution_requalified'] or data['general_platform_portability_inferred']):
        raise ValueError('Relocation scope or independent review changed')
    for path,expected in data['source_sha256'].items():
        if sha(path)!=expected:raise ValueError('Relocation review backing file changed')
    suites=list(ElementTree.parse(review/'tests.xml').getroot().iter('testsuite'))
    if [sum(int(s.get(k,'0')) for s in suites) for k in ('tests','failures','errors','skipped')]!=[22,0,0,0]:
        raise ValueError('Focused tests incomplete')
    mapping[evidence/'projected_portability.json']=review/'summary_verified.json'
    mapping[evidence/'projected_portability_tests.xml']=review/'tests.xml'
    mapping[evidence/'projected_portability_adverse.json']=next(Path(p) for p in data['source_sha256'] if p.endswith('/adverse_inputs.json'))
    identities={str(p.relative_to(ROOT)) if p.is_relative_to(ROOT) else str(p):sha(p) for p in set(mapping.values())}
    evidence.mkdir();figures.mkdir()
    for target,source in mapping.items():shutil.copyfile(source,target)
    if any(sha(ROOT/p)!=expected for p,expected in identities.items()):raise ValueError('Input changed while copying')
    manifest=dict(kind='immutable_draft_evidence_snapshot',captured_UTC=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
                  source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                  parent_manifest=pointer['manifest'],parent_manifest_sha256=sha(parent),source_sha256=identities,
                  snapshot_sha256={str(p.relative_to(ROOT)):sha(p) for p in mapping},
                  scope='Prior chemical results plus explicit source/xTB relocation of saved native-Hessian CPU partitions; no new QM or polymer ensemble result')
    (evidence/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return dict(manifest=str((evidence/'manifest.json').relative_to(ROOT)),figures=str(figures.relative_to(ROOT)))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--name',required=True);args=parser.parse_args()
    print(json.dumps(extend(args.name),indent=2))
