#!/usr/bin/env python
"""Full spin/Hessian controls for between-SCF garbage collection."""
import json
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


def checked_controls(root):
    evidence={}
    for mult in (1,2):
        path=Path(root)/'controls'/str(mult)/'status.json'
        row=json.loads(path.read_text())
        if row['status']!='full_original_compact_equivalence_passed' or not row.get('GC_lifecycle_exercised'):
            raise ValueError('Both full spin controls must pass with exercised collection')
        api.verify_sources(row['source_sha256']);api.verify_sources(row['evidence_sha256'])
        evidence[str(path)]=api.digest(path);evidence.update(row['source_sha256'])
    return evidence


def main():
    module=api.load_module('GC_full_spin_template',ROOT/'scripts/validation/qualify_smd_workspace.py')
    old=module.adapted_body
    def body(text):
        text=old(text);key='scripts/pipeline/smd_workspace_execution.py'
        if text.count(key)!=1:raise ValueError('The original replay template changed')
        return text.replace(key,'scripts/pipeline/smd_gc_workspace.py')
    module.adapted_body=body
    module.main()
    out=Path(sys.argv[sys.argv.index('--out')+1]);path=out/'status.json'
    row=json.loads(path.read_text());proof=out/'compact/execution/GC_execution.json';gc=json.loads(proof.read_text())
    if gc['status']!='GC_workspace_finished_requires_control_and_physical_checks' or gc['collections']<1:
        raise ValueError('The control did not exercise the SCF lifecycle hook')
    row['GC_lifecycle_exercised']=True
    row['source_sha256'].update(gc['source_sha256'])
    row['source_sha256'][str(Path(__file__).resolve())]=api.digest(__file__)
    row['evidence_sha256'][str(proof)]=api.digest(proof)
    path.write_text(json.dumps(row,indent=2)+'\n')


if __name__=='__main__':main()
