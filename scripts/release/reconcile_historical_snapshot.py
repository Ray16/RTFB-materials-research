#!/usr/bin/env python
"""Verify older publication source hashes using explicitly archived replacements."""
import argparse
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def reconcile(records,mapping,out):
    replacements=json.loads(mapping.read_text());checked=[]
    for record in records:
        doc=json.loads(record.read_text())
        for original,expected in doc['source_sha256'].items():
            original_path=ROOT/original
            retained=original_path if original_path.exists() and sha(original_path)==expected else ROOT/replacements.get(original,'__missing__')
            if not retained.is_file() or sha(retained)!=expected:
                raise ValueError('No exact historical source for '+original)
            checked.append(dict(record=str(record.relative_to(ROOT)),original=original,
                retained=str(retained.relative_to(ROOT)),sha256=expected))
    if out.exists():raise FileExistsError('Preserve existing historical verification')
    out.parent.mkdir(parents=True,exist_ok=True)
    result=dict(status='historical_publication_inputs_reconciled_byte_for_byte',
        checked_inputs=len(checked),records=checked,mapping_sha256=sha(mapping),
        reviewer_sha256=sha(__file__),not_a_new_scientific_qualification=True)
    out.write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--record',type=Path,action='append',required=True);p.add_argument('--mapping',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    row=reconcile([r.resolve() for r in a.record],a.mapping.resolve(),a.out.resolve());print(json.dumps({'status':row['status'],'checked_inputs':row['checked_inputs']}))
