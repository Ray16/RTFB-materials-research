#!/usr/bin/env python
"""Stage two saved spin controls for an isolated projected-thermal CLI check."""
import argparse
import csv
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import time


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare(source,controls,out):
    source,controls,out=source.resolve(),controls.resolve(),out.resolve()
    if out.exists() or out.is_relative_to(source):
        raise ValueError('Use a new directory outside the source checkout')
    out.mkdir(parents=True)
    ignore=shutil.ignore_patterns('__pycache__','*.pyc')
    for name in ('src','config'):shutil.copytree(source/name,out/name,ignore=ignore)
    scripts=('scripts/pipeline/run_relocated_projected.py',
             'scripts/release/run_relocated_projected_isolated.py',
             'scripts/validation/check_relocated_native_guards.py')
    for relative in scripts:
        target=out/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source/relative,target)
    spec=importlib.util.spec_from_file_location('projected_certificate_closure',source/'scripts/release/certificate_closure.py')
    closure=importlib.util.module_from_spec(spec);spec.loader.exec_module(closure)
    closure.closure(source,out,out/'certificate_closure.json')
    rows=[];states=[]
    for name in ('formaldehyde','methyl_radical'):
        original=controls/(name+'__baseline')/name/'rest';target=out/'calcs/dft'/name/'rest';target.mkdir(parents=True)
        for filename in ('opt.xyz','result.json'):shutil.copyfile(original/filename,target/filename)
        shutil.copytree(original/'native_thermal',target/'native_thermal',ignore=ignore)
        references=list((original/'projected_thermal').glob('*.json'))
        if len(references)!=1:raise ValueError('One exact reference partition required per control')
        states.append(dict(id=name,state='rest',source=str(original),target=str(target),
                           reference_partition=str(references[0]),reference_sha256=sha(references[0])))
        rows.append(dict(id=name,state='rest'))
    (out/'library').mkdir()
    with (out/'library/manifest.csv').open('w') as stream:
        writer=csv.DictWriter(stream,fieldnames=['id','state']);writer.writeheader();writer.writerows(rows)
    record=dict(status='isolated_native_partition_relocation_inputs_staged',created=time.time(),source=str(source),
                checkout=str(out),states=states,producer_sha256=sha(Path(__file__)),
                files_sha256={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()})
    (out/'staging_manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('source','controls','out'):parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args();r=prepare(args.source,args.controls,args.out)
    print(json.dumps(dict(status=r['status'],states=len(r['states']))))
