#!/usr/bin/env python3
"""Freeze direct sampling-analysis inputs without changing source paths or hashes."""
import argparse
import json
import os
from pathlib import Path
import shutil

from redox.core.hashing import sha256_file
from redox.validation.portable_inputs import PortableInputs


def prepare(folder, output, prefix='', existing=None):
    if output.exists():raise FileExistsError('Preserve earlier portable inputs')
    artifact=lambda name:folder/(prefix+name)
    analysis=artifact('finite_batches.json');saved=json.loads(analysis.read_text())
    expected=saved['source_sha256'];analysis_hash=sha256_file(analysis)
    reader=PortableInputs(existing,expected,analysis_hash) if existing else None
    snapshot_by_hash={sha256_file(p):p for p in folder.iterdir() if p.is_file()}
    sources={}
    for original,digest in expected.items():
        if reader is not None:source=reader.verify(original,digest)
        elif digest in snapshot_by_hash:source=snapshot_by_hash[digest]
        else:source=Path(original)
        if sha256_file(source)!=digest:raise ValueError('Source changed: '+original)
        sources[original]=source
    manifest=json.loads(artifact('manifest.json').read_text())
    for name,digest in manifest['files'].items():
        if sha256_file(artifact(name))!=digest:raise ValueError('Publication artifact changed')
    output.mkdir(parents=True);bindings={}
    for original,source in sources.items():
        digest=expected[original];name='input_'+digest+Path(original).suffix
        if not (output/name).exists():shutil.copyfile(source,output/name)
        if sha256_file(output/name)!=digest:raise ValueError('Copied input changed')
        bindings[original]=dict(file=name,sha256=digest)
    for name in ['manifest.json',*manifest['files'],*(['table.tex'] if artifact('table.tex').exists() else [])]:
        shutil.copyfile(artifact(name),output/name)
    record=dict(schema=1,analysis_sha256=analysis_hash,bindings=bindings,
        scope='Direct source-bound sampling-analysis inputs; no new physical QC or transitive raw archive',
        producer_sha256=sha256_file(__file__))
    (output/'portable_inputs.json').write_text(json.dumps(record,indent=2)+'\n')
    return record


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--folder',type=Path,required=True);parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--prefix',default='');parser.add_argument('--existing',type=Path)
    args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    result=prepare(args.folder.resolve(),args.out.resolve(),args.prefix,args.existing)
    print(json.dumps(dict(bindings=len(result['bindings']),analysis_sha256=result['analysis_sha256'])))
