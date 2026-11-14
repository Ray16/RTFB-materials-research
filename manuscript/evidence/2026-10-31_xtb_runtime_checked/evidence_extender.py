#!/usr/bin/env python3
"""Extend an unchanged draft snapshot while preserving its original bindings."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess

ROOT=Path(__file__).resolve().parents[1]


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def extend(root,pointer,name,additions,scope):
    if not re.fullmatch(r'[A-Za-z0-9_-]+',name):raise ValueError('One safe snapshot directory name required')
    pointer_bytes=pointer.read_bytes();current=json.loads(pointer_bytes);parent=root/current['manifest'];previous=json.loads(parent.read_text())
    for path,expected in previous['snapshot_sha256'].items():
        if digest(root/path)!=expected:raise ValueError('Parent snapshot changed: '+path)
    if len({key for key,_ in additions})!=len(additions) or any(Path(key).name!=key or key in ('manifest.json','evidence_extender.py','..','.') for key,_ in additions):
        raise ValueError('Distinct safe addendum filenames required')
    out=pointer.parent/name
    if out.exists():raise FileExistsError('Preserve the existing snapshot')
    out.mkdir();snapshots=dict(previous['snapshot_sha256']);sources=dict(previous['source_sha256'])
    for original in parent.parent.rglob('*'):
        if not original.is_file() or original==parent:continue
        target=out/original.relative_to(parent.parent);target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copyfile(original,target);snapshots[str(target.relative_to(root))]=digest(target)
    for key,original in additions:
        if (out/key).exists():raise FileExistsError('Addendum would replace inherited evidence: '+key)
        original=original.resolve();expected=digest(original);shutil.copyfile(original,out/key)
        if digest(out/key)!=expected:raise ValueError('Addendum changed during capture')
        sources[str(original.relative_to(root))]=expected;snapshots[str((out/key).relative_to(root))]=expected
    producer=out/'evidence_extender.py';shutil.copyfile(__file__,producer);snapshots[str(producer.relative_to(root))]=digest(producer)
    record=dict(kind='immutable_draft_evidence_snapshot',captured_UTC=datetime.now(timezone.utc).isoformat(),
        source_commit=subprocess.run(['git','rev-parse','HEAD'],cwd=root,capture_output=True,text=True,check=True).stdout.strip(),
        parent_manifest=str(parent.relative_to(root)),parent_manifest_sha256=digest(parent),source_sha256=sources,
        snapshot_sha256=snapshots,inherited_bound_references_sha256=previous['snapshot_sha256'],scope=scope,
        extension_producer_sha256=digest(producer),scientific_qualification_inferred_from_capture=False)
    (out/'manifest.json').write_text(json.dumps(record,indent=2)+'\n')
    for path,expected in snapshots.items():
        if digest(root/path)!=expected:raise ValueError('Snapshot changed before publication: '+path)
    if pointer.read_bytes()!=pointer_bytes:raise ValueError('Current evidence pointer changed during extension')
    current['manifest']=str((out/'manifest.json').relative_to(root));temporary=pointer.with_suffix('.tmp')
    temporary.write_text(json.dumps(current,indent=2)+'\n');temporary.replace(pointer)
    return dict(snapshot=str(out/'manifest.json'),bound_files=len(snapshots))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--name',required=True);p.add_argument('--add',action='append',default=[])
    p.add_argument('--scope',required=True);a=p.parse_args();additions=[]
    for item in a.add:
        key,value=item.split('=',1);additions.append((key,Path(value)))
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    import os
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    print(json.dumps(extend(ROOT,ROOT/'manuscript/evidence/current.json',a.name,additions,a.scope)))
