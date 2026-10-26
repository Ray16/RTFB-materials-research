#!/usr/bin/env python
"""Stage a local hash-qualified xTB prefix for isolated relocation tests.

This is a local test artifact, not an environment lock or redistributable
binary package. Existing qualified files are never modified.
"""
import argparse
import errno
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

from redox.qm.xtb_runtime import QUALIFIED_XTB_BINARY_SHA256,validate_prefix


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def parse_ldd(text,source_prefix):
    dependencies={};system={}
    for line in text.splitlines():
        if 'not found' in line:raise ValueError('Unresolved xTB dependency: '+line.strip())
        match=re.match(r'^\s*(\S+)\s+=>\s+(/\S+)',line)
        if not match:continue
        name,path=match.groups()
        if '/' in name or name in ('.','..'):raise ValueError('Invalid library name')
        resolved=Path(path).resolve()
        if not resolved.is_file():raise ValueError('Missing dependency: '+path)
        target=dependencies if resolved.is_relative_to(source_prefix) else system
        if name in target and target[name]!=resolved:raise ValueError('Conflicting dependency resolution')
        target[name]=resolved
    return dependencies,system


def prepare(source,out):
    from redox.qm import qrrho
    source=source.resolve();out=out.resolve()
    project=Path(__file__).resolve().parents[2]
    if out.is_relative_to(source) or out.is_relative_to(project):
        raise ValueError('Staged prefix must be outside the source environment and project')
    if out.exists():raise FileExistsError('Preserve existing staged prefixes')
    runtime=validate_prefix(source,qrrho.LIBRARY_SHA256)
    dependencies={};system={}
    for path in [Path(runtime['binary']),Path(runtime['library'])]:
        result=subprocess.run(['ldd',str(path)],capture_output=True,text=True,timeout=30,check=True)
        deps,syslibs=parse_ldd(result.stdout,source)
        for name,value in deps.items():
            if name in dependencies and dependencies[name]!=value:raise ValueError('Conflicting dependency across artifacts')
            dependencies[name]=value
        system.update({name:str(value) for name,value in syslibs.items()})
    dependencies['libxtb.so.6.7.1']=Path(runtime['library'])
    out.mkdir(parents=True);(out/'lib').mkdir();(out/'bin').mkdir()
    entries=[]
    try:
        for name,original in sorted(dependencies.items()):
            target=out/'lib'/name
            try:os.link(original,target);method='hardlink'
            except OSError as error:
                if error.errno!=errno.EXDEV:raise
                shutil.copy2(original,target);method='copy'
            if sha(target)!=sha(original):raise ValueError('Staged dependency changed')
            entries.append(dict(path=str(target.relative_to(out)),source=str(original),sha256=sha(target),method=method))
        target=out/'bin/xtb'
        try:os.link(runtime['binary'],target);method='hardlink'
        except OSError as error:
            if error.errno!=errno.EXDEV:raise
            shutil.copy2(runtime['binary'],target);method='copy'
        entries.append(dict(path='bin/xtb',source=runtime['binary'],sha256=sha(target),method=method))
        bound=validate_prefix(out,qrrho.LIBRARY_SHA256)
        if bound['binary_sha256']!=QUALIFIED_XTB_BINARY_SHA256:raise ValueError('Wrong staged executable')
        env=dict(os.environ,LD_LIBRARY_PATH=str(out/'lib'))
        for path in [out/'bin/xtb',out/'lib/libxtb.so.6.7.1']:
            result=subprocess.run(['ldd',str(path)],env=env,capture_output=True,text=True,timeout=30,check=True)
            source_dependencies,_=parse_ldd(result.stdout,source)
            if source_dependencies:
                raise ValueError('Staged runtime still resolves a source-prefix library')
        record=dict(status='local_qualified_xtb_prefix_staged_for_relocation_test',checked=time.time(),
            source_prefix=str(source),staged_prefix=str(out),runtime=bound,files=entries,system_dependencies=system,
            source_read_only=True,public_binary_redistribution=False,scientific_qualification_extended=False,
            producer_sha256=sha(__file__))
        (out/'manifest.json').write_text(json.dumps(record,indent=2)+'\n');return record
    except BaseException as error:
        (out/'PREPARATION_FAILED.txt').write_text(repr(error)+'\n');raise


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    result=prepare(a.source,a.out);print(json.dumps({'status':result['status'],'files':len(result['files'])}))
