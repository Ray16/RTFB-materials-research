#!/usr/bin/env python3
"""Lock the qualified xTB build and its installed conda runtime closure."""
import argparse
import json
import os
from pathlib import Path
import re
from urllib.parse import urlsplit

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.qm import qrrho
from redox.qm.xtb_runtime import validate_prefix


def dependency_closure(records):
    selected={};pending=['xtb']
    while pending:
        name=pending.pop()
        if name in selected or name.startswith('__'):continue
        if name not in records:raise ValueError('Missing installed runtime dependency: '+name)
        row=records[name];url=row.get('url','');parsed=urlsplit(url)
        if (parsed.scheme!='https' or parsed.hostname!='conda.anaconda.org'
                or not parsed.path.startswith('/conda-forge/') or parsed.username or parsed.password or parsed.query
                or not re.fullmatch('[a-f0-9]{64}',row.get('sha256',''))
                or not re.fullmatch('[a-f0-9]{32}',row.get('md5',''))):
            raise ValueError('Public conda-forge URL and complete artifact hashes required: '+name)
        selected[name]={k:row[k] for k in ('name','version','build','url','md5','sha256','depends')}
        pending.extend(value.split()[0] for value in row['depends'])
    return [selected[name] for name in sorted(selected)]


def prepare(prefix,output):
    runtime=validate_prefix(prefix,qrrho.LIBRARY_SHA256)
    if output.exists():raise FileExistsError('Preserve the existing environment lock')
    metadata={};sources={}
    for path in (prefix/'conda-meta').glob('*.json'):
        data=json.loads(path.read_text());name=data['name']
        if name in metadata:raise ValueError('Duplicate installed package: '+name)
        metadata[name]=data;sources[name]=path
    packages=dependency_closure(metadata);output.mkdir(parents=True)
    lines=['# Qualified xTB thermodynamic runtime; linux-64; all runtime dependencies explicit.','@EXPLICIT']
    (output/'environment.explicit.txt').write_text('\n'.join(lines+[row['url']+'#'+row['md5'] for row in packages])+'\n')
    record=dict(status='qualified_xtb_conda_runtime_artifacts_locked',runtime=runtime,packages=packages,
        artifact_SHA256_check_required=True,virtual_package_requirements=sorted({d for row in packages for d in row['depends'] if d.startswith('__')}),
        source_sha256={str(sources[row['name']]):sha256_file(sources[row['name']]) for row in packages},
        producer_sha256=sha256_file(__file__),full_quantum_environment_qualified=False)
    atomic_write_json(output/'packages.json',record)
    return dict(status=record['status'],packages=len(packages),explicit_file=str(output/'environment.explicit.txt'))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--prefix',type=Path,required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    print(json.dumps(prepare(a.prefix.resolve(),a.out.resolve())))
