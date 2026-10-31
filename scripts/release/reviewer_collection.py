#!/usr/bin/env python3
"""Assemble a local review draft from checked archives and committed source."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json


def assemble(archive,output,commit):
    if output.exists() or output.is_relative_to(ROOT) or output.is_relative_to(archive):
        raise FileExistsError('Use a new shared-storage release draft outside scientific checkouts')
    commit=subprocess.check_output(['git','rev-parse',commit+'^{commit}'],cwd=ROOT,text=True).strip()
    rows=[]
    for name in ('primary_review','first_site','controls','neighbor_cycle','tempo_addendum','local_motif'):
        record=archive/'results/releases'/name/'current.json';data=json.loads(record.read_text())
        source=archive/data['compressed_archive']
        rows.append(dict(id=name,source=str(source),sha256=data['archive_sha256'],record=str(record)))
    for name,record_path in [('native_inputs75',ROOT/'validation/native_cache_2026-10-30/summary.json'),
                             ('exact_viologen_CPU',ROOT/'validation/exact_molecular_CPU_2026-10-31/summary.json')]:
        data=json.loads(record_path.read_text())
        source=(Path(data['archive']) if name=='native_inputs75' else ROOT.parent/'redox_release_candidates/2026-10-31_exact_viologen_CPU_final.tar.zst')
        rows.append(dict(id=name,source=str(source),sha256=data['archive_sha256'],record=str(record_path)))
    if any(sha256_file(r['source'])!=r['sha256'] for r in rows):raise ValueError('A checked package changed before assembly')
    estimated=sum(Path(r['source']).stat().st_size for r in rows)
    if shutil.disk_usage(output.parent).free<32*2**30+2*estimated+2*2**30:
        raise RuntimeError('Shared-volume margin unavailable for the bounded review bundle')
    output.mkdir(parents=True);packages=output/'packages';packages.mkdir();records=output/'package_records';records.mkdir()
    files={}
    for row in rows:
        target=packages/(row['id']+'.tar.zst');shutil.copyfile(row['source'],target)
        if sha256_file(target)!=row['sha256']:raise ValueError('Copied checked archive differs')
        shutil.copyfile(row['record'],records/(row['id']+'.json'))
        files[str(target.relative_to(output))]=row['sha256']
    source=output/'source.tar.zst'
    producer=subprocess.Popen(['git','-c','gc.auto=0','archive','--format=tar',commit],cwd=ROOT,stdout=subprocess.PIPE)
    try:
        subprocess.run(['zstd','-T2','-3','-q','-o',str(source)],stdin=producer.stdout,check=True)
    finally:
        producer.stdout.close()
    if producer.wait()!=0:raise RuntimeError('Committed source export failed')
    files[source.name]=sha256_file(source)
    for row in rows:subprocess.run(['zstd','-t','-q',str(packages/(row['id']+'.tar.zst'))],check=True)
    subprocess.run(['zstd','-t','-q',str(source)],check=True)
    text='''# Local reviewer collection — draft

This collection contains a committed source/evidence snapshot and eight previously
checked raw-data/replay packages. Extract source.tar.zst into an empty directory;
each package has its own README, locked environment and replay entry point.
Package scopes overlap, so their record counts must not be added as unique states.
The source includes a versioned ten-state NumPy CPU example and the standard-Python
71-state electronic/144-column derivative replay under validation/.

The recorded SHA-256 manifest checks every nested archive before scientific replay.
Cached replay does not repeat SCF or establish new physical qualification. Ongoing
full-Hessian/charged-polymer pilots are excluded from accepted scientific records.
Their raw shared-storage records are separate from this committed evidence epoch.

This is a review draft. Author metadata, data-license decisions, final scientific
completion and public DOI deposition remain pending. Experimental papers are linked
and identified under their original rights; no project license is assigned to them.
No public deposition or submission readiness is claimed.
'''
    (output/'README.md').write_text(text)
    metadata=dict(status='local_reviewer_collection_prepared_checked_archives_not_deposited',source_commit=commit,
        packages=rows,files_sha256={str(p.relative_to(output)):sha256_file(p) for p in output.rglob('*') if p.is_file()},
        nested_archives_integrity_checked=True,public_deposition_complete=False,data_license=None,author_metadata=None,
        submission_ready=False,producer_sha256=sha256_file(__file__))
    atomic_write_json(output/'manifest.json',metadata)
    return metadata


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--archive',type=Path,required=True)
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--commit',default='HEAD');args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(4,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(4))
    result=assemble(args.archive.resolve(),args.out.resolve(),args.commit)
    print(json.dumps({k:result[k] for k in ('status','source_commit','nested_archives_integrity_checked','public_deposition_complete')}))
