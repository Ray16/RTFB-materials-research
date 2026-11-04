#!/usr/bin/env python3
"""Stream, verify and partially replay a bounded local review collection."""
import argparse
import hashlib
import json
import os
from pathlib import Path,PurePosixPath
import subprocess
import sys
import tarfile

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json


def finalize(folder,archive,extracted):
    if archive.exists() or extracted.exists():raise FileExistsError('Preserve earlier bundle/archive checks')
    metadata=json.loads((folder/'manifest.json').read_text())
    for name,wanted in metadata['files_sha256'].items():
        if sha256_file(folder/name)!=wanted:raise ValueError('Changed collection input: '+name)
    producer=subprocess.Popen(['tar','-cf','-','-C',str(folder.parent),folder.name],stdout=subprocess.PIPE)
    try:subprocess.run(['zstd','-T2','-3','-q','-o',str(archive)],stdin=producer.stdout,check=True)
    finally:producer.stdout.close()
    if producer.wait()!=0:raise RuntimeError('Collection tar export failed')
    expected={folder.name+'/'+p:h for p,h in metadata['files_sha256'].items()}
    expected[folder.name+'/manifest.json']=sha256_file(folder/'manifest.json');found={}
    decoder=subprocess.Popen(['zstd','-dc',str(archive)],stdout=subprocess.PIPE)
    with tarfile.open(fileobj=decoder.stdout,mode='r|') as stream:
        for member in stream:
            if member.isdir():continue
            if not member.isfile() or member.name not in expected:raise ValueError('Unexpected archive member')
            digest=hashlib.sha256();source=stream.extractfile(member)
            while block:=source.read(2**20):digest.update(block)
            found[member.name]=digest.hexdigest()
    decoder.stdout.close()
    if decoder.wait()!=0 or found!=expected:raise ValueError('Streamed outer archive differs from checked inputs')
    extracted.mkdir(parents=True)
    prefix=PurePosixPath('validation/screening_chemical_impact_2026-10-31')
    decoder=subprocess.Popen(['zstd','-dc',str(folder/'source.tar.zst')],stdout=subprocess.PIPE)
    files=0
    with tarfile.open(fileobj=decoder.stdout,mode='r|') as stream:
        for member in stream:
            name=PurePosixPath(member.name)
            if '..' in name.parts or name.is_absolute():raise ValueError('Unsafe source archive path')
            if not name.is_relative_to(prefix):continue
            if member.isdir():continue
            if not member.isfile():raise ValueError('Only actual cached replay files may be extracted')
            target=extracted.joinpath(*name.relative_to(prefix).parts);target.parent.mkdir(parents=True,exist_ok=True)
            with stream.extractfile(member) as source,target.open('xb') as output:
                while block:=source.read(2**20):output.write(block)
            files+=1
    decoder.stdout.close()
    if decoder.wait()!=0:raise RuntimeError('Committed source archive extraction failed')
    code='''import json,pathlib,runpy,sys
blocked=[pathlib.Path(p).resolve() for p in sys.argv[2:]]
def deny(event,args):
 if event=="open" and args and isinstance(args[0],(str,bytes)):
  path=pathlib.Path(args[0]).resolve()
  if any(path==base or path.is_relative_to(base) for base in blocked):raise PermissionError("Original project reads blocked")
sys.addaudithook(deny)
runpy.run_path(sys.argv[1],run_name="__main__")
'''
    result=subprocess.run([sys.executable,'-I','-S','-c',code,str(extracted/'replay.py'),str(ROOT),str(ROOT.parent/'archive'),
        str(ROOT.parent/'redox_kernel_controls'),str(ROOT.parent/'redox_candidate_continuations')],capture_output=True,text=True,check=True)
    replay=json.loads(result.stdout)
    if replay['status']!='stdlib_screening_evidence_replay_passed':raise ValueError('Extracted numerical replay failed')
    return dict(status='reviewer_collection_archive_integrity_and_extracted_stdlib_scope_checked',
        collection=str(folder),archive=str(archive),archive_sha256=sha256_file(archive),archive_bytes=archive.stat().st_size,
        nested_archives=9,outer_members_checked=len(found),extracted_numerical_files=files,source_commit=metadata['source_commit'],
        extracted_replay=replay,original_project_archive_and_raw_paths_blocked=True,
        all_nested_scientific_replays_repeated=False,public_deposition_complete=False,data_license=None,author_metadata=None,
        submission_ready=False,producer_sha256=sha256_file(__file__))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('collection','archive','extract','out'):parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(4,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(4))
    if args.out.exists():raise FileExistsError('Preserve previous collection checks')
    result=finalize(args.collection.resolve(),args.archive.resolve(),args.extract.resolve());atomic_write_json(args.out,result)
    print(json.dumps({k:result[k] for k in ('status','archive_bytes','outer_members_checked','extracted_replay','public_deposition_complete')}))
