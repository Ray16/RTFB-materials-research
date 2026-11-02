#!/usr/bin/env python3
"""Verify every archived control byte and replay after extraction with source reads blocked."""
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


def check(package,archive,extract):
    if extract.exists():raise FileExistsError('Preserve earlier extraction/replay records')
    data=json.loads((package/'manifest.json').read_text());expected={package.name+'/'+name:digest for name,digest in data['files_sha256'].items()}
    expected[package.name+'/manifest.json']=sha256_file(package/'manifest.json')
    extract.mkdir(parents=True);found={};decoder=subprocess.Popen(['zstd','-dc',str(archive)],stdout=subprocess.PIPE)
    with tarfile.open(fileobj=decoder.stdout,mode='r|') as stream:
        for member in stream:
            name=PurePosixPath(member.name)
            if name.is_absolute() or '..' in name.parts or not name.is_relative_to(PurePosixPath(package.name)):raise ValueError('Unsafe archive member')
            if member.isdir():continue
            if not member.isfile() or member.name not in expected or member.name in found:raise ValueError('Unexpected or duplicate archive file')
            target=extract.joinpath(*name.parts);target.parent.mkdir(parents=True,exist_ok=True)
            digest=hashlib.sha256()
            with stream.extractfile(member) as source,target.open('xb') as output:
                for block in iter(lambda:source.read(2**20),b''):digest.update(block);output.write(block)
            found[member.name]=digest.hexdigest()
    decoder.stdout.close()
    if decoder.wait()!=0 or found!=expected:raise ValueError('Archive byte inventory differs from checked package')
    blocked=[str(ROOT),str(ROOT.parent/'archive'),str(ROOT.parent/'redox_kernel_controls'),str(ROOT.parent/'redox_candidate_continuations'),str(package)]
    code='''import importlib.metadata as m,json,pathlib,runpy,sys
root=pathlib.Path(sys.argv[1]);blocked=[pathlib.Path(p).resolve() for p in sys.argv[2:]]
def deny(event,args):
 if event=="open" and args and isinstance(args[0],(str,bytes)):
  path=pathlib.Path(args[0]).resolve()
  if any(path==p or path.is_relative_to(p) for p in blocked):raise PermissionError("Original paths blocked")
sys.addaudithook(deny)
namespace=runpy.run_path(str(root/"code/replay.py"),run_name="portable_review")
r=namespace["replay"](root)
r["replay_versions"]={name:m.version(name) for name in ("numpy","h5py","pyscf")}
r["GPU_libraries_imported"]=any(name in sys.modules for name in ("cupy","gpu4pyscf"))
assert not r["GPU_libraries_imported"]
print(json.dumps(r))
'''
    env=dict(os.environ)
    for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[key]='1'
    result=subprocess.run([sys.executable,'-I','-c',code,str(extract/package.name),*blocked],cwd='/',env=env,capture_output=True,text=True,check=True)
    replay=json.loads(result.stdout)
    return dict(status='component_archive_all_bytes_and_relocated_CPU_replay_checked',archive=str(archive),archive_sha256=sha256_file(archive),archive_bytes=archive.stat().st_size,
        members_verified=len(found),extracted_package=str(extract/package.name),replay=replay,original_scientific_and_prepared_package_reads_blocked=True,
        public_deposition_complete=False,data_license=None,author_metadata=None,producer_sha256=sha256_file(__file__))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('package','archive','extract','out'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(4,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(4))
    if a.out.exists():raise FileExistsError('Preserve previous archive checks')
    report=check(a.package.resolve(),a.archive.resolve(),a.extract.resolve());atomic_write_json(a.out,report)
    print(json.dumps(dict(status=report['status'],members_verified=report['members_verified'],replay=report['replay'])))
