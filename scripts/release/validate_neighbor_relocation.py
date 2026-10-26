#!/usr/bin/env python3
"""Extract/hash-check the neighbor archive and replay with original reads blocked."""
import argparse
import hashlib
import json
import os
from pathlib import Path,PurePosixPath
import runpy
import shutil
import subprocess
import sys
import tarfile
import time


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def worker(package,out,forbidden):
    blocked=[]
    def audit(event,args):
        if event=='open' and isinstance(args[0],(str,bytes)):
            value=os.fsdecode(args[0]);path=Path(value).resolve()
            if path.is_relative_to(forbidden):
                blocked.append(str(path));raise RuntimeError('Original project reads are forbidden')
    sys.addaudithook(audit)
    result=out/'arithmetic.json';sys.argv=[str(package/'reproduce.py'),'--data',str(package),'--out',str(result)]
    runpy.run_path(str(package/'reproduce.py'),run_name='__main__')
    loaded=[n for n in sys.modules if n.split('.')[0] in {'numpy','redox','pyscf','gpu4pyscf','cupy','torch','ase'}]
    if blocked or loaded:raise RuntimeError('Replay attempted original reads or loaded scientific/production packages')
    (out/'isolation.json').write_text(json.dumps(dict(status='stdlib_neighbor_replay_passed_original_reads_blocked',
        checked=time.time(),attempted_original_reads=blocked,forbidden_modules_loaded=loaded,
        site_imports_disabled=bool(sys.flags.no_site),isolated_python=bool(sys.flags.isolated),
        interpreter=sys.executable,CPU_affinity=sorted(os.sched_getaffinity(0)),
        result_sha256=sha(result),validator_sha256=sha(__file__)),indent=2)+'\n')


def validate(source,archive,out,python):
    root=Path(__file__).resolve().parents[2]
    if out.is_relative_to(root):raise ValueError('Extract outside the original checkout')
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    data=json.loads((source/'manifest.json').read_text())
    expected={**data['files_sha256'],'manifest.json':sha(source/'manifest.json')}
    out.mkdir(parents=True,exist_ok=False);extracted=out/'extracted';extracted.mkdir();members={}
    process=subprocess.Popen(['zstd','-dc',str(archive)],stdout=subprocess.PIPE)
    try:
        with tarfile.open(fileobj=process.stdout,mode='r|') as stream:
            for item in stream:
                parts=PurePosixPath(item.name).parts
                if (not parts or parts[0]!=source.name or '..' in parts or PurePosixPath(item.name).is_absolute()
                        or not (item.isdir() or item.isfile())):raise ValueError('Unsafe archive member')
                target=extracted.joinpath(*parts)
                if item.isdir():target.mkdir(parents=True,exist_ok=True);continue
                relative=str(PurePosixPath(*parts[1:]))
                if relative not in expected or relative in members:raise ValueError('Unexpected or duplicate archive file')
                target.parent.mkdir(parents=True,exist_ok=True)
                with stream.extractfile(item) as incoming,target.open('xb') as outgoing:shutil.copyfileobj(incoming,outgoing)
                value=sha(target)
                if value!=expected[relative]:raise ValueError('Archived bytes differ')
                members[relative]=value
        process.stdout.close()
        if process.wait()!=0 or members!=expected:raise ValueError('Incomplete archive inventory')
    finally:
        if process.poll() is None:process.terminate();process.wait()
    helper=out/'validate_source.py';shutil.copyfile(Path(__file__).resolve(),helper)
    env=dict(os.environ);env.pop('PYTHONPATH',None);env.pop('LD_LIBRARY_PATH',None)
    for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):env[key]='1'
    args=[str(python),'-B','-I','-S',str(helper),'--worker','--package',str(extracted/source.name),
          '--out',str(out),'--forbid-root',str(root)]
    with (out/'replay.log').open('x') as log:subprocess.run(args,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    result=json.loads((out/'arithmetic.json').read_text());isolation=json.loads((out/'isolation.json').read_text())
    if result['states']!=14 or result['balanced_reactions']!=14 or not isolation['site_imports_disabled'] or not isolation['isolated_python']:
        raise ValueError('Relocated replay did not meet its declared scope')
    proof=dict(status='neighbor_archive_inventory_and_relocated_stdlib_replay_checked',checked=time.time(),
        archive=str(archive),archive_sha256=sha(archive),archive_bytes=archive.stat().st_size,
        archive_regular_files=len(members),package_manifest_sha256=sha(source/'manifest.json'),
        raw_states=14,balanced_reactions=14,rotor_cutoffs_cm=[25,50,100],isolation=isolation,
        arithmetic=result,validator_sha256=sha(__file__),public_archive_deposited=False,
        nuclear_partition_reintegrated=False,conformer_convergence_established=False)
    (out/'validation.json').write_text(json.dumps(proof,indent=2)+'\n');return proof


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--worker',action='store_true')
    for name in ('package','archive','out','forbid-root','stdlib-python'):p.add_argument('--'+name,type=Path)
    a=p.parse_args()
    if a.worker:worker(a.package.resolve(),a.out.resolve(),a.forbid_root.resolve())
    else:
        if not all([a.package,a.archive,a.out]):p.error('Supply package, archive and out')
        result=validate(a.package.resolve(),a.archive.resolve(),a.out.resolve(),a.stdlib_python or Path(sys.executable))
        print(json.dumps({k:result[k] for k in ['status','archive_regular_files','raw_states','balanced_reactions']}))


if __name__=='__main__':main()
