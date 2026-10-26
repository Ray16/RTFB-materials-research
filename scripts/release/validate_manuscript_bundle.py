#!/usr/bin/env python
"""Extract a reviewer bundle and replay its scoped cached results in isolation."""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path,PurePosixPath
import shutil
import subprocess
import sys
import tarfile
import time

ROOT=Path(__file__).resolve().parents[2]


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()


def safe_extract(stream,destination,expected_root=None):
    files=set();roots=set()
    for item in stream:
        path=PurePosixPath(item.name);parts=path.parts
        if (not parts or path.is_absolute() or '..' in parts or not (item.isdir() or item.isfile())
                or (expected_root is not None and parts[0]!=expected_root)):
            raise ValueError('Unsafe or unexpected archive member: '+item.name)
        roots.add(parts[0]);target=destination.joinpath(*parts)
        if not target.resolve().is_relative_to(destination.resolve()):raise ValueError('Archive path escapes destination')
        if item.isdir():target.mkdir(parents=True,exist_ok=True);continue
        if item.name in files or item.size>1024**3:raise ValueError('Duplicate or oversized archive member')
        target.parent.mkdir(parents=True,exist_ok=True)
        with stream.extractfile(item) as incoming,target.open('xb') as outgoing:shutil.copyfileobj(incoming,outgoing)
        files.add(item.name)
    return files,roots


def extract_zstd(archive,out,expected_root=None):
    proc=subprocess.Popen(['zstd','-dc',str(archive)],stdout=subprocess.PIPE)
    try:
        with tarfile.open(fileobj=proc.stdout,mode='r|') as stream:files,roots=safe_extract(stream,out,expected_root)
        proc.stdout.close()
        if proc.wait()!=0:raise ValueError('Compressed archive could not be decoded')
        return files,roots
    finally:
        if proc.poll() is None:proc.terminate();proc.wait()


BOOTSTRAP='''import os,runpy,sys\nforbidden=os.path.realpath(sys.argv[1])\ndef guard(event,args):\n    if event=="open" and isinstance(args[0],(str,bytes)):\n        path=os.path.realpath(os.fsdecode(args[0]))\n        if path==forbidden or path.startswith(forbidden+os.sep):\n            raise RuntimeError("Original project read blocked: "+path)\nsys.addaudithook(guard)\nsys.argv=sys.argv[2:]\nrunpy.run_path(sys.argv[0],run_name="__main__")\n'''


def validate(archive,bundle_name,out,python):
    archive=archive.resolve();out=out.resolve();python=python.resolve()
    if out.is_relative_to(ROOT):raise ValueError('Reviewer extraction must be outside the original checkout')
    if out.exists():raise FileExistsError('Preserve previous replay results')
    if not bundle_name or PurePosixPath(bundle_name).name!=bundle_name:raise ValueError('One bundle directory name required')
    out.mkdir(parents=True);payload=out/'payload';payload.mkdir()
    files,_=extract_zstd(archive,payload,bundle_name);bundle=payload/bundle_name
    manifest=json.loads((bundle/'manifest.json').read_text())
    if (manifest['status']!='local_manuscript_reviewer_bundle_prepared_not_public'
            or manifest['public_archive_deposited'] or manifest['public_DOI_created']
            or manifest['experimental_accuracy_inferred'] or manifest['conformer_convergence_inferred']):
        raise ValueError('Bundle scientific/publication scope changed')
    expected={bundle_name+'/'+p for p in manifest['files_sha256']}|{bundle_name+'/manifest.json'}
    if files!=expected:raise ValueError('Bundle file inventory differs')
    for relative,value in manifest['files_sha256'].items():
        if sha(bundle/relative)!=value:raise ValueError('Bundle file changed: '+relative)
    project=out/'project';project.mkdir()
    with tarfile.open(bundle/'project-source.tar',mode='r:') as stream:safe_extract(stream,project)
    pointer=json.loads((project/'manuscript/evidence/current.json').read_text())
    paper_manifest=json.loads((project/pointer['manifest']).read_text())
    for name,value in paper_manifest['snapshot_sha256'].items():
        if sha(project/name)!=value:raise ValueError('Relocated manuscript evidence differs: '+name)
    data=out/'data';data.mkdir();packages={}
    for item in manifest['packages']:
        label=item['id'];target=data/label;target.mkdir()
        archive_path=bundle/item['archive']
        if sha(archive_path)!=item['archive_sha256']:raise ValueError('Nested archive hash differs')
        _,roots=extract_zstd(archive_path,target)
        if len(roots)!=1:raise ValueError('Nested archive needs exactly one package root')
        packages[label]=target/next(iter(roots))
    if set(packages)!={'primary_review','first_site','numerical_controls','neighbor_cycle'}:
        raise ValueError('Expected four distinct reviewer packages')
    replay=out/'replay';replay.mkdir();env={'PATH':'/usr/bin:/bin','LC_ALL':'C.UTF-8'}
    for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[key]='1'
    jobs=[('primary',packages['primary_review']/'code/scripts/release/reproduce_review_package.py',[str(packages['primary_review']),'--out',str(replay/'primary.json')]),
        ('first_site',packages['first_site']/'reproduce.py',['--data',str(packages['first_site']),'--out',str(replay/'first_site.json')]),
        ('controls',packages['numerical_controls']/'example/reproduce.py',['--package',str(packages['numerical_controls']/'example'),'--electronic-only','--out',str(replay/'controls.json')]),
        ('neighbor',packages['neighbor_cycle']/'reproduce.py',['--data',str(packages['neighbor_cycle']),'--out',str(replay/'neighbor.json')]),
        ('manuscript',project/'manuscript/check_claims.py',[])]
    for name,script,args in jobs:
        command=[str(python),'-B','-I','-S','-c',BOOTSTRAP,str(ROOT),str(script),*args]
        completed=subprocess.run(command,cwd=out,env=env,capture_output=True,text=True,timeout=300)
        (replay/(name+'.log')).write_text(completed.stdout+completed.stderr)
        if completed.returncode:raise ValueError('Relocated replay failed: '+name)
    neighbor=json.loads((replay/'neighbor.json').read_text());first=json.loads((replay/'first_site.json').read_text());controls=json.loads((replay/'controls.json').read_text())
    if (neighbor['states']!=14 or neighbor['balanced_reactions']!=14 or abs(neighbor['by_cutoff']['50']['deltaDeltaG_kJ_mol']+1.8616197516972124)>1e-8
            or first['states']!=13 or controls['quantum_calculations_run'] or controls['experimental_accuracy_inferred']):
        raise ValueError('Relocated output coverage or scope differs')
    report=dict(status='consolidated_bundle_relocated_stdlib_replays_checked',checked=time.time(),
        archive=str(archive),archive_sha256=sha(archive),bundle_manifest_sha256=sha(bundle/'manifest.json'),
        source_commit=manifest['source_commit'],nested_packages=sorted(packages),replayed_stages=[name for name,_,_ in jobs],
        neighbor_states=neighbor['states'],neighbor_balanced_reactions=neighbor['balanced_reactions'],
        neighbor_deltaDeltaG_kJ_mol=neighbor['by_cutoff']['50']['deltaDeltaG_kJ_mol'],
        original_checkout_reads_blocked=True,scientific_packages_required=False,quantum_calculations_repeated=False,
        full_native_reintegration_inferred=False,public_archive_deposited=False,experimental_accuracy_inferred=False,
        validator_sha256=sha(__file__))
    (out/'validation.json').write_text(json.dumps(report,indent=2)+'\n');return report


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--archive',type=Path,required=True);p.add_argument('--bundle-name',required=True);p.add_argument('--out',type=Path,required=True);p.add_argument('--python',type=Path,default=Path('/usr/bin/python3'));a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    result=validate(a.archive,a.bundle_name,a.out,a.python);print(json.dumps({k:result[k] for k in ('status','source_commit','nested_packages','replayed_stages')}))
