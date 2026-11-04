"""Verify the committed source export by rebasing and hiding original user data."""
import hashlib
import json
import os
from pathlib import Path,PurePosixPath
import subprocess
import sys
import tarfile

ROOT=Path(__file__).resolve().parents[2];FOLDER=Path(__file__).parent
export_path=FOLDER/'source_export.json';export=json.loads(export_path.read_text())
archive=Path(export['archive']);source=ROOT.parent/('redox_release_candidates/2026-11-03_source_'+export['source_commit'][:8]+'_extracted')
out=FOLDER/'committed_source_validation.json'
if source.exists() or out.exists():raise FileExistsError('Preserve prior committed-source checks')
def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(2**20),b''):h.update(block)
    return h.hexdigest()
if digest(archive)!=export['archive_sha256']:raise ValueError('Committed source archive changed')
source.mkdir();decoder=subprocess.Popen(['zstd','-dc',str(archive)],stdout=subprocess.PIPE);count=0
try:
    with tarfile.open(fileobj=decoder.stdout,mode='r|') as stream:
        for member in stream:
            name=PurePosixPath(member.name)
            if name.is_absolute() or '..' in name.parts:raise ValueError('Unsafe archive path')
            if member.isdir():continue
            if not member.isfile():raise ValueError('Unexpected source member')
            target=source.joinpath(*name.parts);target.parent.mkdir(parents=True,exist_ok=True)
            with stream.extractfile(member) as input_file,target.open('xb') as dest:
                for block in iter(lambda:input_file.read(2**20),b''):dest.write(block)
            count+=1
    decoder.stdout.close()
    if decoder.wait()!=0:raise ValueError('Source extraction failed')
except BaseException:
    decoder.stdout.close();decoder.terminate();decoder.wait();raise
(source/'logs/resource_locks').mkdir(parents=True,exist_ok=True)
home=ROOT.parent
argv=['bwrap','--die-with-parent','--unshare-net','--ro-bind','/','/','--tmpfs','/tmp',
      '--bind',str(source),'/tmp/work','--bind',str(ROOT/'logs/resource_locks'),'/tmp/work/logs/resource_locks',
      '--tmpfs',str(home),'--ro-bind',str(home/'miniforge3'),str(home/'miniforge3'),
      '--proc','/proc','--dev-bind','/dev','/dev','--chdir','/tmp/work','--']
code='''from pathlib import Path
import runpy,sys
if sorted(p.name for p in Path(sys.argv[1]).iterdir())!=['miniforge3']:raise RuntimeError('Original user data visible')
sys.path[:0]=['/tmp/work/manuscript','/tmp/work/src']
runpy.run_path('/tmp/work/manuscript/build.py',run_name='__main__')
'''
argv += [str(Path(sys.executable).resolve()),'-I','-c',code,str(home)]
run=subprocess.run(argv,cwd='/',env=dict(os.environ,CUDA_VISIBLE_DEVICES='',XDG_CACHE_HOME='/tmp/work/calcs/isolated_cache'),
                   capture_output=True,text=True)
(FOLDER/'committed_source_validation.stdout.log').write_text(run.stdout)
(FOLDER/'committed_source_validation.stderr.log').write_text(run.stderr)
report=dict(status='committed_source_full_validation_build_passed' if run.returncode==0 else 'committed_source_validation_failed_preserved',
    returncode=run.returncode,source_commit=export['source_commit'],archive_sha256=export['archive_sha256'],
    extracted_source=str(source),extracted_files=count,namespace_source='/tmp/work',
    original_user_data_hidden=True,network_disabled=True,CUDA_VISIBLE_DEVICES='',
    existing_checked_runtime_reused=True,new_quantum_validation=False,
    argv=argv,source_sha256={str(p):digest(p) for p in (export_path,Path(__file__))})
if run.returncode==0:
    built=json.loads((source/'manuscript/build_check.json').read_text())
    report['documents']={k:v['pages'] for k,v in built['documents'].items()}
out.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(dict(status=report['status'],returncode=run.returncode,error_tail=run.stderr[-1400:])))
raise SystemExit(run.returncode)
