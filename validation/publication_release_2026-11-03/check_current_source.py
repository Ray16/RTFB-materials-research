"""Build both scientific evidence and shared prose from a fresh in-memory export.

Durable outputs use explicit existing directories. The extracted checkout lives
only in a private tmpfs namespace and disappears when the check ends. Original
user data and networking are hidden; the installed runtime is reused.
"""
import argparse
import base64
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(2**20),b''):h.update(block)
    return h.hexdigest()


CHILD = r'''
from pathlib import Path, PurePosixPath
import base64, contextlib, io, json, os, re, runpy, subprocess, sys, tarfile, zipfile
home=Path(sys.argv[1])
assert sorted(p.name for p in home.iterdir())==['miniforge3'], 'Original user data visible'
root=Path('/tmp/work');root.mkdir();seen=set()
decoder=subprocess.Popen(['zstd','-dc','/tmp/source.tar.zst'],stdout=subprocess.PIPE)
with tarfile.open(fileobj=decoder.stdout,mode='r|') as stream:
 for item in stream:
  name=PurePosixPath(item.name)
  if name.is_absolute() or '..' in name.parts:raise ValueError('Unsafe source member')
  if item.isdir():continue
  if not item.isfile() or item.name in seen:raise ValueError('Unexpected source member')
  seen.add(item.name);target=root.joinpath(*name.parts);target.parent.mkdir(parents=True,exist_ok=True)
  with stream.extractfile(item) as source,target.open('xb') as out:
   for block in iter(lambda:source.read(2**20),b''):out.write(block)
  target.chmod(item.mode & 0o777)
decoder.stdout.close();assert decoder.wait()==0
(root/'logs').mkdir(exist_ok=True);(root/'logs/resource_locks').symlink_to('/tmp/shared_locks',target_is_directory=True)
os.chdir(root);sys.path[:0]=[str(root/'manuscript'),str(root/'src')]
with contextlib.redirect_stdout(io.StringIO()):
 runpy.run_path(str(root/'manuscript/build.py'),run_name='__main__')
 scientific=json.loads((root/'manuscript/build_check.json').read_text())
 sys.path.insert(0,str(root/'manuscript/overleaf/tools'))
 runpy.run_path(str(root/'manuscript/overleaf/tools/build_documents.py'),run_name='__main__')
shared=json.loads((root/'manuscript/overleaf/narrative_build_check.json').read_text())
from pypdf import PdfReader
normalize=lambda text:re.sub(r'\s+',' ',text).strip()
matched={}
for name,alias in [('main','manuscript'),('supporting_information','supporting_information')]:
 actual=PdfReader(Path(shared['build_directory'])/(name+'.pdf'))
 expected=PdfReader(root/'manuscript/review'/(alias+'.pdf'))
 matched[name]=(len(actual.pages)==len(expected.pages) and all(normalize(a.extract_text())==normalize(b.extract_text()) for a,b in zip(actual.pages,expected.pages)))
 assert matched[name], 'Compiled shared document text changed: '+name
 assert not shared['documents'][name]['overfull_boxes']
logs=io.BytesIO()
with zipfile.ZipFile(logs,'w',compression=zipfile.ZIP_DEFLATED) as output:
 for label,report in [('scientific',scientific),('shared',shared)]:
  output.writestr(label+'/build_check.json',json.dumps(report,indent=2))
  folder=Path(report.get('build_directory',report.get('build_logs')))
  for p in folder.iterdir():
   if p.is_file() and p.suffix in ('.log','.blg'):output.write(p,label+'/'+p.name)
print(json.dumps(dict(status='committed_scientific_and_shared_builds_passed_in_private_namespace',
 extracted_files=len(seen),scientific_pages={k:v['pages'] for k,v in scientific['documents'].items()},
 shared_pages={k:v['pages'] for k,v in shared['documents'].items()},shared_page_text_matches=matched,
 logs_base64=base64.b64encode(logs.getvalue()).decode())))
'''


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ('record','out','logs'):parser.add_argument('--'+key,type=Path,required=True)
    args=parser.parse_args()
    if args.out.exists() or args.logs.exists():raise FileExistsError('Preserve previous checks')
    record=json.loads(args.record.read_text());archive=Path(record['archive'])
    if sha(archive)!=record['archive_sha256']:raise ValueError('Source archive changed')
    home=ROOT.parent
    command=['bwrap','--die-with-parent','--unshare-net','--ro-bind','/','/','--tmpfs','/tmp',
        '--ro-bind',str(archive),'/tmp/source.tar.zst',
        '--bind',str(ROOT/'logs/resource_locks'),'/tmp/shared_locks',
        '--tmpfs',str(home),'--ro-bind',str(home/'miniforge3'),str(home/'miniforge3'),
        '--proc','/proc','--dev-bind','/dev','/dev','--chdir','/tmp','--',
        str(Path(sys.executable).resolve()),'-I','-c',CHILD,str(home)]
    env=dict(os.environ,CUDA_VISIBLE_DEVICES='',XDG_CACHE_HOME='/tmp/cache')
    env.pop('REDOX_MANUSCRIPT_BUILD_ROOT',None)
    run=subprocess.run(command,cwd='/',env=env,capture_output=True,text=True)
    if run.returncode:
        args.out.write_text(json.dumps(dict(status='isolated_source_build_failed_preserved',source_commit=record['source_commit'],returncode=run.returncode,stderr=run.stderr,stdout=run.stdout),indent=2)+'\n')
        print(run.stderr[-2000:]);raise SystemExit(run.returncode)
    result=json.loads(run.stdout);logs=base64.b64decode(result.pop('logs_base64'));args.logs.write_bytes(logs)
    if sha(archive)!=record['archive_sha256']:raise ValueError('Source archive changed during check')
    result.update(source_commit=record['source_commit'],archive_sha256=record['archive_sha256'],
        original_user_data_hidden=True,network_disabled=True,no_CUDA_devices=True,
        existing_runtime_reused=True,fresh_environment_installed=False,new_quantum_calculations=False,
        persistent_extracted_directory_created=False,logs=str(args.logs),logs_sha256=sha(args.logs),
        producer_sha256=sha(Path(__file__)))
    args.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))


if __name__=='__main__':main()
