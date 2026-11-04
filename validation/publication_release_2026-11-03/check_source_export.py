"""Extract committed source and replay key claims with original paths blocked."""
import json
import os
from pathlib import Path,PurePosixPath
import subprocess
import sys
import tarfile

from redox.core.hashing import sha256_file
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

with allocation_guard():
    host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
ROOT=Path(__file__).resolve().parents[2];FOLDER=Path(__file__).resolve().parent
record=FOLDER/'source_export.json';export=json.loads(record.read_text());archive=Path(export['archive'])
out=FOLDER/'source_replay.json'
extract=ROOT.parent/'redox_release_candidates/2026-11-03_publication_source_extracted'
if extract.exists() or out.exists():raise FileExistsError('Preserve previous extraction/replay')
if sha256_file(archive)!=export['archive_sha256']:raise ValueError('Source archive changed')
extract.mkdir();members=0;decoder=subprocess.Popen(['zstd','-dc',str(archive)],stdout=subprocess.PIPE)
try:
    with tarfile.open(fileobj=decoder.stdout,mode='r|') as stream:
        for item in stream:
            name=PurePosixPath(item.name)
            if name.is_absolute() or '..' in name.parts:raise ValueError('Unsafe source member')
            if item.isdir():continue
            if not item.isfile():raise ValueError('Unexpected non-file source member')
            target=extract.joinpath(*name.parts);target.parent.mkdir(parents=True,exist_ok=True)
            with stream.extractfile(item) as source,target.open('xb') as dest:
                for block in iter(lambda:source.read(2**20),b''):dest.write(block)
            members+=1
    decoder.stdout.close()
    if decoder.wait()!=0:raise ValueError('Source extraction failed')
except BaseException:
    decoder.stdout.close();decoder.terminate();decoder.wait();raise
blocked=[ROOT,ROOT.parent/'archive',ROOT.parent/'redox_candidate_continuations',
         ROOT.parent/'redox_counterion_controls',ROOT.parent/'redox_polymer_continuations',
         ROOT.parent/'redox_kernel_controls',ROOT.parent/'redox_motif_interactions']
code='''import json,pathlib,runpy,sys
root=pathlib.Path(sys.argv[1]);kind=sys.argv[2];blocked=[pathlib.Path(p).resolve() for p in sys.argv[3:]]
def deny(event,args):
 if event=='open' and args and isinstance(args[0],(str,bytes)):
  p=pathlib.Path(args[0]).resolve()
  if any(p==b or p.is_relative_to(b) for b in blocked):raise PermissionError('Original scientific paths blocked')
sys.addaudithook(deny)
if kind=='phenothiazine':
 sys.path.insert(0,str(root/'src'))
 pointer=json.loads((root/'manuscript/evidence/current.json').read_text())
 checker=runpy.run_path(str(root/'manuscript/check_phenothiazine_redox.py'))
 print(json.dumps(checker['audit']((root/pointer['manifest']).parent)))
else:
 path=root/('examples/finite_viologen_replay/replay.py' if kind=='finite_viologen' else 'examples/exact_viologen_corrected_replay/replay.py')
 sys.argv=[str(path)]+(['replay'] if kind=='corrected_viologen' else [])
 runpy.run_path(str(path),run_name='__main__')
'''
results={}
for kind in ('phenothiazine','finite_viologen','corrected_viologen'):
    run=subprocess.run([sys.executable,'-I','-c',code,str(extract),kind,*map(str,blocked)],
                       cwd='/',capture_output=True,text=True,check=True)
    results[kind]=json.loads(run.stdout)
report=dict(status='committed_source_three_replays_passed_with_original_paths_blocked',
    source_commit=export['source_commit'],archive_sha256=export['archive_sha256'],
    extracted_source=str(extract),extracted_files=members,replays=results,
    original_scientific_paths_blocked=True,quantum_calculations_repeated=False,
    all_manuscript_claims_replayed=False,public_deposition_complete=False,
    source_sha256={str(p):sha256_file(p) for p in (record,Path(__file__))})
out.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(dict(status=report['status'],files=members,replays=list(results))))
