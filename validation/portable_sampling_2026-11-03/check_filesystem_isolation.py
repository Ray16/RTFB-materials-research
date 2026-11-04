"""Build with original scientific directories hidden at the filesystem level."""
import json
import os
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[2];FOLDER=Path(__file__).parent
source=ROOT.parent/'redox_release_candidates/2026-11-03_publication_source_extracted'
out=FOLDER/'filesystem_isolation_02.json'
if out.exists():raise FileExistsError('Preserve prior isolation evidence')
blocked=[ROOT,ROOT.parent/'archive',*[ROOT.parent/name for name in
    ('redox_candidate_continuations','redox_kernel_controls','redox_polymer_continuations',
     'redox_counterion_controls','redox_motif_interactions')]]
aliases=[]
for path in list(blocked):
    alias=Path('/homes/rzhu')/path.name
    if alias.is_dir():
        aliases.append(alias)
        if alias.resolve()!=path.resolve():blocked.append(alias)
# Share only the CPU allocation lock; never an original scientific data file.
locks=ROOT/'logs/resource_locks';target=source/'logs/resource_locks'
target.mkdir(parents=True,exist_ok=True)
argv=['bwrap','--die-with-parent','--unshare-net','--ro-bind','/','/',
      '--bind',str(source),str(source),'--bind',str(locks),str(target)]
for path in blocked:argv.extend(['--tmpfs',str(path)])
argv.extend(['--tmpfs','/tmp','--proc','/proc','--dev-bind','/dev','/dev','--chdir',str(source),'--'])
code='''import pathlib,runpy,sys
root=pathlib.Path(sys.argv[1]);blocked=[pathlib.Path(x) for x in sys.argv[2:]]
if any(p.is_dir() and list(p.iterdir()) for p in blocked):raise RuntimeError('Original data directory is visible')
sys.path[:0]=[str(root/'manuscript'),str(root/'src')]
runpy.run_path(str(root/'manuscript/build.py'),run_name='__main__')
'''
argv.extend([sys.executable,'-I','-c',code,str(source),*map(str,[*blocked,*aliases])])
env=dict(os.environ,CUDA_VISIBLE_DEVICES='',XDG_CACHE_HOME=str(source/'calcs/isolated_cache'))
run=subprocess.run(argv,cwd='/',env=env,capture_output=True,text=True)
(FOLDER/'filesystem_isolation_02.stdout.log').write_text(run.stdout)
(FOLDER/'filesystem_isolation_02.stderr.log').write_text(run.stderr)
report=dict(status='filesystem_isolated_full_manuscript_build_passed' if run.returncode==0 else 'filesystem_isolated_build_failed_preserved',
    returncode=run.returncode,source=str(source),hidden_original_directories=list(map(str,blocked)),
    original_aliases_checked=list(map(str,aliases)),
    filesystem_isolation_includes_native_libraries_and_child_processes=True,
    network_namespace_isolated=True,CUDA_VISIBLE_DEVICES='',original_CPU_lock_shared_only=True,
    new_quantum_validation=False,argv=argv)
if run.returncode==0:
    built=json.loads((source/'manuscript/build_check.json').read_text())
    report['documents']={k:v['pages'] for k,v in built['documents'].items()}
out.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(dict(status=report['status'],returncode=run.returncode,error_tail=run.stderr[-1600:])))
raise SystemExit(run.returncode)
