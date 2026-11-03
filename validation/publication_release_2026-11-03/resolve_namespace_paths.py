"""Resolve historical logical paths through recorded namespace mounts and Git."""
from collections import defaultdict
import argparse
import json
import os
from pathlib import Path
import subprocess

from redox.core.hashing import sha256_file
from redox.validation.source_inventory import bindings
from redox.workflows.archived_native_namespace import mount_command
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

with allocation_guard():
    host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
ROOT=Path(__file__).resolve().parents[2];FOLDER=Path(__file__).parent
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--source',type=Path,required=True)
parser.add_argument('--out',type=Path,required=True)
parser.add_argument('--summary',type=Path,required=True)
parser.add_argument('--discover-root',type=Path)
args=parser.parse_args();SOURCE=args.source.resolve();OUT=args.out.resolve()
if OUT.exists() or args.summary.exists():raise FileExistsError('Preserve previous resolution')
base=json.loads(SOURCE.read_text());aliases=defaultdict(list);namespace_sources={}
namespace_candidates={r['path']:r['sha256'] for r in base['files'] if Path(r['path']).name=='namespace.json'}
if args.discover_root:
    names=subprocess.check_output(['rg','--files','--hidden','--no-ignore',str(args.discover_root.resolve()),'-g','namespace.json'],text=True).splitlines()
    for name in names:namespace_candidates.setdefault(name,sha256_file(name))
for name,expected in namespace_candidates.items():
    p=Path(name)
    if sha256_file(p)!=expected:raise ValueError('Namespace changed')
    spec=json.loads(p.read_text())
    if not all(k in spec for k in ('original_campaign','attempt','original_raw_parent','raw_parent')):continue
    mount_command(spec)  # Validate declarations only; never execute the mount.
    namespace_sources[str(p)]=expected
    aliases[spec['original_campaign']].append(Path(spec['attempt']))
    if 'original_raw_root' in spec:
        relative=Path(spec['original_raw_root']).relative_to(spec['original_raw_parent'])
        aliases[spec['original_raw_root']].append(Path(spec['raw_parent'])/relative)
    else:
        aliases[spec['original_raw_parent']].append(Path(spec['raw_parent']))

git_path='scripts/validation/analyze_native_conformer_batches.py'
git_commit='b9657248f6fb0969beff5f8d17a4af6580a5fa8b'
snapshot=FOLDER/'historical_analyzer_source.py'
if not snapshot.exists():snapshot.write_bytes(subprocess.check_output(['git','show',git_commit+':'+git_path],cwd=ROOT))
expected_analyzer='b72bf5994b9542e72f186949afd825922cda194fec09ec9421ba758c854b4883'
if sha256_file(snapshot)!=expected_analyzer:raise ValueError('Git snapshot differs from original binding')
known=defaultdict(list);covered=set();files={(r['path'],r['sha256']):r for r in base['files']}
for row in base['files']:
    known[row['sha256']].append(Path(row['path']));covered.add((row['path'],row['sha256']))
for row in base['explicit_archive_relocations']:covered.add((row['original'],row['sha256']))
for p in (snapshot,ROOT/'scripts/validation/run_screening_polymer_control.py'):
    known[sha256_file(p)].append(p)
relocations=list(base['explicit_archive_relocations']);new_files=[];cache={};unresolved=[]
queue=[(r['path'],r['expected_sha256'],r['referenced_by']) for r in base['unresolved']]
visited=set()
def matches(path,digest):
    if not path.is_file():return False
    signature=(path.stat().st_size,path.stat().st_mtime_ns)
    key=str(path)
    if key not in cache or cache[key][0]!=signature:cache[key]=(signature,sha256_file(path))
    return cache[key][1]==digest

while queue:
    name,digest,parents=queue.pop();key=(name,digest)
    if key in covered or key in visited:continue
    visited.add(key);original=Path(name);candidates=[original]
    for prefix in (original,*original.parents):
        candidates.extend(target/original.relative_to(prefix) for target in aliases.get(str(prefix),()))
    candidates.extend(known[digest])
    found=next((p.resolve() for p in candidates if matches(p,digest)),None)
    if found is None:
        unresolved.append(dict(path=name,expected_sha256=digest,referenced_by=parents));continue
    covered.add(key);known[digest].append(found)
    if str(found)!=name:relocations.append(dict(original=name,path=str(found),sha256=digest))
    physical=(str(found),digest)
    if physical in files:continue
    row=dict(path=str(found),sha256=digest,bytes=found.stat().st_size,category='research_evidence_or_source_snapshot')
    files[physical]=row;new_files.append(row)
    if found.suffix=='.json':
        for child,h in bindings(json.loads(found.read_text())):
            child=Path(child);child=child if child.is_absolute() else ROOT/child
            queue.append((str(child),h,[str(found)]))

result=dict(base,status='declared_dependency_inventory_verified' if not unresolved else 'declared_dependency_inventory_incomplete',
    files=list(files.values()),explicit_archive_relocations=relocations,unresolved=unresolved,
    unique_content_bytes=sum({r['sha256']:r['bytes'] for r in files.values()}.values()),
    baseline_inventory=str(SOURCE),baseline_inventory_sha256=sha256_file(SOURCE),
    baseline_raw_files_rehashed=False,new_files_hash_checked=len(new_files),
    namespace_mounts_executed=False,namespace_sources_sha256=namespace_sources,
    recovered_git_source=dict(commit=git_commit,path=git_path,snapshot=str(snapshot),sha256=expected_analyzer))
result['resolution_source_sha256']={str(p):sha256_file(p) for p in (Path(__file__),Path(mount_command.__code__.co_filename),snapshot)}
OUT.write_text(json.dumps(result,indent=2)+'\n')
summary=dict(status=result['status'],inventory=str(OUT),inventory_sha256=sha256_file(OUT),
    verified_files=len(files),unresolved_bindings=len(unresolved),new_files_hash_checked=len(new_files),
    unique_content_bytes=result['unique_content_bytes'],namespace_mounts_executed=False,
    public_deposition_complete=False)
args.summary.write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary))
