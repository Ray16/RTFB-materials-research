"""Inventory declared hash-map dependencies of six checked contact outcomes.

Read-only: verify bytes, classify dependencies and estimate an archival payload.
This does not copy bulk data, assign reuse rights or deposit an archive.
"""
import json
import os
from pathlib import Path
import subprocess

from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
from redox.core.hashing import sha256_file

ROOT=Path(__file__).resolve().parents[2]
ARCHIVE=ROOT.parent/'archive/0_redox'
OUT=Path(__file__).parent/'inventory.json'
if OUT.exists():raise FileExistsError('Preserve prior inventories')
summary_path=ROOT/'validation/counterion_completion_2026-11-03/current.json'
summary=json.loads(summary_path.read_text())
primary_path=ROOT/'validation/ionic_cycle_completion_2026-11-02/primary_cohort_scalar_inputs.json'
primary=json.loads(primary_path.read_text())
tracked=set(subprocess.check_output(['git','ls-files'],cwd=ROOT,text=True).splitlines())
expected=set();roles={};queue=[];seen=set();relocations={};parents={};unresolved={}
def bind(path,digest,role,parent=None):
    original=Path(path).resolve();resolved=original
    # Pre-refactor logical paths are retained in immutable certificates. An
    # explicit archive relocation is accepted only for the exact recorded bytes.
    if original.is_relative_to(ROOT) and (not original.is_file() or sha256_file(original)!=digest):
        candidate=ARCHIVE/original.relative_to(ROOT)
        if candidate.is_file() and sha256_file(candidate)==digest:
            resolved=candidate.resolve();relocations[(str(original),digest)]=dict(original=str(original),path=str(resolved),sha256=digest)
    name=str(resolved)
    key=(name,digest)
    expected.add(key);roles.setdefault(key,set()).add(role)
    if parent is not None:parents.setdefault(key,set()).add(parent)
    if key not in seen and Path(name).suffix=='.json':queue.append(key)

bind(summary_path,sha256_file(summary_path),'contact_collection')
for row in summary['records']:
    bind(row['review'],row['review_sha256'],'checked_contact_review')
    bind(row['ledger'],row['ledger_sha256'],'association_cycle')
bind(primary_path,sha256_file(primary_path),'matched_primary_inputs')
for key in ('graft_ox2','chloride'):
    for path,digest in primary['components'][key]['sources'].items():bind(path,digest,'matched_primary_parent')

while queue:
    key=queue.pop();name,digest=key
    if key in seen:continue
    seen.add(key);path=Path(name)
    actual=sha256_file(path) if path.is_file() else None
    if actual!=digest:
        unresolved[key]=dict(expected_sha256=digest,actual_sha256=actual,
            reason='missing' if actual is None else 'content_changed',children_traversed=False)
        continue
    value=json.loads(path.read_text())
    if not isinstance(value,dict):continue
    for field in ('source_sha256','evidence_sha256','producer_sha256'):
        bindings=value.get(field)
        if not isinstance(bindings,dict):continue
        for source,digest in bindings.items():
            if not isinstance(digest,str) or len(digest)!=64 or any(c not in '0123456789abcdef' for c in digest):continue
            source=Path(source)
            if not source.is_absolute():source=ROOT/source
            bind(source,digest,'nested_'+field,parent=name)

files=[];totals={};unique={}
for name,digest in sorted(expected):
    key=(name,digest)
    path=Path(name)
    actual=sha256_file(path) if path.is_file() else None
    if actual!=digest:
        unresolved[key]=dict(expected_sha256=digest,actual_sha256=actual,
            reason='missing' if actual is None else 'content_changed',children_traversed=key in seen and key not in unresolved)
        continue
    relative=str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else None
    if relative in tracked:category='versioned_project_artifact'
    elif any(part in name for part in ('/miniforge3/','/site-packages/','/environments/')):category='runtime_dependency'
    elif '/literature/' in name and path.suffix in ('.pdf','.xml','.html'):category='third_party_publication'
    else:category='unversioned_research_evidence_or_source_snapshot'
    size=path.stat().st_size;totals.setdefault(category,dict(files=0,bytes=0));totals[category]['files']+=1;totals[category]['bytes']+=size
    unique.setdefault((category,digest),size)
    files.append(dict(path=name,sha256=digest,bytes=size,category=category,roles=sorted(roles[key])))
report=dict(status='declared_dependency_inventory_incomplete' if unresolved else 'declared_dependency_inventory_verified',source_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
    input_collection=str(summary_path),checked_contact_outcomes=6,declared_scope='Six primary contact outcomes and their matched graft/chloride parents',
    files=files,totals=totals,explicit_archive_relocations=list(relocations.values()),
    unresolved=[dict(path=key[0],**row,referenced_by=sorted(parents.get(key,()))) for key,row in sorted(unresolved.items())],
    traversal_fields=['source_sha256','evidence_sha256','producer_sha256'],
    traversal_scope='Top-level path-to-SHA256 mappings in verified JSON only; unresolved JSON stops its dependency branch; aggregate provenance may include other chemistry',
    all_implicit_dependencies_inventoried=False,
    unique_content_bytes_by_category={c:sum(n for (k,h),n in unique.items() if k==c) for c in totals},
    raw_files_copied=False,public_deposition_complete=False,data_license=None,
    third_party_material_relicensed=False,full_project_coverage_inferred=False,
    quantum_validation_repeated=False,producer_sha256=sha256_file(__file__))
OUT.write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(dict(status=report['status'],verified_files=len(files),unresolved=len(unresolved),totals=totals,unique_content_bytes_by_category=report['unique_content_bytes_by_category'])))
