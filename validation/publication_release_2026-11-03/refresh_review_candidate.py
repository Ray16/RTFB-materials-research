"""Append a checked current review manifest inside the existing release folder."""
import hashlib
import json
import os
from pathlib import Path
import shutil

from redox.core.hashing import sha256_file
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

ROOT=Path(__file__).resolve().parents[2]
FOLDER=Path(__file__).parent
with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
BUNDLE=ROOT/'results/redox_release_candidates/2026-11-03_jcim_review_candidate_submission'
report_path=FOLDER/'review_candidate_05.json';manifest_path=BUNDLE/'manifest_5f50c61.json'
if report_path.exists() or manifest_path.exists():raise FileExistsError('Preserve earlier review manifests')
old=json.loads((BUNDLE/'manifest.json').read_text())
source=json.loads((FOLDER/'current_source_export.json').read_text())
checked=json.loads((FOLDER/'current_source_validation.json').read_text())
if (checked['status']!='committed_scientific_and_shared_builds_passed_in_private_namespace'
        or source['source_commit']!=checked['source_commit']
        or source['archive_sha256']!=checked['archive_sha256']):
    raise ValueError('Current source and isolated validation must match')
submission=json.loads((BUNDLE/'records/submission_files_5f50c61.json').read_text())
manuscript=json.loads((BUNDLE/'records/manuscript_5f50c61.json').read_text())
rows=[]
for row in old['packages']:
    current=dict(id=row['id'],file=row['file'],sha256=row['sha256'])
    if row['id']=='scientific_source':
        current.update(file=str(Path(source['archive']).relative_to(BUNDLE)),sha256=source['archive_sha256'],source_commit=source['source_commit'])
    elif row['id']=='manuscript_upload':
        current.update(file=str(Path(manuscript['archive']).relative_to(BUNDLE)),sha256=manuscript['archive_sha256'],overleaf_commit=manuscript['overleaf_commit'])
    else:current['scope']='Existing checked package; retains its recorded method and date scope'
    if sha256_file(BUNDLE/current['file'])!=current['sha256']:raise ValueError('A checked payload changed')
    rows.append(current)
files={r['file']:r['sha256'] for r in rows}
for row in submission['files']:
    if sha256_file(BUNDLE/row['file'])!=row['sha256']:raise ValueError('Submission file changed')
    files[row['file']]=row['sha256']
for name in ('scientific_source_6629fcd1.json','manuscript_5f50c61.json','submission_files_5f50c61.json',
             'source_validation_6629fcd1.json','check_source_6629fcd1.py','manuscript_snapshot_5f50c61.json','validation_logs_6629fcd1.zip'):
    files['records/'+name]=sha256_file(BUNDLE/'records'/name)
for row in old['packages']:
    if row['id'] not in ('scientific_source','manuscript_upload'):
        name='records/'+row['id']+'.json';files[name]=sha256_file(BUNDLE/name)
producer=BUNDLE/'records/current_review_packaging_source.py'
if producer.exists():raise FileExistsError(producer)
shutil.copyfile(__file__,producer);files[str(producer.relative_to(BUNDLE))]=sha256_file(producer)
manifest=dict(status='current_shared_manuscript_and_committed_source_review_package_checked',
    source_commit=source['source_commit'],overleaf_commit=manuscript['overleaf_commit'],
    article_pages=15,SI_pages=41,packages=rows,submission_files=submission['files'],files_sha256=files,
    committed_scientific_and_shared_builds_passed=True,shared_compiled_page_text_matches=True,
    full_new_raw_archive_included=False,raw_archive_status_at_capture='creation_or_verification_pending',
    raw_archive_owner_record='validation/publication_release_2026-11-03/full_raw_archive_owner.json',
    raw_archive_is_required_to_be_verified_before_inclusion=True,
    public_deposition_complete=False,data_license=None,author_metadata_handled_by_user=True,submission_ready=False)
manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
for name,h in files.items():
    if sha256_file(BUNDLE/name)!=h:raise ValueError('Current review package changed during verification')
(BUNDLE/'CURRENT.md').write_text('''# Current review files

Use manifest_5f50c61.json for the current 15-page Article, 41-page SI, draft cover
letter, current TOC TIFF, supporting tables, self-contained manuscript source,
and scientific source commit 6629fcd1. The source archive matches every Git blob;
its scientific and shared manuscript builds pass with original user data hidden
and networking disabled, using the existing checked runtime. Compiler logs are
included. Replays and older control/slide packages retain their individual scopes.

The original README.md and manifest.json describe an earlier preserved package.
Their payloads and records have not been replaced. Package state counts overlap.

The full raw-data archive is a separate in-progress artifact. Do not count it as
complete until full_raw_archive_verified.json and the corresponding updated
review manifest explicitly confirm every member. Public deposition, data-license
metadata and final author declarations remain pending. No submission is claimed.
''')
report=dict(status=manifest['status'],candidate=str(BUNDLE),manifest=str(manifest_path),
    manifest_sha256=sha256_file(manifest_path),source_commit=source['source_commit'],
    overleaf_commit=manuscript['overleaf_commit'],article_pages=15,SI_pages=41,
    packages=len(rows),package_bytes=sum((BUNDLE/r['file']).stat().st_size for r in rows),
    checked_files=len(files),submission_files=len(submission['files']),
    committed_scientific_and_shared_builds_passed=True,full_new_raw_archive_included=False,
    public_deposition_complete=False,submission_ready=False,producer_sha256=sha256_file(__file__))
report_path.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
