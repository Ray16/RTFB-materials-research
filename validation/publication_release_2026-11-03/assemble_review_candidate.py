"""Collect current documents/source and existing checked replay packages locally."""
import json
import os
from pathlib import Path
import shutil

from redox.core.hashing import sha256_file
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

with allocation_guard():
    host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
ROOT=Path(__file__).resolve().parents[2];FOLDER=Path(__file__).parent
OUT=ROOT.parent/'redox_release_candidates/2026-11-03_jcim_review_candidate'
if OUT.exists():raise FileExistsError('Preserve earlier review candidates')
rows=[]
for key,path in (
    ('scientific_source',FOLDER/'source_export.json'),
    ('historical_reviews',ROOT/'validation/reviewer_collection_2026-10-31/summary.json'),
    ('electronic_controls',ROOT/'validation/jcim_completion_2026-11-02/component_archive_check.json'),
    ('manuscript_upload',ROOT/'manuscript/overleaf_export.json')):
    data=json.loads(path.read_text())
    rows.append(dict(id=key,source=data['archive'],sha256=data['archive_sha256'],record=str(path)))
index=json.loads((ROOT/'results/INDEX.json').read_text())['entries']['corrected_primary_CPU_example']
rows.append(dict(id='corrected_molecular_replay',source=index['archive'],sha256=index['archive_sha256'],record=str(ROOT/index['current'])))
deck=json.loads((ROOT/'slides/build_check.json').read_text())
rows.append(dict(id='slides',source=str(ROOT/'slides/redox_workflow.pdf'),sha256=deck['pdf_sha256'],record=str(ROOT/'slides/build_check.json')))
required=sum(Path(r['source']).stat().st_size for r in rows)
# Small bundle only: leave the existing live-calculation promises plus reserve.
headroom=json.loads((FOLDER/'storage_headroom.json').read_text())
if shutil.disk_usage(OUT.parent).free<required+headroom['remaining_existing_budget_bytes']+32*2**30:
    raise RuntimeError('Current calculation/storage margin unavailable')
OUT.mkdir();(OUT/'packages').mkdir();(OUT/'records').mkdir()
for row in rows:
    source=Path(row['source'])
    if sha256_file(source)!=row['sha256']:raise ValueError('Checked package changed')
    suffix='.tar.zst' if source.name.endswith('.tar.zst') else source.suffix
    target=OUT/'packages'/(row['id']+suffix);shutil.copyfile(source,target)
    if sha256_file(target)!=row['sha256']:raise ValueError('Copied package differs')
    row['file']=str(target.relative_to(OUT))
    shutil.copyfile(row['record'],OUT/'records'/(row['id']+'.json'))
for name in ('source_replay.json','raw_seed_manifest.json','storage_headroom.json'):
    shutil.copyfile(FOLDER/name,OUT/'records'/name)
(OUT/'README.md').write_text('''# JCIM local review candidate — 3 November 2026

The manuscript upload contains the checked 21-page article and 36-page SI.
The separate slide PDF has 23 pages. The scientific-source archive is Git
commit 7d19efe8; unrelated uncommitted work is excluded. Existing historical
review and corrected-molecular packages retain their own dated method scopes.
Package state counts overlap and must not be added as unique calculations.

For quick numerical review, extract the scientific-source archive and run:

    python -I examples/finite_viologen_replay/replay.py
    python -I examples/exact_viologen_corrected_replay/replay.py replay

The second command requires NumPy under its supplied environment lock. The
source also contains manuscript/check_phenothiazine_redox.py; its audit(folder)
function replays the frozen phenothiazine value with PYTHONPATH=src and the
evidence folder selected by manuscript/evidence/current.json. All three paths
were tested from an extracted source copy with original scientific paths blocked.

The electronic-controls package includes actual SCF checkpoints and its own
CPU replay instructions. Historical packages include earlier native caches,
PMDI/TEMPO/fragment and benchmark evidence; their individual README files define
which calculations and diagnostics their replays cover.

This candidate does not yet contain the full newly inventoried raw dependency
archive. Its path recovery and storage assessment are separate ongoing work;
the supplied source/replay evidence is not a full-raw-release certificate.
Author details, reuse licensing and public deposition remain outstanding.
No journal submission or external publication has been performed.
''')
manifest=dict(status='local_JCIM_review_candidate_full_new_raw_archive_pending',
    source_commit=json.loads((FOLDER/'source_export.json').read_text())['source_commit'],
    packages=rows,files_sha256={str(p.relative_to(OUT)):sha256_file(p) for p in OUT.rglob('*') if p.is_file()},
    full_new_raw_archive_included=False,public_deposition_complete=False,data_license=None,
    author_metadata_handled_by_user=True,submission_ready=False,
    producer_sha256=sha256_file(__file__))
(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
report=dict(status=manifest['status'],candidate=str(OUT),manifest_sha256=sha256_file(OUT/'manifest.json'),
    packages=len(rows),package_bytes=required,checked_files=len(manifest['files_sha256']),
    full_new_raw_archive_included=False,public_deposition_complete=False,submission_ready=False,
    producer_sha256=sha256_file(__file__))
(FOLDER/'review_candidate.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report))
