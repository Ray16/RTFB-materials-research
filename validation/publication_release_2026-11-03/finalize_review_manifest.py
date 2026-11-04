"""Append a final local review manifest only after full raw byte verification."""
import json
import os
from pathlib import Path
import shutil

from redox.core.hashing import sha256_file
from redox.validation.release_gate import verified_archive_entry
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

ROOT=Path(__file__).resolve().parents[2]
FOLDER=Path(__file__).parent


def main():
    verified_path=FOLDER/'full_raw_archive_verified.json'
    if not verified_path.exists():raise RuntimeError('Full raw archive verification is still pending')
    owner=json.loads((FOLDER/'full_raw_archive_owner.json').read_text())
    if owner['stage']!='verified':raise RuntimeError('Archive owner has not closed the verification stage')
    previous=json.loads((FOLDER/'review_candidate_05.json').read_text());bundle=Path(previous['candidate'])
    manifest_path=bundle/'manifest_5f50c61_with_raw.json';report_path=FOLDER/'review_candidate_06.json'
    if manifest_path.exists() or report_path.exists():raise FileExistsError('Preserve earlier finalization')
    old_path=Path(previous['manifest'])
    if sha256_file(old_path)!=previous['manifest_sha256']:raise ValueError('Previous package manifest changed')
    base=json.loads(old_path.read_text());raw_path=bundle/'records/manifest.json'
    raw=json.loads(raw_path.read_text());verified=json.loads(verified_path.read_text())
    inventory=json.loads((FOLDER/'current_raw_archive_inventory.json').read_text())
    entry=verified_archive_entry(verified,raw,sha256_file(raw_path),inventory,bundle)
    if entry['file']!=str(Path(owner['archive']).relative_to(bundle)):
        raise ValueError('Verification does not belong to the declared archive owner')
    for name,digest in base['files_sha256'].items():
        if sha256_file(bundle/name)!=digest:raise ValueError('Previously checked package file changed: '+name)
    files=dict(base['files_sha256']);files[entry['file']]=entry['sha256']
    for name,digest in raw['support_files_sha256'].items():
        path=bundle/'records'/name
        if sha256_file(path)!=digest:raise ValueError('Raw archive support file changed')
        files[str(path.relative_to(bundle))]=digest
    files[str(raw_path.relative_to(bundle))]=verified['manifest_sha256']
    for source,name in [(verified_path,'raw_verification_d8855682.json'),
        (Path(inventory['inventory']),'raw_inventory_d8855682.json'),
        (Path(__file__),'finalize_review_source.py'),
        (ROOT/'src/redox/validation/release_gate.py','release_gate_source.py')]:
        target=bundle/'records'/name
        if target.exists():raise FileExistsError(target)
        shutil.copyfile(source,target)
        if sha256_file(source)!=sha256_file(target):raise ValueError('Finalization copy differs')
        files[str(target.relative_to(bundle))]=sha256_file(target)
    manifest=dict(base,status='local_current_review_package_with_complete_verified_declared_raw_archive',
        packages=[*base['packages'],entry],files_sha256=files,full_new_raw_archive_included=True,
        raw_archive_status_at_capture='all_archive_members_and_hashes_verified',
        raw_archive_verification='records/raw_verification_d8855682.json',
        prior_manifest_sha256=previous['manifest_sha256'],public_deposition_complete=False,
        numerical_or_experimental_qualification_added=False,submission_ready=False)
    manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
    report=dict(status=manifest['status'],candidate=str(bundle),manifest=str(manifest_path),
        manifest_sha256=sha256_file(manifest_path),source_commit=base['source_commit'],
        overleaf_commit=base['overleaf_commit'],article_pages=15,SI_pages=41,
        packages=len(manifest['packages']),package_bytes=sum((bundle/r['file']).stat().st_size for r in manifest['packages']),
        checked_files=len(files),raw_objects=entry['objects'],raw_archive_sha256=entry['sha256'],
        full_new_raw_archive_included=True,public_deposition_complete=False,data_license=None,
        submission_ready=False,producer_sha256=sha256_file(Path(__file__)))
    report_path.write_text(json.dumps(report,indent=2)+'\n')
    (bundle/'CURRENT.md').write_text('''# Current complete local review package

Use manifest_5f50c61_with_raw.json. The seven payloads include scientific source
commit 6629fcd1, the checked 15-page Article and 41-page SI, current submission
components, earlier scoped replay/control packages, and the full declared raw
archive. All 46,789 raw objects and the complete archive member set have passed
byte verification. This adds no physical, experimental or material qualification.

After extracting the raw archive, run `python -I verify.py` from its
publication_raw_evidence directory to verify its objects and logical mappings.
The original absolute paths describe provenance; the extracted verifier reads
local content-addressed objects. Runtime dependencies and third-party publications
remain explicitly external. The inventory covers declared dependency maps, not
implicit dependencies or a live census of every campaign.

Earlier manifests remain preserved. Public deposition, data reuse metadata and
final author declarations remain pending. This is not a journal submission.
''')
    print(json.dumps(report))


if __name__=='__main__':
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    main()
