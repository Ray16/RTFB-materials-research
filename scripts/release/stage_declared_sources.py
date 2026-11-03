"""Stage a fresh local raw-evidence package from the resolved source inventory."""
import argparse
import json
import os
from pathlib import Path
import shutil

from redox.core.hashing import sha256_file
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def stage(inventory, content_resolution, metadata_resolution, output, *, resume=False):
    if output.exists() and not resume:
        raise FileExistsError('Use a fresh release directory')
    if resume and (not output.is_dir() or {p.name for p in output.iterdir()} != {'objects'}):
        raise ValueError('Resume only an unfinished objects-only staging directory')
    inv = json.loads(inventory.read_text())
    content = json.loads(content_resolution.read_text())
    metadata = json.loads(metadata_resolution.read_text())
    if (sha256_file(inventory) != content['inventory_sha256']
            or sha256_file(content_resolution) != metadata['previous_resolution_sha256']
            or metadata['unresolved']):
        raise ValueError('Exact resolved inventory chain required')
    mappings = list(inv['files']) + [dict(path=r['path'], sha256=r['sha256'],
        bytes=r['bytes'], category='unversioned_research_evidence_or_source_snapshot')
        for r in metadata['resolved']]
    excluded = {'runtime_dependency', 'third_party_publication'}
    external = [r for r in mappings if r['category'] in excluded]
    included = [r for r in mappings if r['category'] not in excluded]
    objects = {r['sha256']: r['bytes'] for r in included}
    remaining = sum(size for digest,size in objects.items() if not (output/'objects'/digest).exists())
    if shutil.disk_usage(output.parent).free < remaining + 32 * 2**30:
        raise RuntimeError('Keep the 32-GiB shared-storage margin')
    output.mkdir(exist_ok=resume); (output / 'objects').mkdir(exist_ok=resume)
    for row in included:
        source = Path(row['path']); target = output / 'objects' / row['sha256']
        if source.stat().st_size != row['bytes'] or sha256_file(source) != row['sha256']:
            raise ValueError('Inventoried source changed: ' + str(source))
        if not target.exists():
            shutil.copyfile(source, target)
        if target.is_symlink() or sha256_file(target) != row['sha256']:
            raise ValueError('Copied evidence differs')
    all_relocations = [*inv['explicit_archive_relocations'], *content['resolved'], *metadata['resolved']]
    external_hashes = {r['sha256'] for r in external}
    if any(row['sha256'] not in objects and row['sha256'] not in external_hashes for row in all_relocations):
        raise ValueError('Relocation has no inventoried dependency')
    relocations = [r for r in all_relocations if r['sha256'] in objects]
    external_relocations = [r for r in all_relocations if r['sha256'] not in objects]
    shutil.copyfile(Path(__file__).with_name('verify_source_objects.py'), output / 'verify.py')
    (output / 'README.md').write_text('''# Local raw-evidence package

Run `python -I verify.py` from this directory, or pass the package directory
explicitly. The standard-library verifier checks every content-addressed raw
object and its file/relocation mapping without original scientific paths.
The inventory includes six chloride outcomes and parents, plus broader
historical molecular/PMDI/base evidence reached through aggregate provenance.

This is a local review package, not a public deposit or a numerical replay.
SCF, Hessian, spin and reaction calculations are not rerun by the verifier.
The original records remain byte-identical; path mappings are external to them.
Runtime dependencies and third-party publications are listed in the manifest
but not redistributed or verified here. Technical inventory categories confer
no license. Existing code rights remain applicable; data reuse terms and public
deposition remain pending. No experimental or ensemble qualification is added.
''')
    data = dict(status='local_declared_raw_evidence_staged_not_deposited',
        objects=objects, files=included, relocations=relocations, external_dependencies=external,
        external_relocations=external_relocations,
        inventory_scope=inv['traversal_scope'], implicit_dependencies_covered=False,
        public_deposition_complete=False, data_license=None,
        source_sha256={str(p): sha256_file(p) for p in
                       (inventory, content_resolution, metadata_resolution, Path(__file__))})
    (output / 'manifest.json').write_text(json.dumps(data, indent=2) + '\n')
    return dict(objects=len(objects), bytes=sum(objects.values()), external_dependencies=len(external))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('inventory', 'content-resolution', 'metadata-resolution', 'out'):
        parser.add_argument('--' + key, type=Path, required=True)
    parser.add_argument('--resume', action='store_true', help='Rehash an unfinished objects-only directory')
    args = parser.parse_args()
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4); os.sched_setaffinity(0, cpu_allocation(1))
    print(json.dumps(stage(args.inventory.resolve(), args.content_resolution.resolve(),
                           args.metadata_resolution.resolve(), args.out.resolve(), resume=args.resume)))
