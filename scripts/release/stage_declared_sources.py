"""Stage a fresh local raw-evidence package from the resolved source inventory."""
import argparse
import json
import os
from pathlib import Path
import shutil

from redox.core.hashing import sha256_file
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def object_layout(mappings, all_relocations):
    """Shared inclusion and content identity rules for staged/streamed packages."""
    excluded = {'runtime_dependency', 'third_party_publication'}
    external = [r for r in mappings if r['category'] in excluded]
    included = [r for r in mappings if r['category'] not in excluded]
    objects = {}
    for row in mappings:
        digest=row['sha256']; size=row['bytes']
        if (len(digest)!=64 or any(c not in '0123456789abcdef' for c in digest)
                or type(size) is not int or size<0):
            raise ValueError('Exact content identity and nonnegative size required')
    for row in included:
        if row['sha256'] in objects and objects[row['sha256']] != row['bytes']:
            raise ValueError('One content hash cannot have conflicting sizes')
        objects[row['sha256']] = row['bytes']
    external_hashes = {r['sha256'] for r in external}
    if any(row['sha256'] not in objects and row['sha256'] not in external_hashes for row in all_relocations):
        raise ValueError('Relocation has no inventoried dependency')
    return dict(objects=objects, files=included, external_dependencies=external,
                relocations=[r for r in all_relocations if r['sha256'] in objects],
                external_relocations=[r for r in all_relocations if r['sha256'] not in objects])


def stage(inventory, content_resolution, metadata_resolution, output, *, resume=False):
    if output.exists() and not resume:
        raise FileExistsError('Use a fresh release directory')
    if resume and (not output.is_dir() or {p.name for p in output.iterdir()} != {'objects'}):
        raise ValueError('Resume only an unfinished objects-only staging directory')
    inv = json.loads(inventory.read_text())
    if content_resolution is None and metadata_resolution is None:
        if inv.get('status') != 'declared_dependency_inventory_verified' or inv.get('unresolved') != []:
            raise ValueError('Completely resolved declared inventory required')
        content = dict(resolved=[]); metadata = dict(resolved=[])
    else:
        if content_resolution is None or metadata_resolution is None:
            raise ValueError('Supply both resolution records or neither')
        content = json.loads(content_resolution.read_text())
        metadata = json.loads(metadata_resolution.read_text())
        if (sha256_file(inventory) != content['inventory_sha256']
                or sha256_file(content_resolution) != metadata['previous_resolution_sha256']
                or metadata['unresolved']):
            raise ValueError('Exact resolved inventory chain required')
    mappings = list(inv['files']) + [dict(path=r['path'], sha256=r['sha256'],
        bytes=r['bytes'], category='unversioned_research_evidence_or_source_snapshot')
        for r in metadata['resolved']]
    layout=object_layout(mappings,[*inv['explicit_archive_relocations'], *content['resolved'], *metadata['resolved']])
    included=layout['files'];objects=layout['objects'];external=layout['external_dependencies']
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
    shutil.copyfile(Path(__file__).with_name('verify_source_objects.py'), output / 'verify.py')
    shutil.copyfile(Path(__file__), output / 'stage_source.py')
    (output / 'README.md').write_text('''# Local raw-evidence package

Run `python -I verify.py` from this directory, or pass the package directory
explicitly. The standard-library verifier checks every content-addressed raw
object and its file/relocation mapping without original scientific paths.
The manifest identifies the selected source inventory and its traversal scope.
Aggregate provenance can include broader historical calculations; inclusion
does not promote a failed or unfinished record into an accepted result.

This is a local review package, not a public deposit or a numerical replay.
SCF, Hessian, spin and reaction calculations are not rerun by the verifier.
The original records remain byte-identical; path mappings are external to them.
Runtime dependencies and third-party publications are listed in the manifest
but not redistributed or verified here. Technical inventory categories confer
no license. Existing code rights remain applicable; data reuse terms and public
deposition remain pending. No experimental or ensemble qualification is added.
''')
    data = dict(status='local_declared_raw_evidence_staged_not_deposited',
        **layout,
        inventory_scope=inv['traversal_scope'], implicit_dependencies_covered=False,
        public_deposition_complete=False, data_license=None,
        support_files_sha256={name:sha256_file(output/name) for name in ('verify.py','stage_source.py','README.md')},
        source_sha256={str(p): sha256_file(p) for p in
                       (inventory, content_resolution, metadata_resolution, Path(__file__)) if p is not None})
    (output / 'manifest.json').write_text(json.dumps(data, indent=2) + '\n')
    return dict(objects=len(objects), bytes=sum(objects.values()), external_dependencies=len(external))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ('inventory', 'out'):
        parser.add_argument('--' + key, type=Path, required=True)
    for key in ('content-resolution', 'metadata-resolution'):
        parser.add_argument('--' + key, type=Path)
    parser.add_argument('--resume', action='store_true', help='Rehash an unfinished objects-only directory')
    args = parser.parse_args()
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4); os.sched_setaffinity(0, cpu_allocation(1))
    print(json.dumps(stage(args.inventory.resolve(),
                           args.content_resolution.resolve() if args.content_resolution else None,
                           args.metadata_resolution.resolve() if args.metadata_resolution else None,
                           args.out.resolve(), resume=args.resume)))
