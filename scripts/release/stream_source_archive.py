#!/usr/bin/env python3
"""Stream resolved evidence into a checked archive without an uncompressed copy."""
import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import tarfile

from redox.core.hashing import sha256_file
from stage_declared_sources import object_layout


PREFIX='publication_raw_evidence'


class CheckedReader:
    def __init__(self, stream, destination):
        self.stream=stream;self.destination=destination;self.digest=hashlib.sha256();self.count=0

    def read(self, size):
        block=self.stream.read(size);self.digest.update(block);self.count+=len(block)
        if self.count//(64*2**20)!=(self.count-len(block))//(64*2**20):
            if shutil.disk_usage(self.destination).free<32*2**30:
                raise RuntimeError('Shared storage fell below the 32-GiB reserve')
        return block


def create(inventory_path, archive, metadata):
    if archive.exists() or metadata.exists():raise FileExistsError('Preserve earlier release outputs')
    inv=json.loads(inventory_path.read_text())
    if inv.get('status')!='declared_dependency_inventory_verified' or inv.get('unresolved')!=[]:
        raise ValueError('Complete declared-map inventory required')
    layout=object_layout(inv['files'],inv['explicit_archive_relocations'])
    # Use the uncompressed size as a conservative upper bound; no guessed ratio.
    if shutil.disk_usage(archive.parent).free<sum(layout['objects'].values())+32*2**30:
        raise RuntimeError('Insufficient worst-case archive headroom and shared reserve')
    metadata.mkdir()
    here=Path(__file__).parent
    for source,name in [(here/'verify_source_objects.py','verify.py'),
                        (Path(__file__),'archive_source.py'),
                        (here/'stage_declared_sources.py','stage_declared_sources.py')]:
        shutil.copyfile(source,metadata/name)
    (metadata/'README.md').write_text('''# Publication raw evidence — local release candidate

After extraction, run `python -I verify.py` from this directory. This verifies
all included object bytes and logical-path mappings without original paths.
The manifest describes the frozen publication inputs and their declared
dependency maps. Aggregate provenance includes historical records; failed and
unfinished records retain their original status. Inclusion is not physical QC.

Runtime dependencies and third-party publications remain external and are
listed separately. Numerical replays, experimental validation, implicit runtime
dependencies, reuse licensing and public deposition retain separate scopes.
The archived producer sources preserve packaging logic; code/data rights are
not reassigned by this local bundle. No fresh quantum calculation is performed.
''')
    support={p.name:sha256_file(p) for p in metadata.iterdir()}
    manifest=dict(status='streamed_publication_raw_evidence_not_deposited',**layout,
        inventory_scope=inv['traversal_scope'],implicit_dependencies_covered=False,
        support_files_sha256=support,inventory_sha256=sha256_file(inventory_path),
        inventory_source=str(inventory_path),public_deposition_complete=False,data_license=None)
    manifest_path=metadata/'manifest.json'
    manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
    sources={}
    for row in layout['files']:sources.setdefault(row['sha256'],Path(row['path']))
    compressor=subprocess.Popen(['zstd','-T1','-3','-q','-o',str(archive)],stdin=subprocess.PIPE)

    def add(stream,name,path,size,expected):
        before=path.stat()
        if before.st_size!=size:raise ValueError('Source size changed: '+str(path))
        info=tarfile.TarInfo(PREFIX+'/'+name);info.size=size;info.mode=0o644;info.mtime=0
        with path.open('rb') as source:
            reader=CheckedReader(source,archive.parent);stream.addfile(info,reader)
        after=path.stat()
        if (reader.digest.hexdigest()!=expected or
                (before.st_size,before.st_mtime_ns,before.st_ctime_ns)!=
                (after.st_size,after.st_mtime_ns,after.st_ctime_ns)):
            raise ValueError('Source bytes changed during archiving: '+str(path))

    try:
        with tarfile.open(fileobj=compressor.stdin,mode='w|') as stream:
            for name in sorted([*support,'manifest.json']):
                path=metadata/name;add(stream,name,path,path.stat().st_size,sha256_file(path))
            for index,(digest,size) in enumerate(sorted(layout['objects'].items()),1):
                if shutil.disk_usage(archive.parent).free<32*2**30:
                    raise RuntimeError('Shared storage reserve reached')
                add(stream,'objects/'+digest,sources[digest],size,digest)
                if index%1000==0:print(json.dumps(dict(objects_written=index)),flush=True)
        compressor.stdin.close()
        if compressor.wait()!=0:raise RuntimeError('Compression failed; partial output preserved')
    except BaseException:
        try:compressor.stdin.close()
        except BrokenPipeError:pass
        compressor.terminate();compressor.wait();raise
    return dict(status='source_checked_archive_created_verification_pending',
        archive=str(archive),objects=len(layout['objects']),bytes=sum(layout['objects'].values()),
        compressed_bytes=archive.stat().st_size,manifest_sha256=sha256_file(manifest_path),
        uncompressed_copy_created=False,public_deposition_complete=False)


def verify_archive(archive, manifest_path):
    manifest=json.loads(manifest_path.read_text())
    expected={PREFIX+'/objects/'+h:h for h in manifest['objects']}
    expected.update({PREFIX+'/'+name:digest for name,digest in manifest['support_files_sha256'].items()})
    expected[PREFIX+'/manifest.json']=sha256_file(manifest_path)
    found=set();decoder=subprocess.Popen(['zstd','-dc',str(archive)],stdout=subprocess.PIPE)
    try:
        with tarfile.open(fileobj=decoder.stdout,mode='r|') as stream:
            for member in stream:
                name=PurePosixPath(member.name)
                if (not member.isfile() or name.is_absolute() or '..' in name.parts
                        or member.name not in expected or member.name in found):
                    raise ValueError('Unexpected, unsafe or duplicate archive member')
                digest=hashlib.sha256()
                with stream.extractfile(member) as source:
                    for block in iter(lambda:source.read(2**20),b''):digest.update(block)
                if digest.hexdigest()!=expected[member.name]:raise ValueError('Archive member bytes differ')
                found.add(member.name)
        decoder.stdout.close()
        if decoder.wait()!=0 or found!=set(expected):raise ValueError('Incomplete archive')
    except BaseException:
        decoder.stdout.close();decoder.terminate();decoder.wait();raise
    return dict(status='all_streamed_archive_members_verified',archive=str(archive),
        archive_sha256=sha256_file(archive),archive_bytes=archive.stat().st_size,
        members_verified=len(found),objects=len(manifest['objects']),
        manifest_sha256=sha256_file(manifest_path),uncompressed_copy_created=False,
        numerical_replay_performed=False,public_deposition_complete=False)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['create','verify'])
    parser.add_argument('--inventory',type=Path)
    parser.add_argument('--archive',type=Path,required=True)
    parser.add_argument('--metadata',type=Path,required=True)
    parser.add_argument('--out',type=Path)
    args=parser.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():
        host_capacity(2,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(2))
    if args.action=='create':
        if args.inventory is None:parser.error('--inventory required for create')
        result=create(args.inventory.resolve(),args.archive.resolve(),args.metadata.resolve())
    else:
        if args.out is None or args.out.exists():parser.error('--out must name a fresh verification record')
        result=verify_archive(args.archive.resolve(),args.metadata.resolve()/'manifest.json')
        args.out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result))
