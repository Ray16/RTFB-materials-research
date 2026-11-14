#!/usr/bin/env python3
"""Package native inputs or test an extracted cache in a fresh NumPy environment."""
import argparse
import json
import os
from pathlib import Path, PurePosixPath
import shutil
import subprocess
import sys
import tarfile
import time

from native_cache import local, sha


BOOTSTRAP = '''import os,sys,runpy,json,importlib.util
forbidden=json.loads(sys.argv[1])
def guard(event,args):
 if event=="open" and isinstance(args[0],(str,bytes)):
  path=os.path.realpath(os.fsdecode(args[0]))
  if any(path==p or path.startswith(p+os.sep) for p in forbidden):
   raise RuntimeError("Original scientific path blocked: "+path)
sys.addaudithook(guard)
for name in ("redox","ase","pyscf","torch"):
 assert importlib.util.find_spec(name) is None, name
sys.argv=sys.argv[2:]
runpy.run_path(sys.argv[0],run_name="__main__")
'''


def inventory(root):
    manifest = json.loads((root/'manifest.json').read_text())
    actual = {str(p.relative_to(root)) for p in root.rglob('*') if p.is_file()}
    if actual != set(manifest['files_sha256'])|{'manifest.json'}:
        raise ValueError('Raw package inventory differs')
    for relative, expected in manifest['files_sha256'].items():
        path = local(root, relative)
        if path.is_symlink() or sha(path) != expected:
            raise ValueError('Raw package content differs: '+relative)
    return manifest


def compress(package, archive):
    if archive.exists():
        raise FileExistsError('Preserve earlier compressed packages')
    manifest_hash = sha(package/'manifest.json')
    inventory(package)
    command = ['tar', '--sort=name', '--mtime=@0', '--owner=0', '--group=0', '--numeric-owner',
               '-cf', '-', '-C', str(package.parent), package.name]
    producer = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        result = subprocess.run(['zstd', '-T1', '-3', '-q', '-o', str(archive)],
                                stdin=producer.stdout, capture_output=True, text=True)
        producer.stdout.close()
        _, error = producer.communicate()
        if result.returncode or producer.returncode:
            raise RuntimeError('Package compression failed: '+error.decode()+result.stderr)
    finally:
        if producer.poll() is None:
            producer.terminate()
            producer.wait()
    inventory(package)
    if sha(package/'manifest.json') != manifest_hash:
        raise ValueError('Raw cache manifest changed during compression')
    return dict(status='native_cache_compressed', archive=str(archive), archive_sha256=sha(archive),
                archive_bytes=archive.stat().st_size, manifest_sha256=manifest_hash,
                compression_threads=1, deterministic_tar_metadata=True)


def safe_extract(stream, destination, expected_root):
    if not expected_root or PurePosixPath(expected_root).name != expected_root:
        raise ValueError('One expected archive root name is required')
    files, total = set(), 0
    for item in stream:
        path = PurePosixPath(item.name)
        if (not path.parts or path.is_absolute() or '..' in path.parts or path.parts[0] != expected_root
                or not (item.isfile() or item.isdir())):
            raise ValueError('Unsafe or unexpected native-cache archive member')
        target = local(destination, str(path))
        if item.isdir():
            target.mkdir(parents=True, exist_ok=True)
            continue
        total += item.size
        if item.name in files or item.size > 2**30 or total > 2**31 or len(files) >= 10000:
            raise ValueError('Duplicate or oversized native-cache archive')
        target.parent.mkdir(parents=True, exist_ok=True)
        with stream.extractfile(item) as incoming, target.open('xb') as outgoing:
            shutil.copyfileobj(incoming, outgoing)
        files.add(item.name)
    return len(files)


def validate(archive, name, output, python, forbidden):
    if output.exists():
        raise FileExistsError('Preserve previous isolated validations')
    output.mkdir(parents=True)
    archive_hash = sha(archive)
    payload = output/'payload'
    payload.mkdir()
    decoder = subprocess.Popen(['zstd', '-dc', str(archive)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    try:
        with tarfile.open(fileobj=decoder.stdout, mode='r|') as stream:
            count = safe_extract(stream, payload, name)
        decoder.stdout.close()
        _, error = decoder.communicate()
        if decoder.returncode:
            raise ValueError('Native cache archive cannot be decoded: '+error.decode())
    finally:
        if decoder.poll() is None:
            decoder.terminate()
            decoder.wait()
    package = payload/name
    manifest = inventory(package)
    environment = output/'venv'
    env = dict(PATH='/usr/bin:/bin', LC_ALL='C.UTF-8', **{key:'1' for key in (
        'OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'XTB_THREADS')})
    with (output/'environment.log').open('x') as log:
        subprocess.run([str(python), '-m', 'venv', '--copies', str(environment)], env=env,
                       stdout=log, stderr=subprocess.STDOUT, check=True, timeout=180)
        executable = environment/'bin/python'
        subprocess.run([str(executable), '-m', 'pip', 'install', '--require-hashes', '--no-deps',
                        '-r', str(package/'requirements.lock.txt')], env=env,
                       stdout=log, stderr=subprocess.STDOUT, check=True, timeout=300)
        subprocess.run([str(executable), '-m', 'pip', 'check'], env=env,
                       stdout=log, stderr=subprocess.STDOUT, check=True, timeout=30)
    blocked = [str(p.resolve()) for p in forbidden]
    (output/'isolation_bootstrap.py').write_text(BOOTSTRAP)
    command = [str(executable), '-B', '-I', '-c', BOOTSTRAP, json.dumps(blocked),
               str(package/'replay.py'), 'replay', '--package', str(package), '--out', str(output/'replay.json')]
    with (output/'replay.log').open('x') as log:
        subprocess.run(command, cwd='/', env=env, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=600)
    replay = json.loads((output/'replay.json').read_text())
    if replay['states'] != len(manifest['states']) or replay['quantum_calculations_run']:
        raise ValueError('Extracted replay coverage or scope differs')
    if sha(archive) != archive_hash:
        raise ValueError('Archive changed during isolated validation')
    result = dict(status='extracted_native_cache_fresh_numpy_replay_checked', checked=time.time(),
                  archive=str(archive), archive_sha256=archive_hash, archive_files=count,
                  package_manifest_sha256=sha(package/'manifest.json'), states=replay['states'],
                  maximum_frequency_error_cm=max(r['frequency_error_cm'] for r in replay['records']),
                  maximum_nuclear_error_eV=max(max(r['nuclear_errors_eV'].values()) for r in replay['records']),
                  blocked_original_paths=blocked, original_project_installed=False,
                  scientific_requirements_sha256=sha(package/'requirements.lock.txt'),
                  scientific_dependency='numpy==2.2.6', quantum_calculations_run=False,
                  full_historical_source_closure_revalidated=False, public_archive_deposited=False)
    (output/'validation.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('compress', 'validate'))
    for name in ('package', 'archive', 'out'):
        parser.add_argument('--'+name, type=Path)
    parser.add_argument('--name')
    parser.add_argument('--python', type=Path, default=Path(sys.executable))
    parser.add_argument('--forbid', type=Path, action='append', default=[])
    args = parser.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard():
        slots = 2 if args.action == 'compress' else 1
        host_capacity(slots, minimum_ram_gib=4)
        os.sched_setaffinity(0, cpu_allocation(slots))
    if args.action == 'compress':
        result = compress(args.package.resolve(), args.archive.resolve())
        with args.out.open('x') as stream:
            json.dump(result, stream, indent=2)
            stream.write('\n')
    else:
        result = validate(args.archive.resolve(), args.name, args.out.resolve(), args.python.resolve(), args.forbid)
    print(json.dumps(result))
