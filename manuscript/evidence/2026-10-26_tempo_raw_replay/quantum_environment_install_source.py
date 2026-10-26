#!/usr/bin/env python
"""Restore recorded third-party packages and an exact local project snapshot.

This CPU-only operation creates a new prefix. It neither imports CUDA packages
nor modifies the qualified environment. Quantum qualification is a separate step.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import socket
import subprocess
import sys
import tarfile
import time


CAPS = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS',
        'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'XTB_THREADS')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def split_lock(text):
    """Never obtain the local project from an unrelated package-index namesake."""
    packages = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        match = re.fullmatch(r'([A-Za-z0-9_.-]+)==([^\s;]+)', line)
        if not match:
            raise ValueError('Every recorded distribution needs one exact version')
        name, version = match.groups()
        name = re.sub(r'[-_.]+', '-', name).lower()
        if name in packages:
            raise ValueError('Duplicate normalized distribution')
        packages[name] = version
    if packages.pop('redox', None) != '0.1.0':
        raise ValueError('Expected the exact local project distribution')
    torch = packages.pop('torch', None)
    if torch != '2.8.0+cu128':
        raise ValueError('Recorded CUDA-specific torch build changed')
    return [f'{name}=={version}' for name, version in sorted(packages.items())], 'torch==' + torch


def atomic_json(path, data):
    temporary = path.with_name(path.name + '.tmp')
    temporary.write_text(json.dumps(data, indent=2) + '\n')
    temporary.replace(path)


def restore(source, out, conda, pip_lock=None, cache_root=None):
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    if any(os.environ.get(key) != '1' for key in CAPS):
        raise ValueError('CPU installer requires all native caps set to one')
    with allocation_guard():
        resources = host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    source, out, conda = source.resolve(), out.resolve(), conda.resolve()
    if out.exists():
        raise FileExistsError('Preserve existing restoration attempts; use a fresh output directory')
    out.mkdir(parents=True)
    prefix, checkout = out / 'environment', out / 'source'
    owner = dict(host=socket.gethostname(), pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19],
        CPU_affinity=sorted(os.sched_getaffinity(0)), thread_caps={k: os.environ[k] for k in CAPS})
    state = dict(status='preparing_fresh_installation', owner=owner, started=time.time(),
        prefix=str(prefix), source_checkout=str(checkout), resources=resources,
        original_environment_modified=False, GPU_calculation_performed=False,
        restoration_producer_sha256=digest(__file__))
    (out / 'restoration_producer.py').write_bytes(Path(__file__).read_bytes())
    atomic_json(out / 'status.json', state)
    cache = cache_root.resolve() if cache_root is not None else out
    environment = dict(os.environ, PYTHONNOUSERSITE='1', CONDA_PKGS_DIRS=str(cache / 'conda_packages'),
        PIP_CACHE_DIR=str(cache / 'pip_cache'), PIP_DISABLE_PIP_VERSION_CHECK='1',
        MAX_JOBS='1', CMAKE_BUILD_PARALLEL_LEVEL='1', MAKEFLAGS='-j1')
    environment.pop('PYTHONPATH', None)
    environment.pop('PYTHONHOME', None)
    commands = []

    def run(label, argv):
        state.update(status=label, updated=time.time())
        commands.append(dict(stage=label, argv=[str(x) for x in argv]))
        atomic_json(out / 'commands.json', commands)
        atomic_json(out / 'status.json', state)
        with (out / (label + '.log')).open('x') as log:
            subprocess.run([str(x) for x in argv], check=True, cwd=out,
                env=environment, stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT)

    try:
        commit = subprocess.check_output(['git', '-C', str(source), 'rev-parse', 'HEAD'], text=True).strip()
        state['source_commit'] = commit
        archive = out / 'source.tar'
        subprocess.run(['git', '-C', str(source), 'archive', '--format=tar', '--output=' + str(archive),
            commit, 'src', 'config', 'scripts', 'library', 'environments',
            'pyproject.toml', 'README.md', 'requirements.txt'], check=True)
        checkout.mkdir()
        with tarfile.open(archive) as contents:
            for member in contents.getmembers():
                if member.issym() or member.islnk() or not (checkout / member.name).resolve().is_relative_to(checkout):
                    raise ValueError('Unsafe project snapshot member')
            contents.extractall(checkout)
        closure_source = Path(__file__).with_name('certificate_closure.py')
        spec = importlib.util.spec_from_file_location('restoration_certificate_closure', closure_source)
        certificate = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(certificate)
        certificate.closure(source, checkout, out / 'certificate_closure.json')
        (out / 'certificate_closure_producer.py').write_bytes(closure_source.read_bytes())
        state['certificate_closure_sha256'] = digest(out / 'certificate_closure.json')
        runtime = json.loads((checkout / 'environments/quantum_runtime.json').read_text())
        for relative, expected in runtime['files_sha256'].items():
            if digest(checkout / relative) != expected:
                raise ValueError('Recorded environment inputs changed')
        recorded_lock = checkout / 'environments/quantum-pip.lock.txt'
        lock = recorded_lock
        if pip_lock is not None:
            lock = out / 'selected-pip.lock.txt'
            lock.write_bytes(pip_lock.resolve().read_bytes())
        third_party, torch = split_lock(lock.read_text())
        external = out / 'third_party.lock.txt'
        external.write_text('\n'.join(third_party) + '\n')
        state.update(source_archive_sha256=digest(archive), pip_lock_sha256=digest(lock),
            recorded_pip_lock_sha256=digest(recorded_lock), selected_pip_lock=str(lock),
            selected_lock_differs_from_recorded=digest(lock) != digest(recorded_lock),
            third_party_distributions=len(third_party) + 1, local_distribution='redox==0.1.0 from source_commit')
        run('conda_install', [conda, 'create', '--yes', '--prefix', prefix, '--file',
            checkout / 'environments/quantum-linux-64.explicit.txt'])
        python = prefix / 'bin/python'
        pip = [python, '-s', '-m', 'pip', 'install', '--no-deps', '--no-build-isolation']
        run('torch_install', [*pip, '--index-url', 'https://download.pytorch.org/whl/cu128',
            '--report', out / 'torch_install.json', torch])
        run('third_party_install', [*pip, '--index-url', 'https://pypi.org/simple',
            '--report', out / 'third_party_install.json', '-r', external])
        run('project_install', [*pip, '--report', out / 'project_install.json', '--editable', checkout])
        run('dependency_check', [python, '-s', '-m', 'pip', 'check'])
        library = prefix / 'lib/libxtb.so.6.7.1'
        if digest(library) != runtime['libxtb_sha256']:
            raise ValueError('Restored xTB library differs from the qualified bytes')
        run('pip_freeze', [python, '-s', '-m', 'pip', 'freeze', '--all'])
        artifacts = []
        for name in ('torch_install', 'third_party_install'):
            report = json.loads((out / (name + '.json')).read_text())
            for item in report['install']:
                artifacts.append(dict(name=item['metadata']['name'], version=item['metadata']['version'],
                    download=item['download_info']))
        atomic_json(out / 'installed_artifacts.json', artifacts)
        state.update(status='fresh_environment_installed_GPU_qualification_pending',
            libxtb_sha256=digest(library), artifact_records=len(artifacts),
            dependency_check_passed=True)
    except BaseException as error:
        state.update(status='fresh_environment_installation_failed_preserved', error=repr(error))
        raise
    finally:
        state['finished'] = time.time()
        atomic_json(out / 'status.json', state)
    return state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--conda', type=Path, required=True)
    parser.add_argument('--pip-lock', type=Path, help='Explicit prospective compatibility lock; never rewrites the recorded lock')
    parser.add_argument('--cache-root', type=Path, help='Reuse task-owned downloaded artifacts, without cloning an installed environment')
    args = parser.parse_args()
    print(json.dumps(restore(args.source, args.out, args.conda, args.pip_lock, args.cache_root)))


if __name__ == '__main__':
    main()
