"""Filesystem-isolated continuation of an unchanged archived native worker.

The archive is mounted read-only at its original checkout path. Only a fresh
attempt directory and its fresh raw-output parent are writable. No stored plan,
protocol, source hash or input path is rewritten. The shared GPU gate and host
resource controls remain outside this filesystem adaptation.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.archived_sources import ArchivedSources


def mount_command(spec):
    prior = Path(spec['prior_root'])
    campaign = Path(spec['original_campaign'])
    raw_parent = Path(spec['original_raw_parent'])
    for field in ('prior_root', 'archive', 'original_campaign', 'original_raw_parent', 'attempt', 'raw_parent'):
        value = Path(spec[field])
        if not value.is_absolute() or '..' in value.parts or value == Path('/'):
            raise ValueError('Absolute non-root paths without traversal required')
    if (not campaign.is_relative_to(prior/'results/polymer_models/campaigns')
            or not raw_parent.is_relative_to(prior/'calcs/polymer_models')
            or not Path(spec['archive']).is_absolute()):
        raise ValueError('Only exact archived polymer campaign/output paths supported')
    attempt = Path(spec['attempt'])
    raw = Path(spec['raw_parent'])
    if (attempt.is_relative_to(prior) or raw.is_relative_to(prior)
            or attempt.is_relative_to(Path(spec['archive'])) or raw.is_relative_to(Path(spec['archive']))):
        raise ValueError('Fresh writable outputs must be outside the archived and original checkout roots')
    return ['bwrap', '--die-with-parent', '--ro-bind', '/', '/',
            '--ro-bind', spec['archive'], spec['prior_root'],
            '--bind', str(attempt), str(campaign),
            '--bind', str(raw), str(raw_parent),
            '--proc', '/proc', '--dev-bind', '/dev', '/dev', '--tmpfs', '/tmp',
            '--chdir', str(prior), '--']


def prepare(archive, prior, campaign, destination):
    reader = ArchivedSources(prior, archive, relative_to_prior=True)
    source = reader.path(campaign)
    plan_path = source/'plan.json'
    plan = json.loads(plan_path.read_text())
    task = plan['tasks'][0]
    if (len(plan['tasks']) != 1 or (task['charge'], task['mult']) != (0, 1)
            or plan['experimental_model']['solvent'] != 'acetonitrile'
            or plan['execution_driver'] != str(prior/'scripts/pipeline/smooth_polymer_native.py')
            or plan['maximum_host_RSS_GiB'] != 80 or plan['maximum_peak_GPU_MiB'] != 31744):
        raise ValueError('Only the unchanged, already declared neutral-singlet native class is supported')
    owner_path = source/'residency_guard/owner.json'
    if not owner_path.exists() or not json.loads(owner_path.read_text()).get('guard_stopped'):
        raise ValueError('Recorded source residency guard must have closed; inspect liveness separately')
    stage_path = source/'experimental_native/status.json'
    stage = json.loads(stage_path.read_text())
    if stage['status'] != 'experimental_native_failed_preserved' or 'finished' not in stage:
        raise ValueError('Preserve completed states and live/unknown attempts')
    reader.verify_all(plan['producer_sha256'])
    reader.verify_all(plan['input_sha256'])
    if sha256_file(source/'input.xyz') != plan['scored_geometry_sha256']:
        raise ValueError('Declared geometry differs')
    if destination.exists():
        raise FileExistsError('Use a fresh continuation destination')
    attempt, raw = destination/'campaign', destination/'raw'
    attempt.mkdir(parents=True)
    raw.mkdir()
    for name in ('plan.json', 'input.xyz'):
        shutil.copyfile(source/name, attempt/name)
    spec = dict(status='prepared_not_submitted', archive=str(archive), prior_root=str(prior),
                original_campaign=str(campaign), original_raw_parent=str(Path(plan['raw_root']).parent),
                original_raw_root=plan['raw_root'], attempt=str(attempt), raw_parent=str(raw),
                source_owner=stage['owner'], source_status_sha256=sha256_file(stage_path),
                source_guard_sha256=sha256_file(owner_path), plan_sha256=sha256_file(plan_path),
                geometry_sha256=plan['scored_geometry_sha256'], method_id=plan['experimental_method_id'],
                key=task['key'], maximum_host_RSS_GiB=80, minimum_free_GiB=32,
                historical_files_modified=False, scientific_settings_changed=False,
                complete_native_G_supplied=False)
    mount_command(spec)
    reader.recheck()
    atomic_write_json(destination/'namespace.json', spec)
    return spec


def verify_spec(path):
    spec = json.loads(path.read_text())
    attempt = Path(spec['attempt'])
    if sha256_file(attempt/'plan.json') != spec['plan_sha256']:
        raise ValueError('Preserved plan bytes changed')
    if sha256_file(attempt/'input.xyz') != spec['geometry_sha256']:
        raise ValueError('Preserved input geometry changed')
    return spec


def probe(path):
    spec = verify_spec(path)
    original = spec['original_campaign']
    # All source and input hashes are read through the mounted original paths.
    code = ('import hashlib,json,pathlib; c=pathlib.Path('+repr(original)+');'
            ' p=json.loads((c/"plan.json").read_text());'
            ' m={**p["producer_sha256"],**p["input_sha256"]};'
            ' bad=[k for k,v in m.items() if hashlib.sha256(pathlib.Path(k).read_bytes()).hexdigest()!=v];'
            ' assert not bad,bad;'
            ' print(json.dumps({"status":"original_path_namespace_hashes_match","bindings":len(m)}))')
    command = mount_command(spec)+['/usr/bin/python3', '-I', '-S', '-c', code]
    result = subprocess.run(command, capture_output=True, text=True, check=True, timeout=180)
    report = json.loads(result.stdout)
    profile = json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text())
    dependency_code = (
        'import importlib.metadata as m,json; '
        'import ase,rdkit,numpy,psutil; '
        'wanted='+repr(profile['versions'])+'; '
        'actual={name:m.version(name) for name in wanted}; '
        'assert wanted==actual,(wanted,actual); '
        'print(json.dumps({"versions":actual,"CPU_dependencies_imported":["ase","rdkit","numpy","psutil"]}))')
    dependency = subprocess.run(mount_command(spec)+[profile['python'], '-c', dependency_code],
                                capture_output=True, text=True, check=True, timeout=60)
    report['CPU_dependency_preflight'] = json.loads(dependency.stdout)
    report.update(namespace_sha256=sha256_file(path), quantum_jobs_submitted=0,
                  host=os.uname().nodename.split('.')[0])
    binary = Path(shutil.which('bwrap') or '')
    if not binary.is_file():
        raise ValueError('Namespace binary identity unavailable')
    report['namespace_binary'] = dict(path=str(binary.resolve()), sha256=sha256_file(binary))
    atomic_write_json(path.with_name('namespace_probe.json'), report)
    return report


def execute(path):
    spec = verify_spec(path)
    proof_path = path.with_name('namespace_probe.json')
    proof = json.loads(proof_path.read_text())
    if proof['namespace_sha256'] != sha256_file(path) or proof['status'] != 'original_path_namespace_hashes_match':
        raise ValueError('Original-path hash probe required before GPU execution')
    from redox.workflows.resources import verify_gate_owner
    verify_gate_owner()
    prior, campaign = Path(spec['prior_root']), Path(spec['original_campaign'])
    profile_path = Path(spec['archive'])/'environments/cutensor22/runtime.json'
    profile = json.loads(profile_path.read_text())
    environment = dict(os.environ, PYTHONPATH=str(prior/'src'),
                       LD_LIBRARY_PATH=profile['library_path'], PYTHONDONTWRITEBYTECODE='1')
    argv = mount_command(spec)+['env', 'PYTHONPATH='+str(prior/'src'),
        profile['python'], str(prior/'scripts/fleet/exclusive_python_worker.py'),
        '--guard-out', str(campaign/'residency_guard'),
        str(prior/'scripts/fleet/bounded_host_python_worker.py'), '--out', str(campaign/'host_guard'),
        '--maximum-rss-gib', '80', '--minimum-free-gib', '32',
        str(prior/'scripts/fleet/tensor_runtime_worker.py'), '--profile', str(prior/'environments/cutensor22/runtime.json'),
        '--out', str(campaign/'tensor_runtime'), str(prior/'scripts/pipeline/smooth_polymer_native.py'),
        'run', '--campaign', str(campaign)]
    os.execvpe(argv[0], argv, environment)


def launch_reserved(path, gpu, out):
    """Put the unchanged gate inside the mount namespace, before the worker.

    This preserves the strict immediate-parent lease check. Only the shared
    reservation directory is additionally writable; archived science stays RO.
    """
    spec = verify_spec(path)
    proof = json.loads(path.with_name('namespace_probe.json').read_text())
    if (proof['namespace_sha256'] != sha256_file(path)
            or proof.get('host') != os.uname().nodename.split('.')[0]):
        raise ValueError('Namespace requires its exact source/dependency probe')
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity, free_slots, gpu_suitability, GATE
    prior, case = Path(spec['prior_root']), Path(spec['original_campaign'])
    profile = json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text())
    if out.exists() or (Path(spec['attempt'])/'residency_guard').exists():
        raise FileExistsError('Preserve prior launches and scientific attempts')
    out.mkdir(parents=True)
    lock_root = GATE.parent.parent/'.gpu_locks'
    mount = mount_command(spec)
    mount[-1:-1] = ['--bind', str(lock_root), str(lock_root)]
    worker = [profile['python'], str(prior/'scripts/fleet/exclusive_python_worker.py'),
              '--guard-out', str(case/'residency_guard'),
              str(prior/'scripts/fleet/bounded_host_python_worker.py'), '--out', str(case/'host_guard'),
              '--maximum-rss-gib', '80', '--minimum-free-gib', '32',
              str(prior/'scripts/fleet/tensor_runtime_worker.py'),
              '--profile', str(prior/'environments/cutensor22/runtime.json'),
              '--out', str(case/'tensor_runtime'), str(prior/'scripts/pipeline/smooth_polymer_native.py'),
              'run', '--campaign', str(case)]
    command = mount+[str(GATE), 'run', str(gpu), '--label', 'redox-archived-native', '--', *worker]
    environment = dict(os.environ, PYTHONPATH=str(prior/'src'), LD_LIBRARY_PATH=profile['library_path'],
                       PYTHONDONTWRITEBYTECODE='1')
    for key in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'XTB_THREADS'):
        environment[key] = '2'
    with allocation_guard():
        capacity = host_capacity(2, minimum_ram_gib=112)
        promised = remaining_host_promises()
        if capacity['available_RAM_bytes']-promised < 112*2**30:
            raise RuntimeError('RAM already promised to live bounded workers; no attempt started')
        cpus = cpu_allocation(2)
        suitability = gpu_suitability(gpu)
        host = os.uname().nodename.split('.')[0]
        if (host, gpu) not in free_slots([host]):
            raise RuntimeError('Device no longer available; no scientific attempt started')
        launcher_bytes = Path(__file__).read_bytes()
        (out/'launcher_source.py').write_bytes(launcher_bytes)
        with (out/'command.log').open('x') as log:
            child = subprocess.Popen(['taskset', '--cpu-list', ','.join(map(str, cpus)), *command],
                                     stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
                                     env=environment, start_new_session=True,
                                     preexec_fn=lambda: os.sched_setaffinity(0, cpus))
        record = dict(status='namespace_gate_started_scientific_start_unverified', host=host,
                      pid=child.pid, process_start=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')', 1)[1].split()[19],
                      GPU=gpu, CPU_affinity=cpus, host_capacity=capacity, gpu_preflight=suitability,
                      namespace=str(path), namespace_sha256=sha256_file(path), argv=command,
                      producer_sha256=sha256_file(out/'launcher_source.py'), started=time.time(),
                      quantum_result_qualified=False, gate_inside_namespace=True)
        record['remaining_other_host_promises_bytes'] = promised
        atomic_write_json(out/'launch.json', record)
    return record


def remaining_host_promises():
    """Reserve the unused portion of every live bounded Python worker's RAM."""
    import psutil
    total = 0
    for process in psutil.process_iter(['pid', 'cmdline', 'memory_info']):
        try:
            args = process.info['cmdline'] or []
            if (not args or 'python' not in Path(args[0]).name
                    or '--maximum-rss-gib' not in args):
                continue
            bound = float(args[args.index('--maximum-rss-gib')+1])*2**30
            rss = process.info['memory_info'].rss
            if not 0 < bound < float('inf'):
                raise ValueError('Live worker has an invalid declared RAM bound')
            total += max(0, bound-rss)
        except psutil.NoSuchProcess:
            continue
    return int(total)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('prepare', 'probe', 'execute', 'launch'))
    p.add_argument('--spec', type=Path)
    p.add_argument('--gpu', type=int)
    p.add_argument('--out', type=Path)
    for name in ('archive', 'prior-root', 'campaign', 'destination'):
        p.add_argument('--'+name, type=Path)
    a = p.parse_args()
    if a.action == 'prepare':
        if not all((a.archive, a.prior_root, a.campaign, a.destination)):
            p.error('Explicit archive, prior root, campaign and fresh destination required')
        result = prepare(a.archive.resolve(), a.prior_root.resolve(), a.campaign,
                         a.destination.resolve())
    elif a.action == 'probe':
        result = probe(a.spec.resolve())
    elif a.action == 'launch':
        if a.gpu is None or a.out is None:
            p.error('Explicit GPU and fresh launch directory required')
        result = launch_reserved(a.spec.resolve(), a.gpu, a.out.resolve())
    else:
        execute(a.spec.resolve())
        return
    print(json.dumps(result))
    if a.action == 'launch':
        # bwrap's death binding requires this controller to remain alive.
        # Submit this CLI detached, retaining its stdout/stderr on shared storage.
        _, status = os.waitpid(result['pid'], 0)
        result.update(status='namespace_process_returned_physical_review_required',
                      returncode=os.waitstatus_to_exitcode(status), finished=time.time())
        atomic_write_json(a.out.resolve()/'completion.json', result)


if __name__ == '__main__':
    main()
