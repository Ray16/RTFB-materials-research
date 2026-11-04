"""Create new native states from checked namespace energy receipts.

Only provenance paths are relocated into a derived parent receipt. Original
records remain byte-identical; numerical fields, state and method are unchanged.
The archived native preparer and worker supply every scientific calculation.
"""
import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.namespace_sources import NamespaceSources, read_checked_review
from redox.workflows.archived_native_namespace import remaining_host_promises, verify_spec


def portable_parent(report, spec, original_report):
    original_hash = sha256_file(original_report)
    if json.loads(original_report.read_text()) != report:
        raise ValueError('Only the exact independently checked report may be relocated')
    reader = NamespaceSources(spec)
    reader.verify_all(report['source_sha256'])
    sources = {str(reader.path(path)): digest for path, digest in report['source_sha256'].items()}
    sources[str(original_report)] = original_hash
    result = dict(report, campaign=spec['attempt'], source_sha256=sources)
    result['path_relocation'] = dict(
        original_review=str(original_report), original_review_sha256=original_hash,
        original_campaign=report['campaign'], relocated_campaign=spec['attempt'],
        permitted_changes=['campaign path', 'source-binding paths', 'this explicit relocation provenance'],
        numerical_or_chemical_fields_changed=False)
    if ({k: v for k, v in result.items() if k not in ('campaign', 'source_sha256', 'path_relocation')}
            != {k: v for k, v in report.items() if k not in ('campaign', 'source_sha256')}):
        raise ValueError('Portable parent changed scientific metadata')
    reader.recheck()
    if sha256_file(original_report) != original_hash:
        raise ValueError('Independent parent changed during relocation')
    return result


def namespace_command(spec):
    prior, archive = Path(spec['prior_root']), Path(spec['archive'])
    batch = Path(spec['original_campaign']).parent
    raw_batch = Path(spec['original_raw_parent'])
    for field in ('prior_root', 'archive', 'original_campaign', 'original_raw_parent', 'attempt', 'raw_parent'):
        value = Path(spec[field])
        if not value.is_absolute() or '..' in value.parts or value == Path('/'):
            raise ValueError('Absolute non-root namespace paths without traversal required')
    for field in ('attempt', 'raw_parent'):
        value = Path(spec[field]).resolve()
        if value.is_relative_to(prior.resolve()) or value.is_relative_to(archive.resolve()):
            raise ValueError('New native outputs must be outside the historical and current checkouts')
    if (not batch.is_relative_to(prior/'results/polymer_models/campaigns')
            or not raw_batch.is_relative_to(prior/'calcs/polymer_models')):
        raise ValueError('Explicit polymer campaign/raw batches required')
    command = ['bwrap', '--die-with-parent', '--ro-bind', '/', '/',
               '--ro-bind', str(archive), str(prior)]
    for relative in ('results/polymer_models/campaigns', 'calcs/polymer_models'):
        original_parent = prior/relative
        command += ['--tmpfs', str(original_parent)]
        for child in sorted((archive/relative).iterdir()):
            if child.is_dir():
                command += ['--ro-bind', str(child), str(original_parent/child.name)]
    command += ['--bind', str(Path(spec['attempt']).parent), str(batch),
                '--bind', spec['raw_parent'], str(raw_batch), '--proc', '/proc',
                '--dev-bind', '/dev', '/dev', '--tmpfs', '/tmp', '--chdir', str(prior), '--']
    return command


def prepare(energy_spec_path, destination):
    report, binding = read_checked_review(energy_spec_path, 'energy')
    energy_spec = json.loads(energy_spec_path.read_text())
    if (destination.resolve().is_relative_to(Path(energy_spec['prior_root']).resolve())
            or destination.resolve().is_relative_to(Path(energy_spec['archive']).resolve())):
        raise ValueError('Use a fresh durable native destination outside the source checkouts')
    if destination.exists():
        raise FileExistsError('Preserve existing derived-native preparation')
    destination.mkdir(parents=True)
    prior, archive = Path(energy_spec['prior_root']), Path(energy_spec['archive'])
    key = energy_spec['key']
    batch_name = ('2026-10-30_occupied_native_from_checked_energy' if report['task'].get('occupied_sites')
                  else '2026-10-30_empty_native_from_checked_energy')
    physical_parent, raw = destination/'campaign_parent', destination/'raw'
    physical_parent.mkdir()
    raw.mkdir()
    native_case = prior/'results/polymer_models/campaigns'/batch_name/key
    spec = dict(status='new_native_dependency_preparation_pending', archive=str(archive), prior_root=str(prior),
                key=key, method_id=energy_spec['method_id'], original_campaign=str(native_case),
                original_raw_parent=str(prior/'calcs/polymer_models'/batch_name),
                original_raw_root=str(prior/'calcs/polymer_models'/batch_name/key),
                attempt=str(physical_parent/key), raw_parent=str(raw),
                maximum_host_RSS_GiB=80, minimum_free_GiB=32,
                source_energy_namespace=str(energy_spec_path), source_energy_namespace_sha256=sha256_file(energy_spec_path),
                source_energy_review=binding, native_G_computed=False, primary_candidate_admission=False)
    parent_path = destination/'portable_electronic_parent.json'
    atomic_write_json(parent_path, portable_parent(report, energy_spec, Path(binding['review'])))
    profile = json.loads((archive/'environments/cutensor22/runtime.json').read_text())
    command = namespace_command(spec)+['env', 'PYTHONPATH='+str(prior/'src'), profile['python'],
        str(prior/'scripts/pipeline/smooth_polymer_native.py'), 'prepare', '--energy-case',
        energy_spec['attempt'], '--review', str(parent_path), '--campaign', str(native_case)]
    env = dict(os.environ, PYTHONPATH=str(prior/'src'), PYTHONDONTWRITEBYTECODE='1')
    with (destination/'preparation.log').open('x') as log:
        subprocess.run(command, env=env, stdout=log, stderr=subprocess.STDOUT, check=True, timeout=180)
    case = Path(spec['attempt'])
    plan = json.loads((case/'plan.json').read_text())
    if (plan['tasks'] != [report['task']] or plan['experimental_method_id'] != spec['method_id']
            or plan['electronic_review'] != str(parent_path)
            or sha256_file(case/'input.xyz') != plan['scored_geometry_sha256']
            or plan['maximum_host_RSS_GiB'] != 80 or plan['maximum_peak_GPU_MiB'] != 31744):
        raise ValueError('Archived native preparation changed exact task, model or resource class')
    spec.update(status='derived_native_prepared_exact_electronic_parent',
                plan_sha256=sha256_file(case/'plan.json'), geometry_sha256=plan['scored_geometry_sha256'],
                portable_parent=str(parent_path), portable_parent_sha256=sha256_file(parent_path))
    atomic_write_json(destination/'namespace.json', spec)
    return spec


def probe(path):
    spec = verify_spec(path)
    prior = Path(spec['prior_root'])
    code = ('import hashlib,json,pathlib; c=pathlib.Path('+repr(spec['original_campaign'])+'); '
            'p=json.loads((c/"plan.json").read_text()); m={**p["producer_sha256"],**p["input_sha256"]}; '
            'bad=[k for k,v in m.items() if hashlib.sha256(pathlib.Path(k).read_bytes()).hexdigest()!=v]; '
            'assert not bad,bad; print(json.dumps({"bindings":len(m)}))')
    result = subprocess.run(namespace_command(spec)+['/usr/bin/python3', '-I', '-S', '-c', code],
                            capture_output=True, text=True, check=True, timeout=180)
    record = json.loads(result.stdout)
    profile = json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text())
    code = ('import importlib.metadata as m,json;import ase,rdkit,numpy,psutil;'
            'wanted='+repr(profile['versions'])+';actual={k:m.version(k) for k in wanted};'
            'assert actual==wanted;print(json.dumps(actual))')
    result = subprocess.run(namespace_command(spec)+[profile['python'], '-c', code],
                            capture_output=True, text=True, check=True, timeout=60)
    record.update(status='derived_native_original_paths_and_CPU_dependencies_checked',
                  versions=json.loads(result.stdout), host=os.uname().nodename.split('.')[0],
                  namespace_sha256=sha256_file(path), native_G_admitted=False)
    atomic_write_json(path.with_name('namespace_probe.json'), record)
    return record


def launch(path, gpu, destination):
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity, gpu_suitability, free_slots, GATE
    spec = verify_spec(path)
    prior, case = Path(spec['prior_root']), Path(spec['original_campaign'])
    proof = json.loads(path.with_name('namespace_probe.json').read_text())
    host = os.uname().nodename.split('.')[0]
    if proof['host'] != host or proof['namespace_sha256'] != sha256_file(path):
        raise ValueError('Exact node-local namespace and dependency probe required')
    if destination.exists() or (Path(spec['attempt'])/'residency_guard').exists():
        raise FileExistsError('Preserve earlier scientific and placement attempts')
    destination.mkdir(parents=True)
    profile = json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text())
    cmd = namespace_command(spec)
    locks = GATE.parent.parent/'.gpu_locks'
    cmd[-1:-1] = ['--bind', str(locks), str(locks)]
    cmd += [str(GATE), 'run', str(gpu), '--label', 'redox-derived-native', '--', profile['python'],
            str(prior/'scripts/fleet/exclusive_python_worker.py'), '--guard-out', str(case/'residency_guard'),
            str(prior/'scripts/fleet/bounded_host_python_worker.py'), '--out', str(case/'host_guard'),
            '--maximum-rss-gib', '80', '--minimum-free-gib', '32',
            str(prior/'scripts/fleet/tensor_runtime_worker.py'), '--profile', str(prior/'environments/cutensor22/runtime.json'),
            '--out', str(case/'tensor_runtime'), str(prior/'scripts/pipeline/smooth_polymer_native.py'),
            'run', '--campaign', str(case)]
    env = dict(os.environ, PYTHONPATH=str(prior/'src'), LD_LIBRARY_PATH=profile['library_path'],
               PYTHONDONTWRITEBYTECODE='1', CUPY_CACHE_DIR='/tmp/cupy_kernel_cache')
    for key in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'XTB_THREADS'):
        env[key] = '2'
    with allocation_guard():
        capacity = host_capacity(2, minimum_ram_gib=112)
        promised = remaining_host_promises()
        if capacity['available_RAM_bytes']-promised < 112*2**30:
            raise RuntimeError('Host RAM already promised; no native calculation started')
        cpus = cpu_allocation(2)
        suitability = gpu_suitability(gpu)
        if (host, gpu) not in free_slots([host]):
            raise RuntimeError('GPU occupied or reserved; no native calculation started')
        shutil.copyfile(Path(__file__), destination/'launcher_source.py')
        with (destination/'command.log').open('x') as log:
            child = subprocess.Popen(['taskset', '--cpu-list', ','.join(map(str, cpus)), *cmd], env=env,
                stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True,
                preexec_fn=lambda: os.sched_setaffinity(0, cpus))
        record = dict(status='derived_native_submitted_scientific_start_unverified', host=host, GPU=gpu,
                      pid=child.pid, process_start=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')', 1)[1].split()[19],
                      CPU_affinity=cpus, host_capacity=capacity, remaining_RAM_promises_bytes=promised,
                      gpu_preflight=suitability, argv=cmd, namespace=str(path), namespace_sha256=sha256_file(path),
                      producer_sha256=sha256_file(destination/'launcher_source.py'), started=time.time(),
                      native_G_admitted=False)
        atomic_write_json(destination/'launch.json', record)
    return child, record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('prepare', 'probe', 'launch'))
    p.add_argument('--energy-spec', type=Path)
    p.add_argument('--spec', type=Path)
    p.add_argument('--out', type=Path)
    p.add_argument('--gpu', type=int)
    a = p.parse_args()
    if a.action == 'prepare':
        if a.energy_spec is None or a.out is None:
            p.error('Exact checked energy namespace and fresh destination required')
        result = prepare(a.energy_spec.resolve(), a.out.resolve())
    elif a.action == 'probe':
        result = probe(a.spec.resolve())
    else:
        if a.gpu is None or a.out is None:
            p.error('Explicit GPU and fresh placement directory required')
        child, result = launch(a.spec.resolve(), a.gpu, a.out.resolve())
        print(json.dumps(result), flush=True)
        result.update(returncode=child.wait(), finished=time.time(),
                      status='derived_native_process_returned_independent_review_required')
        atomic_write_json(a.out.resolve()/'completion.json', result)
        return
    print(json.dumps(result))


if __name__ == '__main__':
    main()
