"""Start declared hybrid-energy RAM experiments without rewriting archived plans."""
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
from redox.validation.archived_sources import ArchivedSources
from redox.workflows.archived_native_namespace import mount_command, verify_spec, remaining_host_promises


def energy_mount(spec):
    """Create absent raw-batch mountpoints in RAM, never in the RO archive."""
    command = mount_command(spec)
    target = str(Path(spec['prior_root'])/'calcs/polymer_models')
    archive = Path(spec['archive'])/'calcs/polymer_models'
    index = command.index('--bind', command.index('--bind')+1)
    overlays = ['--tmpfs', target]
    for child in sorted(archive.iterdir()):
        if child.is_dir():
            overlays += ['--ro-bind', str(child), str(Path(target)/child.name)]
    command[index:index] = overlays
    return command


def probe(path):
    spec = verify_spec(path)
    code = ('import hashlib,json,pathlib; c=pathlib.Path('+repr(spec['original_campaign'])+'); '
            'p=json.loads((c/"plan.json").read_text()); '
            'm={**p["producer_sha256"],**p["input_sha256"]}; '
            'bad=[k for k,v in m.items() if hashlib.sha256(pathlib.Path(k).read_bytes()).hexdigest()!=v]; '
            'assert not bad,bad; print(json.dumps({"bindings":len(m)}))')
    result = subprocess.run(energy_mount(spec)+['/usr/bin/python3', '-I', '-S', '-c', code],
                            capture_output=True, text=True, check=True, timeout=180)
    report = json.loads(result.stdout)
    profile = json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text())
    code = ('import importlib.metadata as m,json; import ase,rdkit,numpy,psutil; '
            'wanted='+repr(profile['versions'])+'; actual={k:m.version(k) for k in wanted}; '
            'assert wanted==actual; print(json.dumps(actual))')
    result = subprocess.run(energy_mount(spec)+[profile['python'], '-c', code],
                            capture_output=True, text=True, check=True, timeout=60)
    report.update(status='original_path_energy_namespace_hashes_match', host=os.uname().nodename.split('.')[0],
                  namespace_sha256=sha256_file(path), versions=json.loads(result.stdout),
                  quantum_jobs_submitted=0, raw_parent_mountpoints='private tmpfs with archived siblings rebound read-only')
    atomic_write_json(path.with_name('namespace_probe.json'), report)
    return report


def validate_plan(plan, prior):
    if (len(plan['tasks']) != 1 or (plan['tasks'][0]['charge'], plan['tasks'][0]['mult']) != (0, 1)
            or plan['experimental_model']['solvent'] != 'acetonitrile'
            or plan['execution_driver'] != str(prior/'scripts/pipeline/smooth_polymer_energy_bounded.py')
            or plan['native_G_computed'] is not False or plan['primary_candidate_admission'] is not False
            or not plan['energy_only_execution'] or plan['maximum_peak_GPU_MiB'] != 31744
            or not math.isfinite(plan['maximum_host_RSS_GiB'])
            or not 80 < plan['maximum_host_RSS_GiB'] <= 128):
        raise ValueError('Only the declared alternate-model neutral-singlet bounded hybrid energy is supported')


def prepare(archive, prior, campaign, destination, inventory_path):
    inventory = json.loads(inventory_path.read_text())
    if len({r['host'] for r in inventory}) != 9 or any(r['processes'] for r in inventory):
        raise ValueError('All-host stopped-predictor inventory required')
    reader = ArchivedSources(prior, archive, relative_to_prior=True)
    source = reader.path(campaign)
    plan = json.loads((source/'plan.json').read_text())
    validate_plan(plan, prior)
    if ((source/'experimental_energy').exists() or (source/'residency_guard').exists()
            or reader.path(plan['raw_root']).exists()):
        raise ValueError('This adapter accepts only the prepared, never-started RAM experiment')
    reader.verify_all(plan['producer_sha256'])
    reader.verify_all(plan['input_sha256'])
    if sha256_file(source/'input.xyz') != plan['scored_geometry_sha256']:
        raise ValueError('Scored geometry differs')
    if destination.exists():
        raise FileExistsError('Preserve previous namespace attempts')
    case, raw = destination/'campaign', destination/'raw'
    case.mkdir(parents=True)
    raw.mkdir()
    for name in ('plan.json', 'input.xyz'):
        shutil.copyfile(source/name, case/name)
    spec = dict(status='prepared_energy_RAM_experiment_not_submitted', archive=str(archive),
                prior_root=str(prior), original_campaign=str(campaign), attempt=str(case),
                original_raw_parent=str(Path(plan['raw_root']).parent), raw_parent=str(raw),
                original_raw_root=plan['raw_root'], plan_sha256=sha256_file(source/'plan.json'),
                geometry_sha256=plan['scored_geometry_sha256'], key=plan['tasks'][0]['key'],
                worker=plan['execution_driver'], maximum_host_RSS_GiB=plan['maximum_host_RSS_GiB'],
                method_id=plan['experimental_method_id'], native_G_computed=False,
                prospective_RAM_experiment=True, historical_files_modified=False,
                inventory=str(inventory_path.resolve()), inventory_sha256=sha256_file(inventory_path))
    mount_command(spec)
    reader.recheck()
    atomic_write_json(destination/'namespace.json', spec)
    return spec


def launch(path, gpu, destination):
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity, gpu_suitability, free_slots, GATE
    spec = verify_spec(path)
    prior, case = Path(spec['prior_root']), Path(spec['original_campaign'])
    plan = json.loads((Path(spec['attempt'])/'plan.json').read_text())
    validate_plan(plan, prior)
    proof = json.loads(path.with_name('namespace_probe.json').read_text())
    host = os.uname().nodename.split('.')[0]
    if proof['host'] != host or proof['namespace_sha256'] != sha256_file(path):
        raise ValueError('Exact node-local input/dependency probe required')
    if destination.exists() or (Path(spec['attempt'])/'residency_guard').exists():
        raise FileExistsError('Preserve previous placement attempts')
    destination.mkdir(parents=True)
    bound = plan['maximum_host_RSS_GiB']
    profile = json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text())
    command = energy_mount(spec)
    locks = GATE.parent.parent/'.gpu_locks'
    command[-1:-1] = ['--bind', str(locks), str(locks)]
    command += [str(GATE), 'run', str(gpu), '--label', 'redox-archived-energy', '--', profile['python'],
                str(prior/'scripts/fleet/exclusive_python_worker.py'), '--guard-out', str(case/'residency_guard'),
                str(prior/'scripts/fleet/bounded_host_python_worker.py'), '--out', str(case/'host_guard'),
                '--maximum-rss-gib', str(bound), '--minimum-free-gib', '32',
                str(prior/'scripts/fleet/tensor_runtime_worker.py'), '--profile', str(prior/'environments/cutensor22/runtime.json'),
                '--out', str(case/'tensor_runtime'), plan['execution_driver'], 'run', '--campaign', str(case)]
    env = dict(os.environ, PYTHONPATH=str(prior/'src'), LD_LIBRARY_PATH=profile['library_path'],
               PYTHONDONTWRITEBYTECODE='1', CUPY_CACHE_DIR='/tmp/cupy_kernel_cache')
    for key in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'XTB_THREADS'):
        env[key] = '2'
    with allocation_guard():
        capacity = host_capacity(2, minimum_ram_gib=bound+32)
        promises = remaining_host_promises()
        if capacity['available_RAM_bytes']-promises < (bound+32)*2**30:
            raise RuntimeError('RAM already promised to other live workers; no scientific attempt started')
        cpus = cpu_allocation(2)
        suitability = gpu_suitability(gpu)
        if (host, gpu) not in free_slots([host]):
            raise RuntimeError('Device occupied or reserved; no scientific attempt started')
        shutil.copyfile(Path(__file__), destination/'launcher_source.py')
        with (destination/'command.log').open('x') as log:
            child = subprocess.Popen(['taskset', '--cpu-list', ','.join(map(str, cpus)), *command],
                                     env=env, stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                     start_new_session=True, preexec_fn=lambda: os.sched_setaffinity(0, cpus))
        record = dict(status='energy_namespace_started_scientific_stage_unverified', host=host, GPU=gpu,
                      pid=child.pid, process_start=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')', 1)[1].split()[19],
                      CPU_affinity=cpus, bound_GiB=bound, host_capacity=capacity, remaining_other_promises_bytes=promises,
                      gpu_preflight=suitability, argv=command, namespace=str(path), namespace_sha256=sha256_file(path),
                      source_sha256=sha256_file(destination/'launcher_source.py'), started=time.time(),
                      native_G_computed=False, energy_qualified=False)
        record['writable_compilation_cache'] = '/tmp/cupy_kernel_cache (private namespace tmpfs)'
        atomic_write_json(destination/'launch.json', record)
    return child, record


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('prepare', 'probe', 'launch'))
    p.add_argument('--spec', type=Path)
    p.add_argument('--gpu', type=int)
    p.add_argument('--out', type=Path)
    for name in ('archive', 'prior-root', 'campaign', 'destination', 'inventory'):
        p.add_argument('--'+name, type=Path)
    a = p.parse_args()
    if a.action == 'prepare':
        if any(getattr(a, name) is None for name in ('archive', 'prior_root', 'campaign', 'destination', 'inventory')):
            p.error('Explicit source roots, exact campaign, inventory and fresh destination required')
        print(json.dumps(prepare(a.archive.resolve(), a.prior_root.resolve(), a.campaign,
                                 a.destination.resolve(), a.inventory.resolve())))
    elif a.action == 'probe':
        print(json.dumps(probe(a.spec.resolve())))
    else:
        if a.gpu is None or a.out is None:
            p.error('Explicit GPU and fresh placement destination required')
        child, record = launch(a.spec.resolve(), a.gpu, a.out.resolve())
        print(json.dumps(record), flush=True)
        record.update(returncode=child.wait(), finished=time.time(),
                      status='energy_namespace_returned_independent_review_required')
        atomic_write_json(a.out.resolve()/'completion.json', record)


if __name__ == '__main__':
    main()
