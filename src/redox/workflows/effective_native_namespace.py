"""Prospective smaller-RSS execution with the original native equations intact.

The original 80-GiB same-PID guard stays active. An outer, stricter guard defines
the effective reservation. Only the resource-admission floor is adapted; native
geometry, SCF, Hessian, mode and thermochemistry checks remain unchanged.
"""
import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time

from rdkit import Chem

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.archived_sources import ArchivedSources
from redox.validation.namespace_sources import NamespaceSources, read_checked_review
from redox.workflows.archived_native_namespace import mount_command, remaining_host_promises, verify_spec
from redox.workflows.derived_native_namespace import namespace_command


def constitution(task):
    molecule = Chem.MolFromSmiles(task['smiles'])
    if molecule is None:
        raise ValueError('Invalid resource-class molecular graph')
    Chem.RemoveStereochemistry(molecule)
    return Chem.MolToSmiles(molecule, isomericSmiles=False)


def validate_class(target, reference, report, bound):
    if (not math.isfinite(bound) or not 0 < bound <= 80
            or target['maximum_host_RSS_GiB'] != 80
            or report['status'] != 'alternate_native_state_independently_reconstructed'
            or report['primary_candidate_admission']
            or target['experimental_method_id'] != report['experimental_method_id']):
        raise ValueError('A checked native reference and stronger effective bound are required')
    for name in ('experimental_native_settings', 'experimental_model', 'quantum_settings', 'runtime_profile'):
        if target[name] != reference[name]:
            raise ValueError('Reference native numerical settings or runtime differ')
    a, b = target['tasks'][0], report['task']
    if (len(target['tasks']) != 1 or reference['tasks'] != [b]
            or (a['charge'], a['mult']) != (0, 1)
            or any(a[k] != b[k] for k in ('atoms', 'active_basis_functions', 'charge', 'mult'))
            or constitution(a) != constitution(b)):
        raise ValueError('Measured resource class has another constitution/state/basis')
    peak = report['peak_host_RSS_bytes']/2**30
    if not 0 < peak*1.25 <= bound or peak+8 > bound:
        raise ValueError('Effective bound needs 25 percent and 8 GiB measured headroom')
    return dict(reference_observed_peak_GiB=peak, effective_bound_GiB=bound,
                sampled_reference_is_not_a_hard_peak=True, profile_qualified=False)


def mount(profile, spec):
    command = mount_command(spec) if profile['layout'] == 'archived' else namespace_command(spec)
    binary = ROOT.parent/'redox_polymer_continuations/runtime/bin/bwrap'
    if binary.is_file():
        command[0] = str(binary)
    return command


def prepare(namespace, reference_path, expected_digest, bound, layout, output, reference_namespace=None):
    spec = verify_spec(namespace)
    if sha256_file(reference_path) != expected_digest:
        raise ValueError('Native resource witness differs from its recorded digest')
    prior, archive = Path(spec['prior_root']), Path(spec['archive'])
    if reference_namespace is None:
        reader = ArchivedSources(prior, archive, relative_to_prior=True)
        report = json.loads(reference_path.read_text())
    else:
        report, checked_binding = read_checked_review(reference_namespace, 'native')
        if checked_binding['review'] != str(reference_path):
            raise ValueError('Native resource namespace has a different checked review')
        reader = NamespaceSources(json.loads(reference_namespace.read_text()))
    reader.verify_all(report['source_sha256'])
    ref_path = str(Path(report['campaign'])/'plan.json')
    reference = json.loads(reader.verify(ref_path, report['source_sha256'][ref_path]).read_text())
    target = json.loads((Path(spec['attempt'])/'plan.json').read_text())
    original = str(prior/'scripts/pipeline/smooth_polymer_native.py')
    if (target['execution_driver'] != original or target['primary_candidate_admission']
            or target['producer_sha256'][original] != reference['producer_sha256'][original]):
        raise ValueError('Only the identical original native scientific worker is supported')
    witness = validate_class(target, reference, report, bound)
    guard_path = str(Path(report['campaign'])/'host_guard/status.json')
    recorded_guard = json.loads(reader.verify(guard_path, report['source_sha256'][guard_path]).read_text())
    if recorded_guard['peak_RSS_bytes'] != report['peak_host_RSS_bytes']:
        raise ValueError('Reported reference RSS differs from its bound monitor')
    if (output.exists() or output.resolve().is_relative_to(prior.resolve())
            or output.resolve().is_relative_to(archive.resolve())):
        raise ValueError('Fresh execution profile outside source checkouts required')
    output.mkdir(parents=True)
    files = {}
    for source, name in ((ROOT/'scripts/fleet/bounded_peak_host_python_worker.py', 'peak_guard.py'),
                         (ROOT/'scripts/pipeline/run_native_effective_budget.py', 'resource_worker.py'),
                         (Path(__file__), 'launcher_source.py')):
        destination = output/name
        shutil.copyfile(source, destination)
        files[str(destination)] = sha256_file(destination)
    files[original] = target['producer_sha256'][original]
    reader.recheck()
    profile = dict(status='prospective_effective_native_memory_profile', native_namespace=str(namespace),
                   native_namespace_sha256=sha256_file(namespace), native_plan_sha256=spec['plan_sha256'],
                   native_settings=target['experimental_native_settings'], original_native_worker=original,
                   original_inner_guard_GiB=80, effective_bound_GiB=bound, minimum_free_GiB=32,
                   reference_review=str(reference_path), reference_review_sha256=expected_digest,
                   reference_namespace=None if reference_namespace is None else str(reference_namespace),
                   resource_derivation=witness, worker_source_sha256=files, layout=layout,
                   scientific_equations_or_settings_changed=False, profile_qualified=False)
    atomic_write_json(output/'profile.json', profile)
    return profile


def verify(path):
    profile = json.loads(path.read_text())
    namespace = Path(profile['native_namespace'])
    if sha256_file(namespace) != profile['native_namespace_sha256']:
        raise ValueError('Native namespace changed after resource-profile preparation')
    spec = verify_spec(namespace)
    reader = NamespaceSources(spec)
    reader.verify_all(profile['worker_source_sha256'])
    reader.recheck()
    return profile, spec


def probe(path):
    profile, spec = verify(path)
    runtime = json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text())
    code = ('import json,pathlib,hashlib,importlib.metadata as m; '
            'import ase,rdkit,numpy,psutil; '
            'c=pathlib.Path('+repr(spec['original_campaign'])+');p=json.loads((c/"plan.json").read_text());'
            's={**p["producer_sha256"],**p["input_sha256"]};'
            'assert all(hashlib.sha256(pathlib.Path(k).read_bytes()).hexdigest()==v for k,v in s.items());'
            'wanted='+repr(runtime['versions'])+';assert {k:m.version(k) for k in wanted}==wanted;'
            'print(json.dumps(dict(bindings=len(s),versions=wanted)))')
    result = subprocess.run(mount(profile, spec)+[runtime['python'], '-c', code],
                            capture_output=True, text=True, check=True, timeout=180)
    proof = dict(status='effective_native_node_sources_checked', host=os.uname().nodename.split('.')[0],
                 profile_sha256=sha256_file(path), result=json.loads(result.stdout), GPU_started=False)
    atomic_write_json(path.with_name('node_probe.json'), proof)
    return proof


def launch(path, gpu, output):
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity, gpu_suitability, free_slots, GATE
    profile, spec = verify(path)
    proof = json.loads(path.with_name('node_probe.json').read_text())
    host = os.uname().nodename.split('.')[0]
    if proof['host'] != host or proof['profile_sha256'] != sha256_file(path):
        raise ValueError('Exact node-local profile/source probe required')
    case, actual = Path(spec['original_campaign']), Path(spec['attempt'])
    if output.exists() or (actual/'residency_guard').exists() or (actual/'effective_host_guard').exists():
        raise FileExistsError('Preserve existing scientific and execution attempts')
    output.mkdir(parents=True)
    prior, bound = Path(spec['prior_root']), profile['effective_bound_GiB']
    runtime = json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text())
    command = mount(profile, spec)
    locks = GATE.parent.parent/'.gpu_locks'
    command[-1:-1] = ['--bind', str(locks), str(locks)]
    command += [str(GATE), 'run', str(gpu), '--label', 'redox-effective-native', '--', runtime['python'],
        str(prior/'scripts/fleet/exclusive_python_worker.py'), '--guard-out', str(case/'residency_guard'),
        str(path.parent/'peak_guard.py'), '--out', str(case/'effective_host_guard'),
        '--maximum-rss-gib', str(bound), '--minimum-free-gib', '32',
        str(prior/'scripts/fleet/bounded_host_python_worker.py'), '--out', str(case/'host_guard'),
        '--maximum-rss-gib', '80', '--minimum-free-gib', '32',
        str(prior/'scripts/fleet/tensor_runtime_worker.py'), '--profile', str(prior/'environments/cutensor22/runtime.json'),
        '--out', str(case/'tensor_runtime'), str(path.parent/'resource_worker.py'),
        '--profile', str(path), '--campaign', str(case)]
    env = dict(os.environ, PYTHONPATH=str(prior/'src'), LD_LIBRARY_PATH=runtime['library_path'],
               PYTHONDONTWRITEBYTECODE='1', CUPY_CACHE_DIR='/tmp/cupy_kernel_cache')
    for name in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'XTB_THREADS'):
        env[name] = '2'
    with allocation_guard():
        capacity = host_capacity(2, minimum_ram_gib=bound+32)
        promises = remaining_host_promises()
        if capacity['available_RAM_bytes']-promises < (bound+32)*2**30:
            raise RuntimeError('Effective native memory already promised; no calculation started')
        cpus, suitable = cpu_allocation(2), gpu_suitability(gpu)
        if (host, gpu) not in free_slots([host]):
            raise RuntimeError('GPU occupied or reserved; no calculation started')
        with (output/'command.log').open('x') as log:
            child = subprocess.Popen(['taskset', '--cpu-list', ','.join(map(str, cpus)), *command], env=env,
                stdout=log, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL, start_new_session=True,
                preexec_fn=lambda: os.sched_setaffinity(0, cpus))
        record = dict(status='effective_native_submitted_start_unverified', host=host, GPU=gpu, pid=child.pid,
                      process_start=Path(f'/proc/{child.pid}/stat').read_text().rsplit(')', 1)[1].split()[19],
                      CPU_affinity=cpus, resource_preflight=capacity, remaining_RAM_promises_bytes=promises,
                      GPU_preflight=suitable, profile=str(path), profile_sha256=sha256_file(path), argv=command,
                      effective_bound_GiB=bound, inner_guard_GiB=80, started=time.time(), native_G_qualified=False)
        atomic_write_json(output/'launch.json', record)
    return child, record


def review(path, output):
    profile, spec = verify(path)
    native, binding = read_checked_review(Path(profile['native_namespace']), 'native')
    case = Path(spec['attempt'])
    guard = json.loads((case/'effective_host_guard/status.json').read_text())
    execution = json.loads((case/'effective_execution.json').read_text())
    if (guard['status'] != 'effective_memory_guard_finished_profile_review_pending'
            or 'finished' not in guard or guard['maximum_RSS_GiB'] != profile['effective_bound_GiB']
            or max(guard['peak_tree_RSS_bytes'], guard['self_VmHWM_bytes']) > profile['effective_bound_GiB']*2**30
            or execution['native_numerical_settings_changed'] or execution['geometry_or_energy_changed']
            or execution['capacity_calls'] != [dict(original_floor_GiB=112,
                    effective_floor_GiB=profile['effective_bound_GiB']+32)]):
        raise ValueError('Effective resource profile is not qualified')
    result = dict(status='effective_native_resource_profile_checked', profile_sha256=sha256_file(path),
                  normal_native_binding=binding, effective_bound_GiB=profile['effective_bound_GiB'],
                  peak_tree_RSS_bytes=guard['peak_tree_RSS_bytes'], self_VmHWM_bytes=guard['self_VmHWM_bytes'],
                  G_by_cutoff_eV=native['G_by_cutoff_eV'], tested_run_profile_qualified=True,
                  experimental_accuracy_inferred=False)
    with output.open('x') as stream:
        json.dump(result, stream, indent=2)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('prepare', 'probe', 'launch', 'review'))
    for name in ('namespace', 'reference-review', 'profile', 'out'):
        p.add_argument('--'+name, type=Path)
    p.add_argument('--reference-digest')
    p.add_argument('--reference-namespace', type=Path)
    p.add_argument('--bound', type=float)
    p.add_argument('--layout', choices=('archived', 'derived'), default='derived')
    p.add_argument('--gpu', type=int)
    a = p.parse_args()
    if a.action == 'prepare':
        result = prepare(a.namespace.resolve(), a.reference_review.resolve(), a.reference_digest,
                         a.bound, a.layout, a.out.resolve(),
                         None if a.reference_namespace is None else a.reference_namespace.resolve())
    elif a.action == 'probe':
        result = probe(a.profile.resolve())
    elif a.action == 'review':
        result = review(a.profile.resolve(), a.out.resolve())
    else:
        child, result = launch(a.profile.resolve(), a.gpu, a.out.resolve())
        print(json.dumps(result), flush=True)
        result.update(returncode=child.wait(), status='effective_native_returned_independent_reviews_required')
        atomic_write_json(a.out/'completion.json', result)
    print(json.dumps(result))


if __name__ == '__main__':
    main()
