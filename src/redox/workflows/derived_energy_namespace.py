"""Measured-class electronic successors using unchanged archived QM workers."""
import argparse
from contextlib import contextmanager
import importlib.util
import json
import os
from pathlib import Path
import shutil
import subprocess
import threading
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.namespace_sources import NamespaceSources, read_checked_review
from redox.workflows import archived_energy_namespace as energy
from redox.workflows.derived_native_namespace import namespace_command, portable_parent
from redox.workflows.archived_native_namespace import verify_spec


def mount(spec):
    command = namespace_command(spec)
    portable = ROOT.parent/'redox_polymer_continuations/runtime/bin/bwrap'
    if portable.is_file():
        command[0] = str(portable)
    return command


def measured_validator():
    source = ROOT/'scripts/pipeline/run_measured_energy_followup.py'
    specification = importlib.util.spec_from_file_location('measured_energy_exact_resource_policy', source)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module.validate_measured


@contextmanager
def policy(archive, spec=None):
    if threading.current_thread() is not threading.main_thread() or threading.active_count() != 1:
        raise RuntimeError('Install filesystem/resource adapters before starting threads')
    original_mount, original_validation = energy.energy_mount, energy.validate_plan
    validate = measured_validator()
    energy.energy_mount = mount
    energy.validate_plan = lambda plan, prior: validate(
        plan, prior, archive, source_reader=None if spec is None else NamespaceSources(spec))
    try:
        yield
    finally:
        energy.energy_mount, energy.validate_plan = original_mount, original_validation


def prepare(configuration_path, destination):
    cfg = json.loads(configuration_path.read_text())
    report, binding = read_checked_review(Path(cfg['energy_namespace']), 'energy')
    source = json.loads(Path(cfg['energy_namespace']).read_text())
    prior, archive = Path(source['prior_root']), Path(source['archive'])
    if (destination.exists() or destination.resolve().is_relative_to(prior.resolve())
            or destination.resolve().is_relative_to(archive.resolve())):
        raise ValueError('Fresh derived outputs outside both source checkouts required')
    label = cfg['class_label']
    if not label or Path(label).name != label or label in ('.', '..'):
        raise ValueError('A plain class label is required')
    destination.mkdir(parents=True)
    campaign_parent, raw = destination/'campaign_parent', destination/'raw'
    campaign_parent.mkdir()
    raw.mkdir()
    logical_parent = prior/'results/polymer_models/campaigns/2026-10-30_measured_energy_classes'/label
    scaffold = dict(prior_root=str(prior), archive=str(archive),
                    original_campaign=str(logical_parent/'placeholder'),
                    attempt=str(campaign_parent/'placeholder'),
                    original_raw_parent=str(prior/'calcs/polymer_models'/('preparation_'+label)),
                    raw_parent=str(raw))
    parent = destination/'portable_resource_review.json'
    atomic_write_json(parent, portable_parent(report, source, Path(binding['review'])))
    config_snapshot = destination/'source_configuration.json'
    shutil.copyfile(configuration_path, config_snapshot)
    driver = ROOT/'scripts/pipeline/prepare_measured_energy_class.py'
    snapshot = destination/'preparation_driver.py'
    shutil.copyfile(driver, snapshot)
    input_cfg = destination/'preparation_inputs.json'
    inner = dict(cfg, prior_root=str(prior), logical_output=str(logical_parent/label),
                 reference=source['attempt'], resource_review=str(parent), primary_candidate_admission=False,
                 source_sha256={str(config_snapshot): sha256_file(config_snapshot),
                                str(snapshot): sha256_file(snapshot), str(parent): sha256_file(parent),
                                str(Path(cfg['energy_namespace'])): sha256_file(cfg['energy_namespace'])})
    atomic_write_json(input_cfg, inner)
    profile = json.loads((archive/'environments/cutensor22/runtime.json').read_text())
    command = mount(scaffold)+['env', 'PYTHONPATH='+str(prior/'src'), profile['python'], str(snapshot),
                               '--configuration', str(input_cfg)]
    env = dict(os.environ, PYTHONDONTWRITEBYTECODE='1')
    for key in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'XTB_THREADS'):
        env[key] = '1'
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4)
        cpus = cpu_allocation(1)
    with (destination/'preparation.log').open('x') as log:
        subprocess.run(['taskset', '--cpu-list', ','.join(map(str, cpus)), *command], env=env,
                       stdout=log, stderr=subprocess.STDOUT, check=True, timeout=1800)
    result = json.loads((campaign_parent/label/'preparation.json').read_text())
    return finalize(configuration_path, destination, result, binding)


def finalize(configuration_path, destination, result, binding=None):
    cfg = json.loads(configuration_path.read_text())
    source = json.loads(Path(cfg['energy_namespace']).read_text())
    prior, archive = Path(source['prior_root']), Path(source['archive'])
    if binding is None:
        _, binding = read_checked_review(Path(cfg['energy_namespace']), 'energy')
    campaign_parent, raw = destination/'campaign_parent', destination/'raw'
    logical_parent = prior/'results/polymer_models/campaigns/2026-10-30_measured_energy_classes'/cfg['class_label']
    namespaces = []
    reuse = {}
    for candidate in cfg.get('prepared_namespace_reuse', []):
        prior_spec = verify_spec(Path(candidate))
        if prior_spec['key'] in reuse:
            raise ValueError('Prepared namespace reuse repeats a conformer identity')
        reuse[prior_spec['key']] = (candidate, prior_spec)
    for entry in result['entries']:
        if entry['status'] != 'prepared':
            continue
        original = Path(entry['campaign'])
        attempt = campaign_parent/original.relative_to(logical_parent)
        plan = json.loads((attempt/'plan.json').read_text())
        raw_parent = raw/Path(plan['raw_root']).parent.name
        raw_parent.mkdir(exist_ok=True)
        spec = dict(status='derived_measured_energy_prepared', archive=str(archive), prior_root=str(prior),
                    key=entry['key'], original_campaign=str(original), attempt=str(attempt),
                    original_raw_parent=str(Path(plan['raw_root']).parent), raw_parent=str(raw_parent),
                    original_raw_root=plan['raw_root'], method_id=plan['experimental_method_id'],
                    plan_sha256=sha256_file(attempt/'plan.json'), geometry_sha256=plan['scored_geometry_sha256'],
                    maximum_host_RSS_GiB=plan['maximum_host_RSS_GiB'], native_G_computed=False,
                    primary_candidate_admission=False, resource_binding=binding)
        measured_validator()(plan, prior, archive, source_reader=NamespaceSources(spec))
        if entry['key'] in reuse:
            path, previous_spec = reuse[entry['key']]
            previous_plan = json.loads((Path(previous_spec['attempt'])/'plan.json').read_text())
            exact = ('tasks', 'scored_geometry_sha256', 'quantum_settings', 'runtime_profile',
                     'experimental_method_id', 'experimental_model', 'maximum_host_RSS_GiB')
            if any(previous_plan[name] != plan[name] for name in exact):
                raise ValueError('Prepared namespace reuse changes the exact scored state or resource class')
            measured_validator()(previous_plan, prior, archive, source_reader=NamespaceSources(previous_spec))
            namespaces.append(path)
            entry['reused_preparation_namespace'] = path
            entry['reused_preparation_namespace_sha256'] = sha256_file(path)
            continue
        folder = destination/'namespaces'/entry['key']
        folder.mkdir(parents=True, exist_ok=True)
        path = folder/'namespace.json'
        if path.exists():
            if json.loads(path.read_text()) != spec:
                raise ValueError('Existing namespace differs; preserve active and completed states')
        else:
            atomic_write_json(path, spec)
        namespaces.append(str(path))
    result.update(namespaces=namespaces, source_configuration=str(configuration_path),
                  source_configuration_sha256=sha256_file(configuration_path), resource_binding=binding)
    previous = destination/'manifest.json'
    if previous.exists():
        history = destination/'manifest_history'
        history.mkdir(exist_ok=True)
        (history/(str(time.time_ns())+'.json')).write_bytes(previous.read_bytes())
    atomic_write_json(previous, result)
    return result


def refresh(destination):
    previous = json.loads((destination/'manifest.json').read_text())
    config = Path(previous['source_configuration'])
    if sha256_file(config) != previous['source_configuration_sha256']:
        raise ValueError('Class configuration changed; preserve the prior preparation')
    cfg = json.loads(config.read_text())
    source = json.loads(Path(cfg['energy_namespace']).read_text())
    prior, archive = Path(source['prior_root']), Path(source['archive'])
    parent = prior/'results/polymer_models/campaigns/2026-10-30_measured_energy_classes'/cfg['class_label']
    scaffold = dict(prior_root=str(prior), archive=str(archive), original_campaign=str(parent/'placeholder'),
                    attempt=str(destination/'campaign_parent/placeholder'),
                    original_raw_parent=str(prior/'calcs/polymer_models'/('preparation_'+cfg['class_label'])),
                    raw_parent=str(destination/'raw'))
    run = destination/'dependency_refreshes'/str(time.time_ns())
    run.mkdir(parents=True)
    script = run/'driver.py'
    shutil.copyfile(ROOT/'scripts/pipeline/prepare_measured_energy_class.py', script)
    profile = json.loads((archive/'environments/cutensor22/runtime.json').read_text())
    command = mount(scaffold)+['env', 'PYTHONPATH='+str(prior/'src'), profile['python'], str(script),
                               '--configuration', str(destination/'preparation_inputs.json'), '--refresh']
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4)
        cpus = cpu_allocation(1)
    with (run/'refresh.log').open('x') as log:
        subprocess.run(['taskset', '--cpu-list', ','.join(map(str, cpus)), *command],
                       stdout=log, stderr=subprocess.STDOUT, check=True, timeout=1800)
    result = json.loads((destination/'campaign_parent'/cfg['class_label']/'preparation.json').read_text())
    return finalize(config, destination, result)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('prepare', 'refresh', 'probe', 'launch'))
    for name in ('configuration', 'spec', 'out'):
        p.add_argument('--'+name, type=Path)
    p.add_argument('--gpu', type=int)
    a = p.parse_args()
    if a.action in ('prepare', 'refresh'):
        result = (prepare(a.configuration.resolve(), a.out.resolve()) if a.action == 'prepare'
                  else refresh(a.out.resolve()))
        print(json.dumps(dict(status=result['status'], namespaces=len(result['namespaces']))))
        return
    spec = json.loads(a.spec.read_text())
    with policy(Path(spec['archive']), spec):
        if a.action == 'probe':
            print(json.dumps(energy.probe(a.spec.resolve())))
            return
        child, result = energy.launch(a.spec.resolve(), a.gpu, a.out.resolve())
    for source in (Path(__file__), ROOT/'scripts/pipeline/run_measured_energy_followup.py'):
        shutil.copyfile(source, a.out/source.name)
    result['filesystem_and_resource_adapters_sha256'] = {
        str(a.out/source.name): sha256_file(a.out/source.name) for source in
        (Path(__file__), ROOT/'scripts/pipeline/run_measured_energy_followup.py')}
    atomic_write_json(a.out/'launch.json', result)
    print(json.dumps(result), flush=True)
    result.update(returncode=child.wait(), status='derived_energy_returned_independent_review_required')
    atomic_write_json(a.out/'completion.json', result)


if __name__ == '__main__':
    main()
