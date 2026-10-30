#!/usr/bin/env python3
"""Advance independently checked energy states into resumable native inputs.

No GPU computation runs here. Existing native namespaces are reused only after
matching the exact task, model, scored geometry and checked electronic values.
Partial preparations remain visible and are never overwritten automatically.
"""
import argparse
import fcntl
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.namespace_sources import NamespaceSources, read_checked_review
from redox.validation.portable_energy_parent import stage
from redox.workflows.archived_native_namespace import verify_spec
from redox.workflows.derived_native_namespace import prepare


def match_native(namespace, electronic):
    spec = verify_spec(namespace)
    reader = NamespaceSources(spec)
    plan = json.loads((Path(spec['attempt'])/'plan.json').read_text())
    if plan['tasks'] != [electronic['task']] or plan['experimental_model'] != electronic['experimental_model']:
        raise ValueError('Existing native namespace belongs to another exact task/model')
    if plan['experimental_method_id'] != electronic['portable_method_id']:
        raise ValueError('Existing native namespace has a different solvent method')
    reader.verify_all(plan['producer_sha256'])
    reader.verify_all(plan['input_sha256'])
    parent = json.loads(reader.path(plan['electronic_review']).read_text())
    for field in ('E_smd_eV', 'E_gas_eV'):
        if abs(parent[field]-electronic[field]) > 1e-9:
            raise ValueError('Existing native state has another electronic parent value')
    expected_geometry = electronic['source_sha256'][str(Path(electronic['campaign'])/'input.xyz')]
    if plan['scored_geometry_sha256'] != expected_geometry:
        raise ValueError('Existing native state has another scored geometry')
    reader.recheck()
    return spec


def run(cohort_path, energy_status, reuse_path, root, output, watch_hours):
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4)
        os.sched_setaffinity(0, cpu_allocation(1))
    allowed = json.loads(cohort_path.read_text())['specs']
    if not allowed or len(allowed) != len(set(allowed)):
        raise ValueError('Distinct electronic source namespaces required')
    first_spec = json.loads(Path(allowed[0]).read_text())
    if any(root.resolve().is_relative_to(Path(first_spec[field]).resolve()) for field in ('prior_root', 'archive')):
        raise ValueError('Native calculations must stay outside both source checkouts')
    reuse = json.loads(reuse_path.read_text())['namespaces']
    bindings = {str(path): sha256_file(path) for path in [cohort_path, reuse_path,
        Path(__file__).resolve(), ROOT/'src/redox/validation/portable_energy_parent.py',
        ROOT/'src/redox/workflows/derived_native_namespace.py', *map(Path, allowed)]}
    root.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    with (root/'preparation.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        status_path = output/'status.json'
        previous = json.loads(status_path.read_text()) if status_path.exists() else {}
        if previous:
            history = output/'history'
            history.mkdir(exist_ok=True)
            shutil.copyfile(status_path, history/(str(time.time_ns())+'.json'))
        finished = {entry['energy_namespace']: entry for entry in previous.get('entries', [])
                    if entry['status'] in ('native_inputs_prepared', 'matching_native_namespace_reused')}
        for row in finished.values():
            if sha256_file(row['native_namespace']) != row['native_namespace_sha256']:
                raise ValueError('Cached native namespace changed; preserve prior preparation')
            verify_spec(Path(row['native_namespace']))
        held = {}
        owner = dict(host=os.uname().nodename, pid=os.getpid(),
                     process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19])
        report = dict(status='checked_energy_native_preparation_running', owner=owner, started=time.time(),
                      source_sha256=bindings, allowed_states=len(allowed), quantum_jobs_submitted=0,
                      native_free_energies_inferred=False)
        deadline = time.time()+watch_hours*3600
        failures = 0
        while True:
            if any(sha256_file(path) != digest for path, digest in bindings.items()):
                raise ValueError('Native preparation source/cohort changed')
            live = json.loads(energy_status.read_text())
            for entry in live['entries']:
                if entry['status'] != 'energy_independently_reviewed':
                    continue
                source = entry['spec']
                if source not in allowed:
                    raise ValueError('Energy collector introduced an undeclared input')
                if source in finished or source in held:
                    continue
                native_report, source_binding = read_checked_review(Path(source), 'energy')
                key = native_report['task']['key']
                folder = root/key
                namespace = Path(reuse[key]) if key in reuse else folder/'native/namespace.json'
                row = dict(key=key, energy_namespace=source, checked_energy_binding=source_binding)
                try:
                    if namespace.exists():
                        match_native(namespace, native_report)
                        row['status'] = 'matching_native_namespace_reused'
                    else:
                        parent = folder/'electronic_parent'
                        if not (parent/'namespace.json').exists():
                            if parent.exists():
                                raise FileExistsError('Partial electronic-parent preparation preserved')
                            stage(Path(source), parent)
                        if (folder/'native').exists():
                            raise FileExistsError('Partial native-input preparation preserved')
                        prepare(parent/'namespace.json', folder/'native')
                        namespace = folder/'native/namespace.json'
                        match_native(namespace, native_report)
                        row['status'] = 'native_inputs_prepared'
                    row.update(native_namespace=str(namespace), native_namespace_sha256=sha256_file(namespace))
                    finished[source] = row
                    failures = 0
                except (OSError, ValueError, subprocess.SubprocessError) as error:
                    row.update(status='preparation_failed_or_partial_preserved', error=repr(error))
                    held[source] = row
                    failures += 1
                report.update(updated=time.time(), entries=[*finished.values(), *held.values()], last_observation=row,
                              prepared_states=len(finished), consecutive_preparation_failures=failures)
                atomic_write_json(status_path, report)
                if failures >= 3:
                    report['status'] = 'preparation_circuit_breaker_existing_outputs_preserved'
                    atomic_write_json(status_path, report)
                    return
            report.update(updated=time.time(), entries=[*finished.values(), *held.values()], prepared_states=len(finished))
            atomic_write_json(status_path, report)
            if not watch_hours or time.time() >= deadline or len(finished) == len(allowed):
                report.update(status='checked_energy_native_preparation_window_complete', finished=time.time())
                atomic_write_json(status_path, report)
                return
            time.sleep(20)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cohort', type=Path, required=True)
    p.add_argument('--energy-status', type=Path, required=True)
    p.add_argument('--reuse', type=Path, required=True)
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--watch-hours', type=float, default=0)
    a = p.parse_args()
    if not 0 <= a.watch_hours <= 48:
        p.error('Use a bounded preparation window up to 48 hours')
    run(a.cohort.resolve(), a.energy_status.resolve(), a.reuse.resolve(), a.root.resolve(), a.out.resolve(), a.watch_hours)
