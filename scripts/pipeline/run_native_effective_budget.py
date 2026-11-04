#!/usr/bin/env python3
"""Run the unchanged native equations under an additional stricter RSS guard."""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

from redox.core.io import atomic_write_json
from redox.workflows import resources


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def capacity_adapter(original, bound, calls):
    def capacity(threads, **kwargs):
        if threads == 2 and kwargs.get('minimum_ram_gib') == 112:
            calls.append(dict(original_floor_GiB=112, effective_floor_GiB=bound+32))
            kwargs['minimum_ram_gib'] = bound+32
        return original(threads, **kwargs)
    return capacity


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--profile', type=Path, required=True)
    p.add_argument('--campaign', type=Path, required=True)
    a = p.parse_args()
    profile = json.loads(a.profile.read_text())
    if any(sha(path) != digest for path, digest in profile['worker_source_sha256'].items()):
        raise ValueError('Effective runtime source changed')
    plan = json.loads((a.campaign/'plan.json').read_text())
    if (sha(a.campaign/'plan.json') != profile['native_plan_sha256']
            or plan['experimental_native_settings'] != profile['native_settings']
            or not 0 < profile['effective_bound_GiB'] <= plan['maximum_host_RSS_GiB'] == 80):
        raise ValueError('Effective execution changed the exact native plan or settings')
    guard = json.loads((a.campaign/'effective_host_guard/status.json').read_text())
    if guard['maximum_RSS_GiB'] != profile['effective_bound_GiB'] or guard['minimum_free_GiB'] != 32:
        raise ValueError('The declared effective memory guard must be active')
    original = resources.host_capacity
    calls = []

    resources.host_capacity = capacity_adapter(original, profile['effective_bound_GiB'], calls)
    source = Path(profile['original_native_worker'])
    try:
        specification = importlib.util.spec_from_file_location('unchanged_native_scientific_worker', source)
        module = importlib.util.module_from_spec(specification)
        specification.loader.exec_module(module)
        module.run(a.campaign)
        if len(calls) != 1:
            raise ValueError('Native resource-admission adaptation did not match its declared scope')
    finally:
        resources.host_capacity = original
        atomic_write_json(a.campaign/'effective_execution.json', dict(
            status='resource_admission_adapter_finished', profile=str(a.profile), profile_sha256=sha(a.profile),
            original_worker=str(source), original_worker_sha256=sha(source), capacity_calls=calls,
            native_numerical_settings_changed=False, geometry_or_energy_changed=False,
            inner_native_guard_GiB=80, effective_outer_guard_GiB=profile['effective_bound_GiB']))


if __name__ == '__main__':
    main()
