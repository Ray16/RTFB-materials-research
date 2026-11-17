#!/usr/bin/env python3
"""Retain native numerical gates and validate the extra execution-wrapper chain.

Only the runtime program-identity predicate is extended. The original SCF,
derivative, geometry, mode, partition and residency tests execute unchanged.
Original runtime records are read, never relabeled or rewritten.
"""
import argparse
import hashlib
import importlib.util
import inspect
import json
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_execution(runtime, plan, campaign, profile_path):
    profile = json.loads(profile_path.read_text())
    if (sha(campaign/'plan.json') != profile['native_plan_sha256']
            or profile['native_settings'] != plan['experimental_native_settings']
            or profile['original_native_worker'] != plan['execution_driver']
            or any(sha(path) != digest for path, digest in profile['worker_source_sha256'].items())):
        raise ValueError('Wrapped runtime has a different scientific plan or source')
    wrapper = str(profile_path.parent/'resource_worker.py')
    if runtime['program'] not in (wrapper, profile['original_native_worker']):
        raise ValueError('Unknown native execution program')
    execution = json.loads((campaign/'effective_execution.json').read_text())
    guard = json.loads((campaign/'effective_host_guard/status.json').read_text())
    if (execution['profile_sha256'] != sha(profile_path)
            or execution['original_worker'] != profile['original_native_worker']
            or execution['original_worker_sha256'] != profile['worker_source_sha256'][profile['original_native_worker']]
            or execution['native_numerical_settings_changed'] or execution['geometry_or_energy_changed']
            or execution['capacity_calls'] != [dict(original_floor_GiB=112,
                    effective_floor_GiB=profile['effective_bound_GiB']+32)]
            or guard['status'] != 'effective_memory_guard_finished_profile_review_pending'
            or 'finished' not in guard or guard['maximum_RSS_GiB'] != profile['effective_bound_GiB']
            or guard['minimum_free_GiB'] != 32
            or not 0 < max(guard['peak_tree_RSS_bytes'], guard['self_VmHWM_bytes']) <= profile['effective_bound_GiB']*2**30):
        raise ValueError('Effective native execution/resource chain failed')
    return True


def adapt_source(body):
    anchor = "rt['program']!=plan['execution_driver']"
    replacement = 'not _checked_execution_program(rt,plan,campaign)'
    if body.count(anchor) != 1:
        raise ValueError('Original native runtime identity predicate changed')
    result = body.replace(anchor, replacement, 1)
    if result.replace(replacement, anchor, 1) != body:
        raise ValueError('Native numerical review body changed')
    return result


def review(campaign, profile_path, output):
    if output.exists():
        raise FileExistsError('Preserve existing independent native review')
    source = ROOT/'scripts/validation/review_smooth_polymer_native.py'
    specification = importlib.util.spec_from_file_location('original_native_numerical_review', source)
    original = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(original)
    body = inspect.getsource(original.review)
    adapted = adapt_source(body)
    folder = campaign/'review_sources'/str(time.time_ns())
    folder.mkdir(parents=True)
    snapshot = folder/'effective_review_function.py'
    snapshot.write_text(adapted)
    scope = dict(vars(original), _checked_execution_program=lambda rt, plan, case:
                 validate_execution(rt, plan, case, profile_path))
    exec(compile(adapted, str(snapshot), 'exec'), scope)
    temporary = folder/'numerical_review.json'
    result = scope['review'](campaign, temporary)
    if result is None:
        result = json.loads(temporary.read_text())
    for path in (Path(__file__).resolve(), snapshot, profile_path, campaign/'effective_execution.json',
                 campaign/'effective_host_guard/status.json'):
        result['source_sha256'][str(path)] = sha(path)
    result['effective_runtime_validation'] = dict(profile=str(profile_path), profile_sha256=sha(profile_path),
        actual_program=json.loads((campaign/'tensor_runtime/status.json').read_text())['program'],
        original_numerical_reviewer=str(source), original_reviewer_sha256=sha(source),
        adapted_review_source=str(snapshot), adapted_review_sha256=sha(snapshot),
        numerical_checks_or_thresholds_changed=False, original_runtime_records_changed=False)
    atomic_write_json(output, result)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--campaign', type=Path, required=True)
    p.add_argument('--profile', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    value = review(a.campaign.resolve(), a.profile.resolve(), a.out.resolve())
    print(json.dumps(dict(status=value['status'], G_eV=value['G_eV'])))
