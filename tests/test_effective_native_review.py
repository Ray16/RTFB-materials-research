import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
specification = importlib.util.spec_from_file_location('effective_native_runtime_review',
    ROOT/'scripts/validation/review_effective_native.py')
reviewer = importlib.util.module_from_spec(specification)
specification.loader.exec_module(reviewer)


def test_only_program_identity_predicate_is_extended():
    anchor = "rt['program']!=plan['execution_driver']"
    source = 'def review():\n    assert error <= .01\n    return '+anchor+'\n'
    changed = reviewer.adapt_source(source)
    assert changed.replace('not _checked_execution_program(rt,plan,campaign)', anchor) == source
    with pytest.raises(ValueError):
        reviewer.adapt_source('another numerical reviewer')


@pytest.fixture
def runtime(tmp_path):
    native = tmp_path/'native.py'
    native.write_text('original native equations')
    wrapper = tmp_path/'resource_worker.py'
    wrapper.write_text('resource wrapper')
    plan = dict(execution_driver=str(native), experimental_native_settings=dict(gradient_limit=.00045))
    (tmp_path/'plan.json').write_text(json.dumps(plan))
    profile = dict(native_plan_sha256=reviewer.sha(tmp_path/'plan.json'),
        native_settings=plan['experimental_native_settings'], original_native_worker=str(native), effective_bound_GiB=40,
        worker_source_sha256={str(p): reviewer.sha(p) for p in (native, wrapper)})
    path = tmp_path/'profile.json'
    path.write_text(json.dumps(profile))
    execution = dict(profile_sha256=reviewer.sha(path), original_worker=str(native),
        original_worker_sha256=reviewer.sha(native), native_numerical_settings_changed=False,
        geometry_or_energy_changed=False, capacity_calls=[dict(original_floor_GiB=112, effective_floor_GiB=72)])
    (tmp_path/'effective_execution.json').write_text(json.dumps(execution))
    guard = tmp_path/'effective_host_guard'
    guard.mkdir()
    (guard/'status.json').write_text(json.dumps(dict(status='effective_memory_guard_finished_profile_review_pending',
        finished=1, maximum_RSS_GiB=40, minimum_free_GiB=32, peak_tree_RSS_bytes=28*2**30, self_VmHWM_bytes=29*2**30)))
    return dict(program=str(wrapper)), plan, tmp_path, path


def test_known_resource_wrapper_requires_full_source_and_peak_chain(runtime):
    assert reviewer.validate_execution(*runtime)


@pytest.mark.parametrize('mutation', ['program', 'source', 'peak', 'settings'])
def test_wrong_program_changed_source_or_failed_profile_cannot_pass(runtime, mutation):
    rt, plan, case, profile = runtime
    if mutation == 'program':rt['program'] = '/unknown/program.py'
    if mutation == 'source':(case/'native.py').write_text('changed equations')
    if mutation == 'peak':
        p = case/'effective_host_guard/status.json'
        data = json.loads(p.read_text());data['self_VmHWM_bytes'] = 41*2**30;p.write_text(json.dumps(data))
    if mutation == 'settings':
        p = case/'effective_execution.json'
        data = json.loads(p.read_text());data['native_numerical_settings_changed'] = True;p.write_text(json.dumps(data))
    with pytest.raises(ValueError):
        reviewer.validate_execution(rt, plan, case, profile)
