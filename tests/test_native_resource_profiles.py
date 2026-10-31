"""Resource admission must preserve the immutable task and avoid head blocking."""
import copy
import json
from pathlib import Path

import pytest

from redox.workflows import native_resource_profiles as profiles


def fits(spec,resources):
    return resources['free_GiB']>=spec['maximum_host_RSS_GiB']+32


def test_original_plan_is_not_mutated_by_the_smaller_scheduling_copy(tmp_path,monkeypatch):
    spec=dict(key='exact',maximum_host_RSS_GiB=80)
    before=copy.deepcopy(spec)
    monkeypatch.setattr(profiles.effective,'verify',lambda _:({'effective_bound_GiB':40},spec))
    assert profiles.admission_spec(spec,None) is spec
    assert profiles.admission_spec(spec,tmp_path/'profile.json')['maximum_host_RSS_GiB']==40
    assert spec==before
    with pytest.raises(ValueError):profiles.admission_spec(dict(spec,key='another'),tmp_path/'profile.json')


def test_a_fitting_small_task_is_not_blocked_by_the_larger_first_task(tmp_path,monkeypatch):
    a,b=tmp_path/'large',tmp_path/'small';large=dict(key='large',maximum_host_RSS_GiB=80);small=dict(key='small',maximum_host_RSS_GiB=80)
    profile=tmp_path/'profile.json'
    monkeypatch.setattr(profiles.effective,'verify',lambda _:({'effective_bound_GiB':40},small))
    pending=[(a,large),(b,small)]
    assert profiles.first_fitting(pending,dict(free_GiB=75),{str(b):profile},fits)==(1,profile)
    assert profiles.first_fitting(pending,dict(free_GiB=120),{str(b):profile},fits)==(0,None)
    assert profiles.first_fitting(pending,dict(free_GiB=70),{str(b):profile},fits)==(None,None)


def test_shared_volume_headroom_controls_new_admission_only(tmp_path,monkeypatch):
    from collections import namedtuple
    usage=namedtuple('usage','total used free')
    monkeypatch.setattr(profiles.shutil,'disk_usage',lambda _:usage(400*2**40,0,63*2**30))
    assert not profiles.output_room(tmp_path)['admit_new_work']
    monkeypatch.setattr(profiles.shutil,'disk_usage',lambda _:usage(400*2**40,0,64*2**30))
    assert profiles.output_room(tmp_path)['admit_new_work']


def fixture(tmp_path,monkeypatch):
    case=tmp_path/'case';case.mkdir();(case/'plan.json').write_text('{}')
    native=tmp_path/'native.json';native.write_text('{}')
    namespace=tmp_path/'namespace.json';namespace.write_text('{}')
    profile=tmp_path/'profile.json';profile.write_text('{}')
    binding=dict(namespace=str(namespace),review=str(native),review_sha256=profiles.sha256_file(native))
    resource=dict(status='effective_native_resource_profile_checked',tested_run_profile_qualified=True,
        experimental_accuracy_inferred=False,normal_native_binding=binding,profile_sha256=profiles.sha256_file(profile),effective_bound_GiB=40)
    path=tmp_path/'resource.json';path.write_text(json.dumps(resource))
    reference=dict(effective_runtime_validation=dict(profile=str(profile)))
    spec=dict(attempt=str(case));settings=dict(native_namespace=str(namespace))
    monkeypatch.setattr(profiles,'read_checked_review',lambda *args:(reference,binding))
    monkeypatch.setattr(profiles.effective,'verify',lambda _: (settings,spec))
    monkeypatch.setattr(profiles.effective,'review',lambda *args:resource)
    return path,resource,native


def test_qualified_reference_is_reconstructed_and_then_source_bound(tmp_path,monkeypatch):
    path,resource,native=fixture(tmp_path,monkeypatch)
    reader=profiles.MatchedNativeProfiles(path)
    reader.verify_sources()
    native.write_text('changed')
    with pytest.raises(ValueError,match='changed'):reader.verify_sources()


@pytest.mark.parametrize('failure',['not_qualified','different_binding','unreconstructed_resource'])
def test_a_passed_label_or_mismatched_reference_cannot_transfer_a_budget(tmp_path,monkeypatch,failure):
    path,resource,native=fixture(tmp_path,monkeypatch)
    if failure=='not_qualified':resource['tested_run_profile_qualified']=False
    if failure=='different_binding':resource['normal_native_binding']['review_sha256']='wrong'
    if failure=='unreconstructed_resource':
        monkeypatch.setattr(profiles.effective,'review',lambda *args:dict(resource,effective_bound_GiB=10))
    path.write_text(json.dumps(resource))
    if failure=='different_binding':
        # The checked reader returns its independently reconstructed digest.
        monkeypatch.setattr(profiles,'read_checked_review',lambda *args:(dict(effective_runtime_validation=dict(profile=str(tmp_path/'profile.json'))),
            dict(resource['normal_native_binding'],review_sha256=profiles.sha256_file(native))))
    with pytest.raises(ValueError):profiles.MatchedNativeProfiles(path)
