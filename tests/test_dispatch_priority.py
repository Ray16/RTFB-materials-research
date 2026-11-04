"""Only a live pending lead can defer one compatible free slot from bulk work."""
import json
from pathlib import Path

import pytest

from redox.core.hashing import sha256_file
from redox.workflows.charged_native_followup import priority_requests,priority_for_slot


@pytest.fixture
def priority_case(tmp_path):
    case=tmp_path/'case';case.mkdir();plan=case/'plan.json'
    plan.write_text(json.dumps(dict(minimum_GPU_memory_MiB=80*1024,minimum_GPU_compute_capability=9,
                                   maximum_host_RSS_GiB=80,maximum_durable_output_bytes=6*2**30)))
    namespace=tmp_path/'namespace.json';namespace.write_text(json.dumps(dict(case=str(case),plan_sha256=sha256_file(plan))))
    status=tmp_path/'status.json';status.write_text(json.dumps(dict(namespace=str(namespace),owner={'pid':1})))
    row=dict(namespace=str(namespace),namespace_sha256=sha256_file(namespace),controller_status=str(status),deferred_slots=1,reason='lead')
    config=tmp_path/'priority.json';config.write_text(json.dumps(dict(schema=1,requests=[row])))
    return config,row,case,status


def test_one_live_lead_defers_only_one_suitable_free_slot_and_never_itself(priority_case):
    config,row,case,status=priority_case;requests=priority_requests(config)
    observed=dict(GPU_MiB=95830,GPU_compute_capability=9.)
    selected=priority_for_slot(requests,{'namespace':'bulk'},observed,set(),lambda o:True,lambda p,o:True)
    assert selected[0]==row
    assert priority_for_slot(requests,{'namespace':'bulk'},observed,{row['namespace']},lambda o:True,lambda p,o:True) is None
    assert priority_for_slot(requests,{'namespace':row['namespace']},observed,set(),lambda o:True,lambda p,o:True) is None


@pytest.mark.parametrize('change',['small_gpu','wrong_architecture','insufficient_RAM','dead','unknown','started','guard','finished'])
def test_no_extra_slot_is_deferred_for_unsuitable_or_inactive_jobs(priority_case,change):
    config,row,case,status=priority_case;requests=priority_requests(config);observed=dict(GPU_MiB=95830,GPU_compute_capability=9.)
    live=True;fits=True
    if change=='small_gpu':observed['GPU_MiB']=32768
    if change=='wrong_architecture':observed['GPU_compute_capability']=8.
    if change=='insufficient_RAM':fits=False
    if change=='dead':live=False
    if change=='unknown':live=None
    if change=='started':(case/'status.json').write_text('{}')
    if change=='guard':(case/'residency_guard').mkdir();(case/'residency_guard/owner.json').write_text('{}')
    if change=='finished':status.write_text(json.dumps(dict(namespace=row['namespace'],owner={'pid':1},finished=1)))
    assert priority_for_slot(requests,{'namespace':'bulk'},observed,set(),lambda o:live,lambda p,o:fits) is None


def test_changed_priority_plan_or_namespace_cannot_silently_hold_capacity(priority_case):
    config,row,case,status=priority_case
    (case/'plan.json').write_text('{}')
    with pytest.raises(ValueError):priority_for_slot([row],{'namespace':'bulk'},dict(GPU_MiB=95830,GPU_compute_capability=9.),set(),lambda o:True,lambda p,o:True)
    Path(row['namespace']).write_text('{}')
    with pytest.raises(ValueError):priority_requests(config)


def test_one_shared_controller_prioritizes_only_its_current_pending_input(priority_case,tmp_path):
    config,first,case,status=priority_case
    other=tmp_path/'other';other.mkdir();(other/'plan.json').write_bytes((case/'plan.json').read_bytes())
    namespace=other/'namespace.json';namespace.write_text(json.dumps(dict(case=str(other),plan_sha256=sha256_file(other/'plan.json'))))
    second=dict(first,namespace=str(namespace),namespace_sha256=sha256_file(namespace))
    requests=[first,second];observed=dict(GPU_MiB=95830,GPU_compute_capability=9.)
    # The two immutable requests share one controller's next-input field.
    assert priority_for_slot(requests,{'namespace':first['namespace']},observed,set(),lambda o:True,lambda p,o:True) is None
    assert priority_for_slot(requests,{'namespace':'bulk'},observed,set(),lambda o:True,lambda p,o:True)[0]==first
    (case/'status.json').write_text('{}')
    status.write_text(json.dumps(dict(namespace=second['namespace'],owner={'pid':1})))
    assert priority_for_slot(requests,{'namespace':'bulk'},observed,set(),lambda o:True,lambda p,o:True)[0]==second
    assert priority_for_slot(requests,{'namespace':second['namespace']},observed,set(),lambda o:True,lambda p,o:True) is None
