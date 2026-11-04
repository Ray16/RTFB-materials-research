"""Queue retries require exact prelaunch classification, cooldown and a limit."""
import importlib.util
import json
import pytest
from redox.core.hashing import sha256_file

from redox.core.common import ROOT


def module(path,name):
    spec=importlib.util.spec_from_file_location(name,ROOT/path)
    value=importlib.util.module_from_spec(spec);spec.loader.exec_module(value);return value


def test_only_pending_or_cooled_proven_prelaunch_refusals_are_eligible():
    queue=module('scripts/fleet/run_corrected_sampling_queue.py','test_corrected_queue')
    assert queue.ready('pending',0,100,60,100)
    assert not queue.ready('prelaunch_refusal_no_scientific_attempt',1,100,60,159)
    assert queue.ready('prelaunch_refusal_no_scientific_attempt',1,100,60,160)
    assert not queue.ready('prelaunch_refusal_no_scientific_attempt',3,100,60,1000)
    for condition in ('ownership_unverified','running','scientific_failure_preserved','stopped_attempt_preserved'):
        assert not queue.ready(condition,0,0,60,1000000)


def test_finished_status_is_explicit_for_the_corrected_workflow(tmp_path):
    watcher=module('scripts/fleet/run_screening_energy_panel.py','test_corrected_status')
    root=tmp_path/'input';case=root/'case';case.mkdir(parents=True)
    (case/'status.json').write_text(json.dumps(dict(status='corrected_sampling_energy_finished_independent_review_pending',finished=10)))
    row=dict(namespace=str(root/'namespace.json'));alive=lambda owner:False
    assert watcher.inspect(row,root/'dispatch_claim',alive)[0]=='scientific_failure_preserved'
    assert watcher.inspect(row,root/'dispatch_claim',alive,
        finished_status='corrected_sampling_energy_finished_independent_review_pending')[0]=='completed_review_pending'


def test_storage_peer_roster_can_extend_but_never_drop_or_change_history(tmp_path):
    queue=module('scripts/fleet/run_corrected_sampling_queue.py','test_corrected_queue_roster')
    first=tmp_path/'first.json';first.write_text('{}');second=tmp_path/'second.json';second.write_text('{}')
    roster=tmp_path/'peers.json';rows=[dict(namespace=str(first),namespace_sha256=sha256_file(first))]
    roster.write_text(json.dumps(rows));_,old=queue.storage_peers(roster,{})
    rows.append(dict(namespace=str(second),namespace_sha256=sha256_file(second)));roster.write_text(json.dumps(rows))
    assert len(queue.storage_peers(roster,old)[0])==2
    roster.write_text(json.dumps(rows[1:]))
    with pytest.raises(ValueError):queue.storage_peers(roster,old)
    roster.write_text(json.dumps(rows));first.write_text('{"changed":true}')
    with pytest.raises(ValueError):queue.storage_peers(roster,old)
