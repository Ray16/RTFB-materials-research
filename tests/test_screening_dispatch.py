"""Do not retry scientific or uncertain starts, including PID-reuse boundaries."""
import importlib.util
import json
from pathlib import Path

import pytest

from redox.core.common import ROOT


@pytest.fixture
def inspect_case(tmp_path):
    spec=importlib.util.spec_from_file_location('tested_screening_dispatch',ROOT/'scripts/fleet/run_screening_energy_panel.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    (tmp_path/'case').mkdir()
    return module.inspect,dict(namespace=str(tmp_path/'namespace.json')),tmp_path/'dispatch_claim',tmp_path


def test_unowned_pending_case_and_uncertain_claim_are_distinct(inspect_case):
    inspect,row,claim,_=inspect_case
    assert inspect(row,claim,lambda owner:False)[0]=='pending'
    claim.mkdir()
    assert inspect(row,claim,lambda owner:False)[0]=='uncertain_submission_preserved'


def test_live_scientific_worker_is_not_reclaimed(inspect_case):
    inspect,row,claim,root=inspect_case
    (root/'case/status.json').write_text(json.dumps(dict(status='screening_energy_control_running',owner=dict(pid=123))))
    assert inspect(row,claim,lambda owner:True)[0]=='running'
    assert inspect(row,claim,lambda owner:None)[0]=='ownership_unverified'
    assert inspect(row,claim,lambda owner:False)[0]=='stopped_attempt_preserved'


def test_dead_controller_with_live_wrapper_cannot_be_retried(inspect_case):
    inspect,row,claim,root=inspect_case
    claim.mkdir();placement=root/'placement';placement.mkdir()
    (claim/'try_001.json').write_text(json.dumps(dict(owner=dict(pid=12),placement=str(placement))))
    (placement/'launch.json').write_text(json.dumps(dict(host='node',pid=34,process_start='expected')))
    assert inspect(row,claim,lambda owner:owner['pid']==34)[0]=='submitted'
    assert inspect(row,claim,lambda owner:False)[0]=='prelaunch_refusal_no_scientific_attempt'


def test_scientific_failure_remains_preserved_even_after_owner_exit(inspect_case):
    inspect,row,claim,root=inspect_case
    (root/'case/status.json').write_text(json.dumps(dict(status='screening_energy_control_failed_preserved',finished=1,owner=dict(pid=12))))
    assert inspect(row,claim,lambda owner:False)[0]=='scientific_failure_preserved'
