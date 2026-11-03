"""A reservation race must not exhaust the lead's scientific/prelaunch retry budget."""
import importlib.util
import json
from pathlib import Path

import pytest

path=Path(__file__).resolve().parents[1]/'scripts/fleet/run_saved_ionic_validation.py'
spec=importlib.util.spec_from_file_location('saved_ionic_retry_accounting',path)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


@pytest.fixture
def refusals(tmp_path):
    (tmp_path/'dispatch_claim').mkdir();(tmp_path/'dispatch_logs').mkdir()
    for index,error in enumerate(['ValueError: storage profile','RuntimeError: GPU occupied; no screening control started',
                                  'RuntimeError: GPU occupied; no screening control started'],1):
        name=f'try_{index:03d}'
        (tmp_path/'dispatch_claim'/(name+'.json')).write_text(json.dumps(dict(owner={'pid':index},placement=str(tmp_path/('placement_'+str(index))))))
        (tmp_path/'dispatch_logs'/(name+'.log')).write_text('Traceback\n'+error+'\n')
    return tmp_path/'namespace.json'


def test_two_closed_gpu_races_leave_the_original_single_environment_failure(refusals):
    assert module.noncontention_prelaunch_failures(refusals,lambda owner:False)==1


@pytest.mark.parametrize('alive',[True,None])
def test_live_or_uncertain_controller_cannot_be_retried(refusals,alive):
    with pytest.raises(ValueError):module.noncontention_prelaunch_failures(refusals,lambda owner:alive)


def test_started_worker_is_preserved_even_with_an_old_contention_message(refusals):
    p=refusals.parent/'case';p.mkdir();(p/'status.json').write_text('{}')
    with pytest.raises(ValueError):module.noncontention_prelaunch_failures(refusals,lambda owner:False)


def test_a_launched_wrapper_is_not_classified_as_the_prelaunch_gpu_refusal(refusals):
    p=refusals.parent/'placement_2';p.mkdir();(p/'launch.json').write_text('{}')
    assert module.noncontention_prelaunch_failures(refusals,lambda owner:False)==2
