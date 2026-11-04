"""A reservation race must not exhaust the lead's scientific/prelaunch retry budget."""
import importlib.util
import ast
import hashlib
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
        if error.startswith('RuntimeError: '):
            source=path.parents[2]/'src/redox/workflows/screening_control.py'
            tree=ast.parse(source.read_text());message=error.split(': ',1)[1]
            line=next(n.lineno for n in ast.walk(tree) if isinstance(n,ast.Raise)
                and isinstance(n.exc,ast.Call) and n.exc.args and isinstance(n.exc.args[0],ast.Constant)
                and n.exc.args[0].value==message)
            text=f'Traceback (most recent call last):\n  File "{source}", line {line}, in launch\n'+error+'\n'
        else:text='Traceback\n'+error+'\n'
        (tmp_path/'dispatch_logs'/(name+'.log')).write_text(text)
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


def test_matching_message_without_the_recorded_pre_spawn_branch_is_not_exempt(refusals):
    log=refusals.parent/'dispatch_logs/try_002.log'
    log.write_text('RuntimeError: GPU occupied; no screening control started\n')
    assert module.noncontention_prelaunch_failures(refusals,lambda owner:False)==2


def test_ram_or_storage_refusal_is_resource_wait_not_an_environment_failure(refusals):
    source=path.parents[2]/'src/redox/workflows/screening_control.py';message='Control RAM or output headroom unavailable'
    tree=ast.parse(source.read_text());line=next(n.lineno for n in ast.walk(tree) if isinstance(n,ast.Raise)
        and isinstance(n.exc,ast.Call) and n.exc.args and isinstance(n.exc.args[0],ast.Constant) and n.exc.args[0].value==message)
    (refusals.parent/'dispatch_logs/try_001.log').write_text(f'Traceback (most recent call last):\n  File "{source}", line {line}, in launch\nRuntimeError: {message}\n')
    assert module.noncontention_prelaunch_failures(refusals,lambda owner:False)==0


def test_prior_frozen_launcher_version_proves_its_original_traceback_lines(refusals):
    from redox.workflows.prelaunch_accounting import prelaunch_accounting
    source=path.parents[2]/'src/redox/workflows/screening_control.py'
    content='# A prior frozen launcher version\n'+source.read_text()
    digest=hashlib.sha256(content.encode()).hexdigest();frozen=refusals.parent/'sources'/digest/'screening_control.py'
    frozen.parent.mkdir(parents=True);frozen.write_text(content)
    message='GPU occupied; no screening control started'
    line=next(n.lineno for n in ast.walk(ast.parse(content)) if isinstance(n,ast.Raise)
        and isinstance(n.exc,ast.Call) and n.exc.args and isinstance(n.exc.args[0],ast.Constant) and n.exc.args[0].value==message)
    (refusals.parent/'dispatch_logs/try_002.log').write_text(f'Traceback (most recent call last):\n  File "{source}", line {line}, in launch\nRuntimeError: {message}\n')
    assert prelaunch_accounting(refusals,lambda owner:False)['prelaunch_failures']==2
    result=prelaunch_accounting(refusals,lambda owner:False,[frozen])
    assert result['prelaunch_failures']==1 and result['records'][1]['launcher_source']==str(frozen)
    frozen.write_text(content+'# changed\n')
    with pytest.raises(ValueError):prelaunch_accounting(refusals,lambda owner:False,[frozen])


@pytest.mark.parametrize('message',['GPU occupied; no molecular job started',
    'Host RAM promised; no molecular job started','Shared output headroom unavailable; no molecular job started'])
def test_molecular_delegate_resource_branch_does_not_consume_failure_budget(refusals,message):
    from redox.workflows.prelaunch_accounting import prelaunch_accounting
    source=path.parents[2]/'src/redox/workflows/molecular_native_namespace.py'
    tree=ast.parse(source.read_text());line=next(n.lineno for n in ast.walk(tree) if isinstance(n,ast.Raise)
        and isinstance(n.exc,ast.Call) and n.exc.args and isinstance(n.exc.args[0],ast.Constant)
        and n.exc.args[0].value==message)
    log=refusals.parent/'dispatch_logs/try_001.log'
    log.write_text(f'Traceback (most recent call last):\n  File "{source}", line {line}, in launch\nRuntimeError: {message}\n')
    result=prelaunch_accounting(refusals,lambda owner:False)
    assert result['prelaunch_failures']==0 and result['resource_refusals']==3
    assert result['records'][0]['launcher_source']==str(source)
    # A started wrapper or guard is never reopened by this accounting exception.
    placement=refusals.parent/'placement_1';placement.mkdir();(placement/'launch.json').write_text('{}')
    assert prelaunch_accounting(refusals,lambda owner:False)['prelaunch_failures']==1
