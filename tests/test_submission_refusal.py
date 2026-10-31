import json
from pathlib import Path
import subprocess

import pytest

from redox.workflows import submission_refusal as refusal


@pytest.fixture(params=['native', 'energy'])
def declined(tmp_path, monkeypatch, request):
    stage = request.param
    policy = refusal.STAGES[stage]
    monkeypatch.setattr(refusal, 'ROOT', tmp_path/'root')
    source = refusal.ROOT/'src/redox/workflows'/policy['module']
    source.parent.mkdir(parents=True)
    message = policy['messages'][0]
    source.write_text('def launch():\n    if busy:\n        raise RuntimeError('+repr(message)+')\n    child=subprocess.Popen(argv)\n')
    case, raw, archive, output = [tmp_path/name for name in ('case','raw','archive','queue')]
    for path in (case,raw,archive,output):path.mkdir()
    spec = dict(key='exact',attempt=str(case),raw_parent=str(raw),archive=str(archive),prior_root='/old',
        original_campaign='/old/results/polymer_models/campaigns/batch/exact',
        original_raw_parent='/old/calcs/polymer_models/batch',original_raw_root='/old/calcs/polymer_models/batch/exact')
    namespace = tmp_path/'namespace.json';namespace.write_text(json.dumps(spec))
    placement = output/'placements/exact/attempt/placement';placement.mkdir(parents=True)
    (placement.parent/'controller.log').write_text('Traceback (most recent call last):\n'+
        f'  File "{source}", line 3, in launch\n    raise RuntimeError({message!r})\nRuntimeError: '+message+'\n')
    claim = case/policy['claim'];claim.mkdir()
    (claim/'owner.json').write_text(json.dumps(dict(queue=str(output))))
    record=dict(status='controller_spawned_scientific_start_unverified',owner=dict(host='node',pid=1,process_start='start'),
        argv=['python','--out',str(placement)],namespace=str(namespace),namespace_sha256=refusal.sha256_file(namespace))
    (claim/'launch.json').write_text(json.dumps(record if stage=='native' else dict(receipts=[record])))
    monkeypatch.setattr(refusal,'producer_inventory',lambda *args: {'node':[]})
    return spec,output,claim,source,placement,namespace,stage,{str(source):refusal.sha256_file(source)}


def release(fixture, alive=lambda _:False):
    spec,output,claim,source,placement,namespace,stage,proof=fixture
    return refusal.release_prelaunch_refusal(spec,output,alive,['node'],stage,proof)


def test_proven_prelaunch_refusal_preserves_history_without_spending_an_attempt(declined):
    assert release(declined)
    spec,output,claim,*_=declined
    assert not claim.exists()
    history=list((output/'prelaunch_refusal_history/exact').glob('*/resolution.json'))
    assert len(history)==1
    result=json.loads(history[0].read_text())
    assert not result['scientific_attempt_consumed'] and result['scientific_outputs_absent']
    assert (history[0].parent/'claim/launch.json').exists()


@pytest.mark.parametrize('mutation',['live','unknown','changed_source','source_after_spawn','raw','guard','launch',
                                    'wrong_namespace','other_producer','unreachable','unrecognized_log'])
def test_uncertain_or_started_attempts_are_preserved(declined,monkeypatch,mutation):
    spec,output,claim,source,placement,namespace,stage,proof=declined
    alive=lambda _:False
    if mutation=='live':alive=lambda _:True
    if mutation=='unknown':alive=lambda _:None
    if mutation=='changed_source':source.write_text(source.read_text()+'\n')
    if mutation=='source_after_spawn':
        message=refusal.STAGES[stage]['messages'][0]
        source.write_text('def launch():\n    child=subprocess.Popen(argv)\n    raise RuntimeError('+repr(message)+')\n')
        proof[str(source)]=refusal.sha256_file(source)
    if mutation=='raw':(Path(spec['raw_parent'])/'exact').mkdir()
    if mutation=='guard':(Path(spec['attempt'])/'residency_guard').mkdir()
    if mutation=='launch':(placement/'launch.json').write_text('{}')
    if mutation=='wrong_namespace':namespace.write_text('{}')
    if mutation=='other_producer':monkeypatch.setattr(refusal,'producer_inventory',lambda *args:{'node':[dict(pid=2)]})
    if mutation=='unreachable':
        def timeout(*args):raise subprocess.TimeoutExpired('ssh',25)
        monkeypatch.setattr(refusal,'producer_inventory',timeout)
    if mutation=='unrecognized_log':(placement.parent/'controller.log').write_text('a different failure')
    assert not release(declined,alive)
    assert claim.exists()


def test_second_owner_observation_timeout_cannot_reopen_the_claim(declined):
    observations=iter([False,subprocess.TimeoutExpired('ssh',20)])
    def alive(owner):
        value=next(observations)
        if isinstance(value,Exception):raise value
        return value
    assert not release(declined,alive)
    assert declined[2].exists()
