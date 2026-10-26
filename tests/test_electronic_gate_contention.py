import importlib.util
import json
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('capacity_v3_test',ROOT/'scripts/pipeline/coordinate_with_capacity_preflight_v3.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def fixture(tmp_path):
    case=tmp_path/'case';attempt=case/'dispatch_attempts/123';attempt.mkdir(parents=True)
    (case/'plan.json').write_text(json.dumps({'raw_root':str(tmp_path/'raw')}))
    handle={'status':'gate_started_scientific_start_pending','GPU':4,'owner':{'host':'node','pid':99,'process_start':'1'},'attempt':str(attempt)}
    (case/'launch.json').write_text(json.dumps(handle))
    (case/'dispatch_intent.json').write_text(json.dumps({'host':'node','GPU':4}))
    (attempt/'worker.log').write_text('gpu_reserve: gpu4 BUSY/reserved -> not running\n')
    return case,handle


def test_explicit_gate_refusal_is_archived_and_reopened(tmp_path):
    case,handle=fixture(tmp_path)
    audit=module.recover_gate_contention(case,handle,owner_confirmed_dead=True)
    assert json.loads(Path(audit).read_text())['scientific_attempts']==0
    assert not (case/'launch.json').exists()
    assert json.loads((case/'dispatch_intent.json').read_text())['outcome']=='not_started'
    assert (Path(handle['attempt'])/'worker.log').is_file()


@pytest.mark.parametrize('mutation',['live_owner','other_error','raw_output','guard_exists','changed_launch'])
def test_ambiguous_or_started_cases_stay_held(tmp_path,mutation):
    case,handle=fixture(tmp_path)
    if mutation=='other_error':(Path(handle['attempt'])/'worker.log').write_text('worker failed after start')
    if mutation=='raw_output':(tmp_path/'raw').mkdir()
    if mutation=='guard_exists':(case/'residency_guard').mkdir()
    if mutation=='changed_launch':(case/'launch.json').write_text('{}')
    assert module.recover_gate_contention(case,handle,owner_confirmed_dead=mutation!='live_owner') is None
    assert (case/'launch.json').is_file() and not (case/'gate_contention_history').exists()


def test_controller_template_compiles_with_gate_recovery():
    upstream=module.api.load_module('gate_recovery_controller_test',ROOT/'scripts/pipeline/coordinate_smooth_electronic_batch.py')
    patched=module.coordinate_with_durable_errors(upstream)
    assert callable(patched) and '_recover_gate_contention' in patched.__code__.co_names


def test_transition_fails_before_signaling_unexpected_owners(tmp_path,monkeypatch):
    path=ROOT/'scripts/pipeline/transition_electronic_gate_recovery.py'
    spec=importlib.util.spec_from_file_location('transition_gate_test',path)
    transition=importlib.util.module_from_spec(spec);spec.loader.exec_module(transition)
    launch=tmp_path/'launch.json';launch.write_text(json.dumps({'jobs':[]}))
    monkeypatch.setattr(transition.os,'killpg',lambda *a:pytest.fail('must not signal'))
    with pytest.raises(ValueError,match='Expected three'):
        transition.transition(launch,tmp_path/'audit')
    assert not (tmp_path/'audit').exists()
