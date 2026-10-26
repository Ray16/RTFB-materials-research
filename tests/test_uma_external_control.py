import csv
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest


source = Path(__file__).resolve().parents[1] / 'scripts/pipeline/uma_external_control.py'
spec = importlib.util.spec_from_file_location('uma_external_control_test', source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def fixture(tmp_path):
    seeds = tmp_path / 'seeds'
    (seeds / 'conformers').mkdir(parents=True)
    xyz = seeds / 'conformers/conf_001.xyz'
    xyz.write_text('1\nfixture\nH 0 0 0\n')
    challenge = tmp_path / 'challenge.json'
    challenge.write_text('{}')
    manifest = seeds / 'manifest.csv'
    fields = ['id','name','family','state','charge','mult_hint','n_e','n_conf','solv_preopt','smiles']
    with manifest.open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields); writer.writeheader()
        for state, charge, mult in [('neu',0,1),('red1',-1,2),('red2',-2,1)]:
            writer.writerow(dict(id=module.GROUP,name='fixture',family='control',state=state,charge=charge,
                                 mult_hint=mult,n_e=charge,n_conf=1,solv_preopt=0,smiles='C'))
    plan = dict(challenge=str(challenge), challenge_sha256=module.digest(challenge),
                proposed_manifest_sha256=module.digest(manifest),source_smiles='C',selected=1,
                entries=[{'path':'conformers/conf_001.xyz','sha256':module.digest(xyz)}],
                quantum_jobs_submitted=0,charged_state_minima_computed=False)
    plan_path = seeds / 'plan.json'; plan_path.write_text(json.dumps(plan))
    review = tmp_path / 'review.json'
    review.write_text(json.dumps(dict(status='external_neutral_seed_geometry_review_passed_not_quantum',
        plan_sha256=module.digest(plan_path), selected=1,
        reviewer_sha256=module.digest(module.ROOT/'scripts/validation/review_external_control_geometries.py'),
        producer_sha256=module.digest(module.ROOT/'scripts/validation/prepare_external_control_geometries.py'))))
    return seeds, review


def test_staging_isolated_and_changed_seed_is_rejected(tmp_path):
    seeds, review = fixture(tmp_path)
    campaign = tmp_path / 'campaign'
    plan = module.prepare(seeds, review, campaign)
    checked, rows = module.validate_campaign(campaign)
    assert checked['group'] == module.GROUP and len(rows) == 3
    assert str(campaign) in plan['output_root']
    staged = campaign / 'library' / module.GROUP / 'conformers/conf_001.xyz'
    staged.write_text('modified')
    with pytest.raises(ValueError, match='changed'):
        module.validate_campaign(campaign)


def test_library_context_restores_state_even_after_failure(tmp_path):
    obj = SimpleNamespace(LIBRARY=Path('/old'), read_manifest=lambda: ['old'])
    original_reader = obj.read_manifest
    with pytest.raises(RuntimeError):
        with module.isolated_library(obj, tmp_path, ['new']):
            assert obj.LIBRARY == tmp_path and obj.read_manifest() == ['new']
            raise RuntimeError('deliberate')
    assert obj.LIBRARY == Path('/old') and obj.read_manifest is original_reader


def test_worker_uses_only_isolated_manifest_and_persists_receipt(tmp_path, monkeypatch):
    from redox.qm import uma, uma_current
    from redox.workflows import resources
    seeds, review = fixture(tmp_path)
    campaign = tmp_path / 'campaign'
    module.prepare(seeds, review, campaign)
    monkeypatch.setattr(resources, 'verify_gate_owner', lambda: None)
    monkeypatch.setattr(module.subprocess, 'check_output', lambda *a, **kw: '7.0')
    monkeypatch.setenv('CUDA_VISIBLE_DEVICES', '0')
    policy = json.loads((module.ROOT/'config/uma_execution.json').read_text())
    monkeypatch.setattr(module.sys, 'prefix', str(Path(policy['legacy_python']).parent.parent))
    original_library, original_reader = uma.LIBRARY, uma.read_manifest

    def fake_run(args, profile, selected_policy):
        assert profile == 'legacy' and args == ['--only', module.GROUP + ':red1', '--device', 'cuda']
        assert uma.LIBRARY == campaign/'library'
        assert [r['state'] for r in uma.read_manifest()] == ['neu', 'red1', 'red2']
        assert selected_policy['legacy_output'] == str(campaign/'output/legacy')
        result = campaign/'output/legacy'/module.GROUP/'red1/result.json'
        result.parent.mkdir(parents=True)
        result.write_text('{"mock": true}')

    monkeypatch.setattr(uma_current, 'run_uma', fake_run)
    module.worker(campaign, 'red1')
    receipt = json.loads((campaign/'completed/red1.json').read_text())
    assert receipt['status'] == 'UMA_control_worker_returned_review_required'
    assert not (campaign/'claims/red1').exists()
    assert uma.LIBRARY == original_library and uma.read_manifest is original_reader


def test_worker_preserves_scientific_failure_and_refuses_retry(tmp_path, monkeypatch):
    from redox.qm import uma_current
    from redox.workflows import resources
    seeds, review = fixture(tmp_path)
    campaign = tmp_path / 'campaign'
    module.prepare(seeds, review, campaign)
    monkeypatch.setattr(resources, 'verify_gate_owner', lambda: None)
    monkeypatch.setattr(module.subprocess, 'check_output', lambda *a, **kw: '7.0')
    monkeypatch.setenv('CUDA_VISIBLE_DEVICES', '0')
    policy = json.loads((module.ROOT/'config/uma_execution.json').read_text())
    monkeypatch.setattr(module.sys, 'prefix', str(Path(policy['legacy_python']).parent.parent))
    def fail(*args):
        raise RuntimeError('synthetic scientific failure')
    monkeypatch.setattr(uma_current, 'run_uma', fail)
    with pytest.raises(RuntimeError, match='synthetic'):
        module.worker(campaign, 'neu')
    assert json.loads((campaign/'failures/neu.json').read_text())['status'] == 'scientific_attempt_failed_preserved'
    assert not (campaign/'claims/neu').exists()
    with pytest.raises(RuntimeError, match='Previous scientific failure'):
        module.worker(campaign, 'neu')


def test_worker_refuses_wrong_interpreter_without_replacing_guard(tmp_path, monkeypatch):
    from redox.workflows import resources
    seeds, review = fixture(tmp_path)
    campaign = tmp_path / 'campaign'
    module.prepare(seeds, review, campaign)
    monkeypatch.setattr(resources, 'verify_gate_owner', lambda: None)
    monkeypatch.setattr(module.subprocess, 'check_output', lambda *a, **kw: '7.0')
    monkeypatch.setattr(module.sys, 'prefix', str(tmp_path / 'wrong-env'))
    monkeypatch.setenv('CUDA_VISIBLE_DEVICES', '0')
    with pytest.raises(RuntimeError, match='before reserving'):
        module.worker(campaign, 'neu')
    assert not (campaign/'claims').exists()
