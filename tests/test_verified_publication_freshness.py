"""A published value must be withdrawn when its scientific evidence changes."""
import hashlib
import importlib.util
import json
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_shared_inventory_includes_selected_equilibrium_and_its_bound_artifacts(tmp_path,monkeypatch):
    from redox.validation.provenance import spin_decision_paths
    spec=importlib.util.spec_from_file_location('ensemble_inventory_overlay',ROOT/'scripts/plotting/candidates/verified_overlay.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    monkeypatch.setattr(module,'ROOT',tmp_path)
    component=tmp_path/'calcs/spincheck/mol/red2_m3/native_thermal/hessian.npy'
    component.parent.mkdir(parents=True);component.write_bytes(b'qualified array')
    policy=tmp_path/'calcs/dft/mol/red2/spin_ensemble_policy.json';policy.parent.mkdir(parents=True)
    policy.write_text(json.dumps({'evidence_sha256':{str(component):hashlib.sha256(component.read_bytes()).hexdigest()}}))
    paths=spin_decision_paths(tmp_path,['mol'])
    assert {policy,component}<=paths
    summary=dict(spin_input_gids=['mol'],spin_input_paths=sorted(str(p) for p in paths),source_sha256={})
    assert not module.input_inventory_changed(summary)
    policy.unlink()
    assert module.input_inventory_changed(summary)


@pytest.fixture
def feed(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location('publication_freshness', ROOT / 'scripts/plotting/candidates/verified_overlay.py')
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    monkeypatch.setattr(module, 'ROOT', tmp_path)
    source = tmp_path / 'config/spin_overrides.json'; source.parent.mkdir(); source.write_text('{}')
    evidence = tmp_path / 'calcs/uma/mol/neu/result.json'; evidence.parent.mkdir(parents=True)
    evidence.write_text('{"spin_gap_eV": 1.0}')
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    summary = dict(arithmetic_discrepancies=0, inputs_changed_during_snapshot=[],
                   source_sha256={str(source): digest(source)},
                   admitted_input_sha256={str(evidence): digest(evidence)},
                   spin_input_gids=['mol'], spin_input_paths=[str(evidence)])
    path = tmp_path / 'feed.json'
    payload = dict(summary=summary, original39=[dict(id='mol',registry_id='mol',solution_potentials_status='verified native calculation')],
                   extra_NDI=[], snapshot_path=str(tmp_path))
    path.write_text(json.dumps(payload))
    return module, path, source, evidence, summary


def test_unchanged_evidence_can_publish(feed):
    module, path, _, _, summary = feed
    assert module.snapshot_inputs_current(summary)
    _, rows, admitted = module.load(path)
    assert len(rows) == len(admitted) == 1


@pytest.mark.parametrize('which', ['source', 'spin'])
def test_changed_spin_or_model_withdraws_old_numeric_points(feed, which):
    module, path, source, evidence, summary = feed
    (source if which == 'source' else evidence).write_text('{"new_ground_state": true}')
    assert not module.snapshot_inputs_current(summary)
    _, rows, admitted = module.load(path)
    assert not admitted
    if which=='source':assert not rows
    # A changed molecule may remain in the raw snapshot; it cannot be a plotted point.


def test_new_spin_report_invalidates_old_snapshot_even_without_file_mutation(feed):
    module, path, _, _, summary = feed
    report = module.ROOT / 'calcs/dft/mol/neu/spin_validation/new.json'
    report.parent.mkdir(parents=True); report.write_text('{"status": "lower_spin_promoted"}')
    assert not module.snapshot_inputs_current(summary)
    assert module.load(path)[1:] == ([], [])


def test_additional_cohort_can_contain_non_ndi_candidates(feed):
    module, path, _, _, _ = feed
    payload = json.loads(path.read_text())
    payload.pop('extra_NDI')
    payload['feed_schema'] = 2
    payload['additional_candidates'] = [dict(id='new_phenothiazine', registry_id='new_phenothiazine',
                                            solution_potentials_status='pending native calculation')]
    path.write_text(json.dumps(payload))
    _, rows, admitted = module.load(path)
    assert [row['id'] for row in rows] == ['mol', 'new_phenothiazine']
    assert [row['id'] for row in admitted] == ['mol']


def test_external_control_dependencies_are_discovered(tmp_path):
    from redox.validation.provenance import referenced_files
    proof = tmp_path / 'proof/status.json'; proof.parent.mkdir(); proof.write_text('{}')
    matrix = tmp_path / 'proof/matrix.npy'; matrix.write_bytes(b'matrix')
    metadata = dict(controls=[dict(relative_path='proof/status.json')],
                    execution=dict(target_evidence_sha256={'proof/matrix.npy':'a'*64}))
    assert referenced_files(metadata,tmp_path) == {proof,matrix}


def test_new_inputs_for_a_pending_other_molecule_do_not_hide_finished_points(feed):
    module,path,_,evidence,summary=feed
    payload=json.loads(path.read_text())
    summary['spin_input_gids'].append('pending')
    summary['spin_input_paths_by_gid']={'mol':[str(evidence)],'pending':[]}
    payload['summary']=summary
    payload['additional_candidates']=[dict(id='pending',registry_id='pending',solution_potentials_status='pending')]
    path.write_text(json.dumps(payload))
    new=module.ROOT/'calcs/spincheck/pending/red2_m3/result.json'
    new.parent.mkdir(parents=True);new.write_text('{}')
    assert not module.snapshot_inputs_current(summary)  # full inventory did change
    assert [r['id'] for r in module.load(path)[2]]==['mol']


def test_collector_can_explicitly_hold_one_changed_molecule_and_publish_others(feed):
    module,path,_,evidence,summary=feed
    payload=json.loads(path.read_text());changed=str(module.ROOT/'calcs/spincheck/held/red2_m3/new.json')
    summary.update(inputs_changed_during_snapshot=[changed],global_input_changed=False,
        input_change_scope={changed:['held']},quarantined_input_gids=['held'],
        spin_input_paths_by_gid={'mol':[str(evidence)]})
    payload['summary']=summary
    # Defend against a stale row status too; the quarantine list takes precedence.
    payload['additional_candidates']=[dict(id='held',registry_id='held',solution_potentials_status='verified native calculation')]
    path.write_text(json.dumps(payload))
    assert module.snapshot_inputs_current(summary)
    assert [r['id'] for r in module.load(path)[2]]==['mol']
    summary['input_change_scope'][changed]=['mol']
    path.write_text(json.dumps(payload))
    assert not module.snapshot_inputs_current(summary)
    assert module.load(path)[2]==[]
