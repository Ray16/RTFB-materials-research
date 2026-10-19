"""A fixed native record cannot conceal changed external qualification evidence."""
import hashlib
import json

from redox.workflows import tables


def test_inventory_detects_changed_external_record_evidence(tmp_path, monkeypatch):
    monkeypatch.setattr(tables, 'ROOT', tmp_path)
    monkeypatch.setattr(tables, 'CALCS', tmp_path / 'calcs')
    directory = tmp_path / 'calcs/dft/case/neu/native_thermal'
    directory.mkdir(parents=True)
    evidence = tmp_path / 'results/private/gradient.npy'
    evidence.parent.mkdir(parents=True)
    evidence.write_bytes(b'original saved gradient')
    record = directory / 'address.json'
    record.write_text(json.dumps({'execution_derivatives': {
        'target_evidence_sha256': {str(evidence): hashlib.sha256(evidence.read_bytes()).hexdigest()}
    }}))
    before = tables.snapshot_inputs(tmp_path / 'results')
    assert str(evidence) in before
    evidence.write_bytes(b'changed saved gradient')
    after = tables.snapshot_inputs(tmp_path / 'results')
    assert before[str(record)] == after[str(record)]
    assert before[str(evidence)] != after[str(evidence)]
