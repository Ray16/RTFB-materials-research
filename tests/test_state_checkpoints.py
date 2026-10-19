import hashlib
import json

import pytest

from redox.workflows.checkpoints import verify_state_checkpoint


def checkpoint(tmp_path, **updates):
    output = tmp_path / 'opt.xyz'
    output.write_text('2\n\nH 0 0 0\nH 0 0 .74\n')
    report = dict(status='state_records_verified', task=dict(key='m__neu'),
                  output_sha256={str(output): hashlib.sha256(output.read_bytes()).hexdigest()})
    report.update(updates)
    status = tmp_path / 'status.json'
    status.write_text(json.dumps(report))
    receipt = tmp_path / 'done.json'
    receipt.write_text(json.dumps(dict(key='m__neu', status_path=str(status),
        status_sha256=hashlib.sha256(status.read_bytes()).hexdigest())))
    return receipt, status, output


def test_completion_requires_unchanged_actual_outputs(tmp_path):
    receipt, _, output = checkpoint(tmp_path)
    assert verify_state_checkpoint(receipt)['task']['key'] == 'm__neu'
    output.write_text('a later geometry replacement')
    with pytest.raises(ValueError, match='output changed'):
        verify_state_checkpoint(receipt)


def test_missing_output_and_changed_state_report_never_pass(tmp_path):
    receipt, status, output = checkpoint(tmp_path)
    output.unlink()
    with pytest.raises(FileNotFoundError):
        verify_state_checkpoint(receipt)
    status.write_text('{}')
    with pytest.raises(ValueError, match='report changed'):
        verify_state_checkpoint(receipt)


@pytest.mark.parametrize('updates', [dict(status='computing'), dict(task=dict(key='other__neu')),
                                    dict(output_sha256={})])
def test_marker_alone_or_wrong_state_is_insufficient(tmp_path, updates):
    receipt, _, _ = checkpoint(tmp_path, **updates)
    with pytest.raises(ValueError):
        verify_state_checkpoint(receipt)
