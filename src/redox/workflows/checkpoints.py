"""Verify state-completion receipts against the actual saved output bytes."""
import hashlib
import json
from pathlib import Path


def verify_state_checkpoint(receipt_path):
    """Return the state report; refuse missing, mutated or incomplete evidence.

    The receipt is a cache of a verified state, never permission to skip checking
    its outputs. This verifies execution evidence, not experimental accuracy or
    completion of subsequent spin-selection and property-publication reviews.
    """
    receipt = json.loads(Path(receipt_path).read_text())
    report_path = Path(receipt['status_path'])
    content = report_path.read_bytes()
    if hashlib.sha256(content).hexdigest() != receipt['status_sha256']:
        raise ValueError('State report changed after completion')
    report = json.loads(content)
    if report.get('status') != 'state_records_verified' or report.get('task', {}).get('key') != receipt['key']:
        raise ValueError('Completion receipt does not identify a verified state')
    outputs = report.get('output_sha256')
    if not isinstance(outputs, dict) or not outputs:
        raise ValueError('Completion receipt lacks actual output evidence')
    for path, digest in outputs.items():
        if hashlib.sha256(Path(path).read_bytes()).hexdigest() != digest:
            raise ValueError('Completed output changed: ' + path)
    return report
