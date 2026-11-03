"""Resolve absent logical paths using already inventoried, identical bytes."""
from collections import defaultdict
import hashlib
import json
from pathlib import Path


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


folder = Path(__file__).resolve().parent
source = folder / 'inventory.json'
output = folder / 'content_resolution.json'
if output.exists():
    raise FileExistsError('Preserve earlier content-resolution records')
inventory = json.loads(source.read_text())
known = defaultdict(list)
for row in inventory['files']:
    known[row['sha256']].append(row['path'])
resolved, remaining = [], []
for row in inventory['unresolved']:
    candidates = sorted(known[row['expected_sha256']])
    if not candidates:
        remaining.append(row)
        continue
    actual = candidates[0]
    if digest(actual) != row['expected_sha256']:
        raise ValueError('Inventoried replacement bytes changed: ' + actual)
    resolved.append(dict(original=row['path'], path=actual,
                         sha256=row['expected_sha256'],
                         referenced_by=row['referenced_by']))
report = dict(status='known_content_relocated_remaining_dependencies_unresolved',
              inventory=str(source), inventory_sha256=digest(source),
              resolved=resolved, unresolved=remaining,
              resolution_rule='Exact recorded SHA256 matches an already traversed inventoried file',
              input_records_modified=False, quantum_validation_repeated=False,
              complete_dependency_closure=False, public_deposition_complete=False,
              producer_sha256=digest(__file__))
output.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(dict(resolved=len(resolved), unresolved=len(remaining))))
