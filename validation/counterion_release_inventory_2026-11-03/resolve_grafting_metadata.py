"""Resolve the ten historical grafting metadata bindings without changing bytes."""
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
FOLDER = Path(__file__).resolve().parent
SOURCE = FOLDER / 'content_resolution.json'
OUT = FOLDER / 'grafting_metadata_resolution.json'
if OUT.exists():
    raise FileExistsError('Preserve earlier resolution reports')
digest = lambda path: hashlib.sha256(Path(path).read_bytes()).hexdigest()
previous = json.loads(SOURCE.read_text())
resolved = []
for row in previous['unresolved']:
    relative = Path(row['path']).relative_to(ROOT / 'results/grafting/campaigns')
    campaign, workspace, calcs, dft, *tail = relative.parts
    if (workspace, calcs, dft) != ('workspace', 'calcs', 'dft'):
        raise ValueError('Unexpected historical grafting path')
    target = ROOT.parent / 'archive/0_redox/calcs/grafting' / campaign
    target = target.joinpath(*tail)
    if digest(target) != row['expected_sha256']:
        raise ValueError('Relocated metadata bytes differ: ' + str(target))
    value = json.loads(target.read_text())
    for field in ('source_sha256', 'evidence_sha256', 'producer_sha256'):
        if isinstance(value.get(field), dict) and value[field]:
            raise ValueError('New dependency branch requires further inventory')
    resolved.append(dict(original=row['path'], path=str(target),
                         sha256=row['expected_sha256'], bytes=target.stat().st_size,
                         referenced_by=row['referenced_by']))
if len(resolved) != 10:
    raise ValueError('Expected exactly the ten recorded metadata bindings')
report = dict(status='all_declared_inventory_bindings_located', resolved=resolved,
              previous_resolution=str(SOURCE), previous_resolution_sha256=digest(SOURCE),
              unresolved=[], new_declared_hash_map_dependencies=0,
              implicit_dependencies_covered=False, public_deposition_complete=False,
              input_records_modified=False, quantum_validation_repeated=False,
              producer_sha256=digest(__file__))
OUT.write_text(json.dumps(report, indent=2) + '\n')
print(json.dumps(dict(resolved=len(resolved), unresolved=0)))
