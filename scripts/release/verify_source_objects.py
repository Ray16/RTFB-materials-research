"""Verify a relocated raw-evidence package with the Python standard library."""
import hashlib
import json
from pathlib import Path
import sys


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(2**20), b''):
            h.update(block)
    return h.hexdigest()


def verify(folder):
    folder = Path(folder).resolve()
    data = json.loads((folder / 'manifest.json').read_text())
    expected = set()
    for digest, size in data['objects'].items():
        if len(digest) != 64 or any(c not in '0123456789abcdef' for c in digest):
            raise ValueError('Invalid object identity')
        path = folder / 'objects' / digest
        if path.is_symlink() or path.stat().st_size != size or sha256(path) != digest:
            raise ValueError('Object differs: ' + digest)
        expected.add(path.name)
    if {p.name for p in (folder / 'objects').iterdir()} != expected:
        raise ValueError('Unexpected or missing object')
    for row in data['files']:
        if row['sha256'] not in data['objects']:
            raise ValueError('Missing mapped evidence object')
    for row in data['relocations']:
        if row['sha256'] not in data['objects']:
            raise ValueError('Missing relocated evidence object')
    return dict(status='raw_evidence_objects_and_mappings_verified',
                verified_objects=len(expected), file_mappings=len(data['files']),
                bytes=sum(data['objects'].values()),
                external_dependencies_verified=False,
                quantum_validation_repeated=False, public_deposition_complete=False)


if __name__ == '__main__':
    print(json.dumps(verify(Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent)))
