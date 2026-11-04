#!/usr/bin/env python
"""Materialize the hash-required backing files of the two native certificates."""
import argparse
import hashlib
import json
from pathlib import Path


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def closure(source, checkout, manifest):
    source, checkout = source.resolve(), checkout.resolve()
    projected = json.loads((checkout / 'config/projected_thermal_validation.json').read_text())
    native = json.loads((checkout / 'config/native_thermal_validation.json').read_text())
    required = dict(projected['evidence_sha256'])
    for row in native['controls']:
        name, expected = row['relative_path'], row['source_sha256']
        if name in required and required[name] != expected:
            raise ValueError('Conflicting certificate dependency identity')
        required[name] = expected
    # Validate the entire closure before copying any missing dependency.
    pending = []
    for relative, expected in required.items():
        target = (checkout / relative).resolve()
        original = (source / relative).resolve()
        if (Path(relative).is_absolute() or not target.is_relative_to(checkout)
                or not original.is_relative_to(source)):
            raise ValueError('Certificate dependency escapes its source tree')
        witness = target if target.exists() else original
        if sha(witness) != expected:
            raise ValueError('Certificate dependency differs from its qualified hash: ' + relative)
        if not target.exists():
            pending.append((original, target))
    for original, target in pending:
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(original.read_bytes())
    if any(sha(checkout / p) != h for p, h in required.items()):
        raise ValueError('Certificate closure changed while materializing')
    result = dict(status='certificate_dependency_closure_checked', source=str(source),
        checkout=str(checkout), dependency_sha256=required,
        closure_producer_sha256=sha(Path(__file__)), scientific_certificate_changed=False)
    if manifest.exists():
        if json.loads(manifest.read_text()) != result:
            raise ValueError('Preserve an existing different closure manifest')
    else:
        manifest.parent.mkdir(parents=True, exist_ok=True)
        manifest.write_text(json.dumps(result, indent=2) + '\n')
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('source', 'checkout', 'manifest'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    result = closure(args.source, args.checkout, args.manifest)
    print(json.dumps(dict(status=result['status'], dependencies=len(result['dependency_sha256']))))


if __name__ == '__main__':
    main()
