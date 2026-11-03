"""Resolve declared evidence hash maps, preserving multiple versions per path.

Only named path-to-SHA256 maps are traversed. This is provenance/integrity
inventory, not physical QC or discovery of implicit runtime dependencies.
"""
from collections import defaultdict
import json
from pathlib import Path

from redox.core.hashing import sha256_file

FIELDS = ('source_sha256', 'evidence_sha256', 'producer_sha256')


def bindings(value):
    if isinstance(value, dict):
        for key, item in value.items():
            if key in FIELDS and isinstance(item, dict):
                for name, digest in item.items():
                    if (isinstance(digest, str) and len(digest) == 64
                            and all(c in '0123456789abcdef' for c in digest)):
                        yield name, digest
            else:
                yield from bindings(item)
    elif isinstance(value, list):
        for item in value:
            yield from bindings(item)


def inventory(root, archive, seeds, known_files=()):
    root, archive = Path(root).resolve(), Path(archive).resolve()
    known = defaultdict(set)
    for row in known_files:
        known[row['sha256']].add(row['path'])
    cache = {}; found = {}; missing = {}; relocations = {}; parents = defaultdict(set)
    seed_items = seeds.items() if isinstance(seeds, dict) else seeds
    queue = [(str(Path(p).resolve()), digest, 'seed') for p, digest in seed_items]
    seen = set()

    def actual(path):
        path = Path(path)
        if not path.is_file():
            return None
        stat = path.stat(); signature = (stat.st_size, stat.st_mtime_ns, stat.st_ctime_ns)
        key = str(path)
        if key not in cache or cache[key][0] != signature:
            digest = sha256_file(path)
            end = path.stat()
            if signature != (end.st_size, end.st_mtime_ns, end.st_ctime_ns):
                raise ValueError('Source changed during inventory: ' + key)
            cache[key] = (signature, digest)
        return cache[key][1]

    while queue:
        name, expected, parent = queue.pop(); key = (name, expected)
        parents[key].add(parent)
        if key in seen:
            continue
        seen.add(key)
        original = Path(name); candidates = [original]
        if original.is_relative_to(root):
            candidates.append(archive / original.relative_to(root))
        candidates.extend(Path(p) for p in sorted(known[expected]))
        resolved = next((p for p in candidates if actual(p) == expected), None)
        if resolved is None:
            missing[key] = dict(path=name, expected_sha256=expected, actual_sha256=actual(original))
            continue
        resolved = resolved.resolve(); known[expected].add(str(resolved))
        if str(resolved) != name:
            relocations[key] = dict(original=name, path=str(resolved), sha256=expected)
        physical = (str(resolved), expected)
        if physical in found:
            continue
        path_text = str(resolved)
        category = ('runtime_dependency' if any(s in path_text for s in ('/miniforge3/', '/site-packages/', '/environments/'))
                    else 'third_party_publication' if '/literature/' in path_text and resolved.suffix in ('.pdf','.xml','.html')
                    else 'research_evidence_or_source_snapshot')
        found[physical] = dict(path=path_text, sha256=expected, bytes=resolved.stat().st_size, category=category)
        if resolved.suffix == '.json':
            value = json.loads(resolved.read_text())
            for child, digest in bindings(value):
                child = Path(child)
                if not child.is_absolute():
                    child = root / child
                queue.append((str(child), digest, path_text))
    # A later sibling may have supplied the exact bytes for an earlier miss.
    late = [(name,digest) for name,digest in missing
            if any(actual(p) == digest for p in known[digest])]
    if late:
        retry = inventory(root, archive, late, [dict(path=p,sha256=h) for h, paths in known.items() for p in paths])
        for row in retry['files']:found[(row['path'],row['sha256'])] = row
        for row in retry['explicit_archive_relocations']:
            relocations[(row['original'],row['sha256'])] = row
        retry_missing={(r['path'],r['expected_sha256']) for r in retry['unresolved']}
        for key in list(missing):
            if key in late and key not in retry_missing:missing.pop(key)
        for row in retry['unresolved']:
            key=(row['path'],row['expected_sha256'])
            missing[key]={k:v for k,v in row.items() if k!='referenced_by'}
            parents[key].update(row['referenced_by'])
    for row in found.values():
        if actual(row['path']) != row['sha256']:
            raise ValueError('Evidence changed before inventory completion')
    unique = {r['sha256']:r['bytes'] for r in found.values()}
    return dict(status='declared_dependency_inventory_verified' if not missing else 'declared_dependency_inventory_incomplete',
        files=list(found.values()), explicit_archive_relocations=list(relocations.values()),
        unresolved=[dict(**row, referenced_by=sorted(parents[key])) for key,row in missing.items()],
        unique_content_bytes=sum(unique.values()),
        traversal_fields=list(FIELDS),
        traversal_scope='Recursive named path-to-SHA256 maps, including nested component records; only verified JSON branches traversed',
        implicit_dependencies_covered=False, quantum_validation_repeated=False,
        public_deposition_complete=False)
