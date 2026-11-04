#!/usr/bin/env python3
"""Publish verified canonical locations in the mutable results index."""
import argparse
import json
from pathlib import Path

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.workspace_index import relocate_index


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--index', type=Path, default=ROOT/'results/INDEX.json')
    parser.add_argument('--mapping', type=Path, default=ROOT/'validation/workspace_relocation.json')
    parser.add_argument('--report', type=Path, required=True)
    args = parser.parse_args()
    if args.report.exists():
        raise FileExistsError('Preserve earlier relocation reports')
    previous = sha256_file(args.index)
    result, changes = relocate_index(json.loads(args.index.read_text()),
                                     json.loads(args.mapping.read_text()), ROOT)
    if sha256_file(args.index) != previous:
        raise ValueError('Results index changed during preparation')
    atomic_write_json(args.index, result)
    report = dict(status='mutable_result_pointers_relocated',changes=changes,
                  previous_index_sha256=previous,current_index_sha256=sha256_file(args.index),
                  mapping_sha256=sha256_file(args.mapping),evidence_files_rewritten=False,
                  scientific_status_changed=False,compatibility_links_removed=False,
                  deferred_roots=json.loads(args.mapping.read_text())['deferred_roots'])
    atomic_write_json(args.report, report)
    print(json.dumps(dict(status=report['status'],changed_pointers=len(changes))))


if __name__ == '__main__':
    main()
