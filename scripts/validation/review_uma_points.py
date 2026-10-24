#!/usr/bin/env python
"""Replay a frozen UMA energy/force comparison from cached records on the CPU."""
import argparse
import json
from pathlib import Path

from redox.core.io import atomic_write_json
from redox.validation.uma_parity import compare_points
from redox.workflows import polymer_campaign as api


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('plan', 'reference', 'candidate', 'out'):
        parser.add_argument('--'+name, type=Path, required=True)
    args = parser.parse_args()
    paths = {k: getattr(args, k).resolve() for k in ('plan', 'reference', 'candidate')}
    data = {k: json.loads(p.read_text()) for k, p in paths.items()}
    plan_hash = api.digest(paths['plan'])
    if any(data[k]['plan_sha256'] != plan_hash for k in ('reference', 'candidate')):
        raise ValueError('Point records do not belong to the exact frozen panel')
    api.verify_sources(data['plan']['input_sha256'])
    api.verify_sources(data['plan']['producer_sha256'])
    report = compare_points(**data)
    report['source_sha256'] = {str(p): api.digest(p) for p in [*paths.values(), Path(__file__).resolve(),
        Path(__import__('redox.validation.uma_parity', fromlist=['x']).__file__)]}
    atomic_write_json(args.out, report)
    print(json.dumps({k: v for k, v in report.items() if k not in ('comparisons', 'source_sha256')}))


if __name__ == '__main__':
    main()
