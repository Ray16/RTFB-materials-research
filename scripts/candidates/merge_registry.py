#!/usr/bin/env python
"""Publish a frozen union of the enumerated pool, evaluations and proposal panels."""
import argparse
import json
from redox.screening.candidate_registry import publish
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--spec', default='config/candidate_registry_sources.json')
p.add_argument('--archive', required=True)
p.add_argument('--out', required=True)
a=p.parse_args()
print(json.dumps(publish(a.spec,a.archive,a.out),indent=2))
