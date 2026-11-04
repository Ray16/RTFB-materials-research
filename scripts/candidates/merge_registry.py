#!/usr/bin/env python
"""Publish a frozen union of the enumerated pool, evaluations and proposal panels."""
import argparse
import json
from redox.screening.candidate_registry import publish
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--spec', default='config/candidate_registry_sources.json')
p.add_argument('--archive', required=True)
p.add_argument('--out', required=True)
p.add_argument('--route-evidence', help='Explicit updated exact-route catalog; original archive remains immutable')
a=p.parse_args()
import os
from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
print(json.dumps(publish(a.spec,a.archive,a.out,route_evidence=a.route_evidence),indent=2))
