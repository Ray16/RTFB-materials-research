"""Prepare a new comparison of two own minima; never select it in production."""
from pathlib import Path
import os
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.validation.spin_ensemble import prepare_policy,evaluate_policy
b=Path(__file__).resolve().parent;sd=ROOT/'calcs/dft/search_pmdi_1864c2c74267/red2';helper=ROOT/'calcs/spincheck/search_pmdi_1864c2c74267/red2_m3'
p=prepare_policy(sd,helper);atomic_write_json(b/'prepared_policy.json',p);r=evaluate_policy(sd,p)
assert not (sd/'spin_ensemble_policy.json').exists()
print({k:v for k,v in r.items() if k not in ['evidence_sha256','chosen_identity']})
