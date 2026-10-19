"""Read-only hydrogen topology audit of every manifest-addressed saved geometry."""
import collections
import hashlib
import json
import os
import time
from pathlib import Path

from redox.core.common import ROOT, read_manifest
from redox.core.io import atomic_write_json
from redox.validation.protonation import hydrogen_topology_status
from redox.workflows.resources import allocation_guard, cpu_allocation

with allocation_guard():
    os.sched_setaffinity(0, cpu_allocation(1))
records = []
for row in read_manifest():
    for stage, file in [('dft', 'opt.xyz'), ('uma', 'relaxed.xyz')]:
        path = ROOT / 'calcs' / stage / row['id'] / row['state'] / file
        if not path.exists():
            continue
        before = hashlib.sha256(path.read_bytes()).hexdigest()
        result = hydrogen_topology_status(row['smiles'], path)
        assert hashlib.sha256(path.read_bytes()).hexdigest() == before
        records.append(dict(id=row['id'], state=row['state'], stage=stage,
                            path=str(path), sha256=before, hydrogen_topology=result))
report = dict(status='stored_geometry_hydrogen_connectivity_audited',
              records=records, geometry_count=len(records),
              by_stage={s:dict(collections.Counter(r['hydrogen_topology'] for r in records if r['stage']==s)) for s in ('dft','uma')},
              energy_or_thermal_accuracy_validated=False, canonical_outputs_modified=False,
              warnings_require_case_review=True, checked=time.time())
out = ROOT / 'results/validation/hydrogen_topology_2026-10-19'
out.mkdir(exist_ok=True)
atomic_write_json(out / 'current_geometry_audit.json', report)
print(json.dumps({k:v for k,v in report.items() if k!='records'}))
