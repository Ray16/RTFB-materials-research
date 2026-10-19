"""Build the nine previously curated graphs in an isolated library, before QM."""
import csv
import hashlib
import json
import os
import pprint
import runpy
import sys
import time
from pathlib import Path

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity

base = Path(__file__).resolve().parent
source = base.parent / 'aryl_viologen_2026-10-18'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
with allocation_guard():
    resources = host_capacity(1)
    os.sched_setaffinity(0, cpu_allocation(1))
assert not (base / 'isolated_library').exists(), 'Preserve prior preparation; inspect before reuse'
groups = json.loads((source / 'proposals.json').read_text())
pool = json.loads((source / 'structure_pool.json').read_text())
assert len(groups) == len(pool) == 9
assert sha(ROOT / 'data/curated/discovery/aryl_viologen_2025.json') == json.loads((source / 'summary.json').read_text())['source_curated_sha256']
config = base / 'prepared_config.py'
config.write_text('SCAFFOLD = "Cc1ccc(C[*:1])cc1"\nGROUPS = ' + pprint.pformat(groups, width=100, sort_dicts=False) + '\n')
sys.argv = ['redox.build.candidates', '--config', str(config), '--library', str(base / 'isolated_library')]
runpy.run_module('redox.build.candidates', run_name='__main__')
from redox.qm.uma import geometry_checks
from redox.validation.protonation import hydrogen_topology_status
with (base / 'isolated_library/manifest.csv').open() as stream:
    rows = list(csv.DictReader(stream))
assert len(rows) == 27
checks = []
for row in rows:
    for seed in sorted((base / 'isolated_library' / row['id'] / 'conformers').glob('conf_*.xyz')):
        assert not geometry_checks(row['smiles'], seed, int(row['charge']), int(row['mult_hint']))
        assert hydrogen_topology_status(row['smiles'], seed) == 'ok'
        checks.append(dict(id=row['id'], state=row['state'], seed=str(seed), sha256=sha(seed), atoms=int(seed.read_text().splitlines()[0])))
atomic_write_json(base / 'preparation.json', dict(status='nine_isolated_aryl_graphs_and_27_state_definitions_verified',
    candidates=9, states=27, seed_checks=len(checks), maximum_atoms=max(c['atoms'] for c in checks),
    checks=checks, source_sha256={str(p):sha(p) for p in [source / 'proposals.json', source / 'structure_pool.json', ROOT / 'data/curated/discovery/aryl_viologen_2025.json', config]},
    host_preflight=resources, registered_candidates=0, solution_properties_computed=0,
    protonation_model='Frozen isolated acid/phenol stoichiometry; not established MeCN equilibrium speciation',
    source_measurements_inherited=False, created=time.time()))
print(f'Verified {len(checks)} seed/state checks; shared registry unchanged')
