"""Freeze isolated UMA tasks and the three-charge 70-atom resource pilot."""
import csv
import json
import os
import time
from pathlib import Path

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.qm import uma

base = Path(__file__).resolve().parent
assert not (base / 'uma_plan.json').exists(), 'Existing declaration must stay immutable'
prepared = json.loads((base / 'preparation.json').read_text())
assert prepared['states'] == 27 and prepared['maximum_atoms'] == 70
with (base / 'isolated_library/manifest.csv').open() as stream:
    rows = list(csv.DictReader(stream))
tasks = []
for row in rows:
    seeds = sorted((base/'isolated_library'/row['id']/'conformers').glob('conf_*.xyz'))
    assert seeds
    tasks.append(dict(key=row['id']+'__'+row['state'], row=row,
                      atoms=int(seeds[0].read_text().splitlines()[0]),
                      seed_sha256={str(p):uma.file_sha256(p) for p in seeds}))
pilot_id = 'study_aryl_c285539322fa'  # AV-3c, maximum-size tetracarboxy parent graft.
pilots = [i for i,t in enumerate(tasks) if t['row']['id']==pilot_id]
assert len(pilots)==3 and all(tasks[i]['atoms']==70 for i in pilots)
os.environ['FAIRCHEM_CACHE_DIR'] = '/nfs/lambda_stor_01/homes/rzhu/.cache/fairchem'
paths = [base/'worker.py',base/'isolated_library/manifest.csv',base/'prepared_config.py',
         Path(uma.__file__),ROOT/'src/redox/validation/protonation.py',
         ROOT/'src/redox/build/geometry.py',ROOT/'src/redox/core/common.py',
         ROOT/'config/project.json']
plan = dict(status='27_frozen_UMA_points_predeclared', tasks=tasks, pilot_indices=pilots,
    model_identity=uma.model_identity(uma.DEFAULT_MODEL), source_sha256={str(p):uma.file_sha256(p) for p in paths},
    max_atoms=70, allocation_fraction=.75, fmax_eV_A=.02, optimizer_steps=300,
    spin_and_conformer_scope='All prepared seeds and parity-consistent UMA multiplicities',
    resource_pilot_required_before_screen=True, resource_pilot_scope='70-atom AV-3c in all three charges on 16-GB T4; full 9-graph output retains per-state memory checks',
    proton_topology_gate=True, speciation_equilibrium_established=False,
    source_measurements_inherited=False, registered_candidates=0, solution_properties_computed=0,
    created=time.time())
atomic_write_json(base/'uma_plan.json', plan)
print('27 states frozen; pilot task indices:', pilots)
