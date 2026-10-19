"""Prepare the initial vinyl-monomer pilot in a private conformer library."""
import csv
import hashlib
import json
import os
from pathlib import Path
import pprint
import runpy
import sys

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

BASE=Path(__file__).resolve().parent
proposal=json.loads((BASE.parent/'proposals.json').read_text())
review=json.loads((BASE.parent/'independent_review.json').read_text());assert review['failed_checks']==0
snapshot=Path(proposal['source_snapshot']);source={}
for name in ('original39.csv','additional_candidates.csv','comparison_controls.csv'):
    with (snapshot/name).open() as stream:
        for row in csv.DictReader(stream):source[row['id']]=row
groups=[];targets={}
for row in proposal['records']:
    if not row['initial_modeling_pilot']:continue
    original=source[row['source_capped_id']]
    groups.append(dict(id=row['id'],name=row['name'],family=row['family'],frag=original['frag'],states=row['states']))
    targets[row['id']]=row['smiles']
assert len(groups)==3 and not (BASE/'isolated_library').exists()
digest=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
manifest_hash=digest(ROOT/'library/manifest.csv')
with allocation_guard():resources=host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
config=BASE/'prepared_config.py'
config.write_text('SCAFFOLD = '+repr(proposal['scaffold'])+'\nGROUPS = '+pprint.pformat(groups,width=100,sort_dicts=False)+'\n')
sys.argv=['redox.build.candidates','--config',str(config),'--library',str(BASE/'isolated_library')]
runpy.run_module('redox.build.candidates',run_name='__main__')
from rdkit import Chem
from redox.qm.uma import geometry_checks
from redox.validation.protonation import hydrogen_topology_status
with (BASE/'isolated_library/manifest.csv').open() as stream:states=list(csv.DictReader(stream))
checks=[]
for row in states:
    assert Chem.MolToSmiles(Chem.MolFromSmiles(row['smiles']))==Chem.MolToSmiles(Chem.MolFromSmiles(targets[row['id']]))
    for seed in sorted((BASE/'isolated_library'/row['id']/'conformers').glob('conf_*.xyz')):
        assert not geometry_checks(row['smiles'],seed,int(row['charge']),int(row['mult_hint']))
        assert hydrogen_topology_status(row['smiles'],seed)=='ok'
        checks.append(dict(id=row['id'],state=row['state'],seed=str(seed),sha256=digest(seed)))
assert len(states)==9 and digest(ROOT/'library/manifest.csv')==manifest_hash
record=dict(status='three_vinyl_graphs_and_nine_state_seed_sets_prepared',states=9,molecules=3,seed_state_checks=len(checks),checks=checks,
    shared_registry_unchanged=True,quantum_jobs_submitted=0,solution_properties_computed=False,resources=resources,
    source_sha256={str(p):digest(p) for p in (Path(__file__),config,BASE.parent/'proposals.json',BASE.parent/'independent_review.json')})
atomic_write_json(BASE/'preparation.json',record)
print(json.dumps({k:v for k,v in record.items() if k not in ('checks','resources','source_sha256')}))
