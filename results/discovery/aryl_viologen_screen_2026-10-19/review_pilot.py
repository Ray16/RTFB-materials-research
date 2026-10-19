"""Independently verify all three 70-atom T4 pilot outputs and terminal owners."""
import hashlib
import json
import time
from pathlib import Path

from redox.core.io import atomic_write_json
from redox.qm import uma
from redox.validation.protonation import hydrogen_topology_status

base=Path(__file__).resolve().parent
plan=json.loads((base/'uma_plan.json').read_text())
observation=json.loads((base/'lambda5_observation.json').read_text())
launch=json.loads((base/'pilot_lambda5_launch.json').read_text())
pilot_pids={w['gate_pid'] for w in launch['workers']}
pilot_gates=[g for g in observation['gates'] if g['pid'] in pilot_pids]
assert len(pilot_gates)==3 and all(g['terminal'] for g in pilot_gates), 'Original pilot gates still live'
for path,digest in plan['source_sha256'].items():
    assert uma.file_sha256(path)==digest
cases=[]
for index in plan['pilot_indices']:
    task=plan['tasks'][index]
    report=json.loads((base/'tasks'/task['key']/'status.json').read_text())
    assert report['status']=='gas_UMA_state_and_hydrogen_topology_verified'
    assert report['task']==task and report['device']=='Tesla T4'
    assert 14*1024**3<=report['total_memory_bytes']<=17*1024**3
    assert task['atoms']==70 and report['allocation_fraction']==.75
    peak=max(s['used_MiB'] for s in report['memory_samples'])
    assert peak==report['peak_device_used_MiB'] and peak*1024**2<.8*report['total_memory_bytes']
    assert report['peak_reserved_bytes']<.8*report['total_memory_bytes']
    row=task['row'];directory=base/'uma'/row['id']/row['state']
    result=json.loads((directory/'result.json').read_text())
    for path,digest in task['seed_sha256'].items():assert uma.file_sha256(path)==digest
    assert not uma.cache_errors(result,row,directory/'relaxed.xyz',list(map(Path,task['seed_sha256'])),plan['model_identity'],.02)
    assert hydrogen_topology_status(row['smiles'],directory/'relaxed.xyz')=='ok'
    assert report['result_sha256']==uma.file_sha256(directory/'result.json')
    assert report['geometry_sha256']==uma.file_sha256(directory/'relaxed.xyz')
    cases.append(dict(key=task['key'],charge=int(row['charge']),peak_device_used_MiB=peak,
        peak_reserved_bytes=report['peak_reserved_bytes'],energy_eV=result['energy_eV'],
        chosen_mult=result['chosen_mult'],spin_gap_eV=result['spin_gap_eV'],fmax=result['fmax'],
        n_conf=result['n_conf'],result_sha256=report['result_sha256'],geometry_sha256=report['geometry_sha256']))
review=dict(status='three_charge_70_atom_T4_resource_pilot_verified',cases=cases,
    terminal_gates=pilot_gates,terminal_observation_checked=observation['checked'],
    plan_sha256=hashlib.sha256((base/'uma_plan.json').read_bytes()).hexdigest(),
    maximum_device_memory_MiB=max(c['peak_device_used_MiB'] for c in cases),
    permits='UMA-s-1p2p1 OMol gas preoptimization for the declared 9-graph, 70-atom aryl-viologen panel; every state retains independent memory and topology gates',
    DFT_or_Hessian_on_16_GB_permitted=False,solution_accuracy_validated=False,checked=time.time())
atomic_write_json(base/'pilot_review.json',review)
print(json.dumps(review))
