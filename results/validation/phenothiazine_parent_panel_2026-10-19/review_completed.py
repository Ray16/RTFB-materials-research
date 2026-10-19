"""Read-only independent review of finished states; never score incomplete couples."""
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import time

from redox.core.common import ROOT, DFT, apply_thermal_policy, free_energy
from redox.core.io import atomic_write_json
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity

base=Path(__file__).resolve().parent
out=base/'reviews'/str(time.time_ns());out.mkdir(parents=True,exist_ok=False)
with allocation_guard():
    resources=host_capacity(2);cpus=cpu_allocation(2);os.sched_setaffinity(0,cpus)
assert os.environ['OPENBLAS_NUM_THREADS']==os.environ['OMP_NUM_THREADS']=='2'
path=ROOT/'scripts/validation/audit/recompute_axes.py'
spec=importlib.util.spec_from_file_location('independent_parent_state_review',path)
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
plan=json.loads((base/'plan.json').read_text());checks=[];reviewed=[];pending=[]
for task in plan['tasks']:
    status_path=base/'tasks'/task['key']/'status.json'
    status=json.loads(status_path.read_text()) if status_path.exists() else {}
    if status.get('status')!='state_records_verified':
        pending.append(dict(key=task['key'],status=status.get('status','not_started'),stage=status.get('stage')))
        continue
    for name,digest in status['output_sha256'].items():
        checks.append(dict(key=task['key'],check='completed_output_sha256',path=name,
                           ok=hashlib.sha256(Path(name).read_bytes()).hexdigest()==digest))
    sd=DFT/task['id']/task['state']
    state=audit.state(task['id'],task['state'])
    checks.append(dict(key=task['key'],check='independent_surface_and_partition_admission',
                       ok=state is not None and state['G'] is not None
                          and state['thermal_definition']=='native_solution_dft_qrrho_spin'))
    raw=json.loads((sd/'result.json').read_text())
    primary,quality=free_energy(apply_thermal_policy(dict(raw,e_smd_eV=state['e_smd']),sd))
    checks.append(dict(key=task['key'],check='independent_total_G_eV',
                       ok=primary is not None and state['G'] is not None and abs(primary-state['G'])<1e-9))
    reviewed.append(dict(key=task['key'],G_eV=state['G'],thermal_eV=state['g_th'],
                         primary_thermal_quality=quality,alternate_spin_ground_state_proven=False))
report=dict(status='partial_finished_state_review',observed=time.time(),states_reviewed=len(reviewed),
            checks=len(checks),failed_checks=sum(not c['ok'] for c in checks),reviewed=reviewed,pending=pending,
            CPU_affinity=cpus,resources=resources,
            experimental_accuracy_validated=False,spin_ground_states_independently_compared=False,
            primary_table_promotion=False,
            next_gate='Complete the adjacent state pairs; check relevant alternate spins and graph integrity before strict experimental scoring',
            reviewer_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                             for p in (Path(__file__),path,path.with_name('native_partition.py'))})
atomic_write_json(out/'checks.json',checks);atomic_write_json(out/'summary.json',report)
atomic_write_json(base/'partial_review_current.json',dict(path=str(out.relative_to(ROOT)),**report))
print(json.dumps(dict(reviewed=len(reviewed),checks=len(checks),failed_checks=report['failed_checks'],path=str(out))))
assert report['failed_checks']==0
