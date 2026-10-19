"""Freeze registered inputs and unchanged electronic/thermal protocols before launch."""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

from redox.core.common import ROOT, read_manifest
from redox.core.io import atomic_write_json
from redox.core.protocol import ACTIVE_SP
from redox.qm.native_thermal import NATIVE_THERMAL, implementation_controls
from redox.qm.uma import cache_errors
from redox.validation.protonation import hydrogen_topology_status
from redox.workflows.resources import allocation_guard, cpu_allocation

base=Path(__file__).resolve().parent
with allocation_guard():os.sched_setaffinity(0,cpu_allocation(1))
assert not (base/'plan.json').exists()
screen=ROOT/'results/discovery/aryl_viologen_screen_2026-10-19'
adoption=json.loads((screen/'registration_prepared/adoption.json').read_text())
ids=set(adoption['registered_ids'])
assert len(ids)==9 and all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in adoption['canonical_sha256'].items())
rows=[r for r in read_manifest() if r['id'] in ids]
assert len(rows)==27
screen_plan=json.loads((screen/'uma_plan.json').read_text())
screen_review=json.loads((screen/'uma_completion.json').read_text())
assert screen_review['verified_states']==27 and screen_review['failed_states']==0
tasks=[]
for row in rows:
    gid,state=row['id'],row['state'];sd=ROOT/'calcs/uma'/gid/state
    seeds=sorted((ROOT/'library'/gid/'conformers').glob('conf_*.xyz'))
    result=json.loads((sd/'result.json').read_text())
    assert not cache_errors(result,row,sd/'relaxed.xyz',seeds,screen_plan['model_identity'],.02)
    assert hydrogen_topology_status(row['smiles'],sd/'relaxed.xyz')=='ok'
    assert not (ROOT/'calcs/dft'/gid/state).exists()
    paths=seeds+[sd/'result.json',sd/'relaxed.xyz']
    tasks.append(dict(key=gid+'__'+state,id=gid,state=state,manifest_row=row,
                      input_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}))
# The previous qualified campaign enumerated the numerical producers and their
# external derivative/partition controls. Reuse that scope with current hashes.
old=json.loads((ROOT/'results/validation/cata_solution_2026-10-19/plan.json').read_text())
assert ACTIVE_SP==old['protocol'] and NATIVE_THERMAL==old['thermal_protocol']
implementation_controls()
paths={Path(p) for p in old['producer_sha256'] if p!=str(ROOT/'results/validation/cata_solution_2026-10-19/worker.py')}
paths.update([base/'worker.py',base/'queue_worker.py',ROOT/'src/redox/validation/protonation.py',ROOT/'src/redox/workflows/checkpoints.py',base/'comparison_plan.json'])
plan=dict(status='27_registered_source_bound_solution_states_prepared',tasks=tasks,
    producer_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(paths)},
    protocol=ACTIVE_SP,thermal_protocol=NATIVE_THERMAL,
    optimization_seeds=old['optimization_seeds'],native_Hessian_response_memory_GiB=24,
    spin_candidates_from_UMA=[c['key'] for c in screen_review['cases'] if c['DFT_spin_comparison_required']],
    native_minimum_repair_scheduled=False,attempt_limit_per_state=1,
    shared_resumable_queue=True,per_worker_consecutive_failure_limit=3,
    hydrogen_topology_gate_after_solution_relaxation=True,physical_thresholds_changed=False,
    parent_measurements_transferred=False,equilibrium_speciation_established=False,
    independent_admission_required=True,git_head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),created=time.time())
atomic_write_json(base/'plan.json',plan)
for name in ('claims','done','failures','workers','tasks'):(base/name).mkdir(exist_ok=True)
print('Prepared 27 new solution states with unchanged native protocol and five predeclared comparisons')
