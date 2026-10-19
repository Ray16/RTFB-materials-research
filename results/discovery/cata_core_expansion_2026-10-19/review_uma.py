"""Inspect completed warm-start contracts and actual terminal processes on lambda13."""
import hashlib,json,os,time
from pathlib import Path
from redox.qm import uma
from redox.core.io import atomic_write_json
base=Path(__file__).resolve().parent;plan=json.loads((base/'uma_plan.json').read_text());launch=json.loads((base/'uma_launch.json').read_text())
os.environ['FAIRCHEM_CACHE_DIR']=plan['model_cache']
assert uma.model_identity(uma.DEFAULT_MODEL)==plan['model_identity']
for path,digest in plan['source_sha256'].items():assert uma.file_sha256(path)==digest,path
records=[];processes={}
for task,allocation in zip(plan['tasks'],launch['workers']):
 row=task['row'];key=row['id']+'__'+row['state'];report=json.loads((base/'uma_workers'/key/'status.json').read_text())
 execution=json.loads((base/f"uma_reserved_{allocation['task']}"/'execution.json').read_text())
 assert report['status']=='gas_UMA_state_verified' and execution['status']=='process_completed' and execution['returncode']==0
 directory=base/'uma'/row['id']/row['state'];result=json.loads((directory/'result.json').read_text())
 seeds=list(map(Path,task['seed_sha256']))
 for path,digest in task['seed_sha256'].items():assert uma.file_sha256(path)==digest
 assert not uma.cache_errors(result,row,directory/'relaxed.xyz',seeds,plan['model_identity'],.02)
 assert uma.file_sha256(directory/'result.json')==report['result_sha256']
 assert report['peak_reserved_bytes']<.8*report['total_memory_bytes']
 for pid in [report['pid'],allocation['controller_pid']]:
  path=Path(f'/proc/{pid}/stat');state=path.read_text().rsplit(')',1)[1].split()[0] if path.exists() else None
  processes[str(pid)]=state;assert state in (None,'Z'),'Live process remains '+str(pid)
 records.append(dict(id=row['id'],state=row['state'],energy_eV=result['energy_eV'],chosen_mult=result['chosen_mult'],spin_gap_eV=result['spin_gap_eV'],converged=result['converged'],fmax_eV_A=result['fmax'],n_conf=result['n_conf'],model_identity=result['model_identity'],result_sha256=report['result_sha256'],geometry_sha256=report['geometry_sha256'],peak_reserved_bytes=report['peak_reserved_bytes'],DFT_alternate_spin_required=result['spin_gap_eV'] is not None and result['spin_gap_eV']<.5))
report=dict(status='all_six_gas_UMA_states_verified_terminal',cases=records,actual_starts=6,complete_UMA_states=6,failed_states=0,canonical_registry_unchanged=True,registered_candidates=0,solution_potentials_computed=0,solution_kinetics_validated=False,terminal_process_states=processes,reviewed=time.time(),note='Warm starts only; UMA spin gaps do not select a DFT/native ground state. No parent label transfers.')
atomic_write_json(base/'uma_completion.json',report);print(json.dumps({k:v for k,v in report.items() if k not in ('cases','terminal_process_states')}))
