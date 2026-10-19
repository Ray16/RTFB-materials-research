"""Preflight and launch all six independent prepared states on free suitable GPUs."""
import csv,json,os,subprocess,sys,time
from pathlib import Path
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.qm import uma
from redox.workflows.resources import allocation_guard,cpu_allocation,free_slots,gpu_suitability,host_capacity
base=Path(__file__).resolve().parent
os.environ['FAIRCHEM_CACHE_DIR']='/nfs/lambda_stor_01/homes/rzhu/.cache/fairchem'
if (base/'uma_launch.json').exists():raise RuntimeError('Existing launch preserved; inspect actual owners before any resume')
with (base/'isolated_library/manifest.csv').open() as f:rows=list(csv.DictReader(f))
assert len(rows)==6
with allocation_guard():
 resources=host_capacity(12,minimum_ram_gib=64)
 slots=[(h,g,gpu_suitability(g)) for h,g in free_slots(['lambda13'])]
 if len(slots)<len(rows):raise RuntimeError('Six suitable free devices required for this six-state batch')
 cpus=cpu_allocation(12)
 model=uma.model_identity(uma.DEFAULT_MODEL)
 sources=[base/'uma_worker.py',base/'isolated_library/manifest.csv',Path(uma.__file__),ROOT/'src/redox/core/common.py',ROOT/'src/redox/core/io.py',ROOT/'src/redox/build/structures.py',ROOT/'config/project.json']
 tasks=[dict(row=row,seed_sha256={str(p):uma.file_sha256(p) for p in sorted((base/'isolated_library'/row['id']/'conformers').glob('conf_*.xyz'))}) for row in rows]
 atomic_write_json(base/'uma_plan.json',dict(tasks=tasks,source_sha256={str(p):uma.file_sha256(p) for p in sources},model_identity=model,scientific_protocol='UMA-s-1p2p1/OMol gas; all prepared conformers; parity-consistent singlet/triplet or doublet scan; fmax .02 eV/A; 300 steps',allocation_fraction=.75,registered_candidates=0,model_cache=os.environ['FAIRCHEM_CACHE_DIR'],source_parent_measurements_transferred=False,created=time.time()))
 record=dict(status='submitting',resources=resources,workers=[],actual_outputs_required=6,quantum_solution_jobs_submitted=0)
 atomic_write_json(base/'uma_launch.json',record)
 env=dict(os.environ,PYTHONPATH=str(ROOT/'src'))
 for key in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS']:env[key]='2'
 for i,((host,gpu,device),task) in enumerate(zip(slots,rows)):
  assigned=cpus[2*i:2*i+2]
  argv=[sys.executable,'-u','-m','redox.workflows.gpu','--gpu',str(gpu),'--cpus',','.join(map(str,assigned)),'--threads','2','--out',str(base/f'uma_reserved_{i}'),'--',sys.executable,'-u',str(base/'uma_worker.py'),'--task',str(i)]
  with (base/f'uma_controller_{i}.log').open('w') as output:
   proc=subprocess.Popen(argv,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,stdout=output,stderr=subprocess.STDOUT,start_new_session=True,preexec_fn=lambda assigned=assigned:os.sched_setaffinity(0,assigned))
  worker=dict(host=host,gpu=gpu,device=device,cpus=assigned,task=i,id=task['id'],state=task['state'],controller_pid=proc.pid,process_start=Path(f'/proc/{proc.pid}/stat').read_text().rsplit(')',1)[1].split()[19],argv=argv,submitted=time.time())
  record['workers'].append(worker);atomic_write_json(base/'uma_launch.json',record)
  time.sleep(1)
 record['status']='six_controllers_submitted';atomic_write_json(base/'uma_launch.json',record)
print(json.dumps(record))
