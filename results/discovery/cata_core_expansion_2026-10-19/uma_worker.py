"""One isolated, source-bound UMA state; one actual attempt, no blind retries."""
import argparse,csv,fcntl,json,os,signal,socket,sys,time
from pathlib import Path
from redox.core.io import atomic_write_json
from redox.workflows.resources import verify_gate_owner
from redox.qm import uma

base=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('--task',type=int,required=True);args=parser.parse_args()
plan=json.loads((base/'uma_plan.json').read_text());task=plan['tasks'][args.task]
key=task['row']['id']+'__'+task['row']['state'];folder=base/'uma_workers'/key;folder.mkdir(parents=True,exist_ok=True)
lock=(folder/'writer.lock').open('a+');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
report=dict(status='starting',key=key,host=socket.gethostname(),pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],CPU_affinity=sorted(os.sched_getaffinity(0)),started=time.time(),attempts=1,solution_properties_validated=False)
def save():atomic_write_json(folder/'status.json',report)
def verify():
 for p,h in {**plan['source_sha256'],**task['seed_sha256']}.items():
  if uma.file_sha256(p)!=h:raise RuntimeError('Prepared UMA source/input changed: '+p)
def interrupted(signum,frame):raise KeyboardInterrupt('signal '+str(signum))
signal.signal(signal.SIGTERM,interrupted)
save()
try:
 verify_gate_owner();verify()
 import torch
 if torch.cuda.device_count()!=1:raise RuntimeError('Expected exactly one reserved device')
 device=torch.cuda.get_device_properties(0)
 torch.cuda.set_per_process_memory_fraction(.75,0);torch.cuda.reset_peak_memory_stats(0)
 if uma.model_identity(uma.DEFAULT_MODEL)!=plan['model_identity']:raise RuntimeError('Model identity changed')
 uma.LIBRARY=base/'isolated_library';uma.OUT=base/'uma'
 with (uma.LIBRARY/'manifest.csv').open() as stream:rows=list(csv.DictReader(stream))
 uma.read_manifest=lambda:rows
 report.update(status='computing',physical_gpu=os.environ['CUDA_VISIBLE_DEVICES'],device=device.name,total_memory_bytes=device.total_memory,model_identity=plan['model_identity'],native_thread_caps={key:os.environ.get(key) for key in ['OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS']});save()
 sys.argv=['redox.qm.uma','--only',task['row']['id']+':'+task['row']['state'],'--device','cuda']
 uma.main()
 directory=uma.OUT/task['row']['id']/task['row']['state'];result=json.loads((directory/'result.json').read_text())
 errors=uma.cache_errors(result,task['row'],directory/'relaxed.xyz',list(map(Path,task['seed_sha256'])),plan['model_identity'],.02)
 if errors:raise RuntimeError('UMA output not qualified: '+str(errors))
 verify();torch.cuda.synchronize()
 report.update(status='gas_UMA_state_verified',result=str(directory/'result.json'),result_sha256=uma.file_sha256(directory/'result.json'),geometry_sha256=uma.file_sha256(directory/'relaxed.xyz'),chosen_mult=result['chosen_mult'],spin_gap_eV=result['spin_gap_eV'],energy_eV=result['energy_eV'],peak_allocated_bytes=torch.cuda.max_memory_allocated(0),peak_reserved_bytes=torch.cuda.max_memory_reserved(0))
 if report['peak_reserved_bytes']>.8*device.total_memory:raise RuntimeError('Insufficient measured memory headroom')
except BaseException as exc:
 report.update(status='failed_preserving_outputs',error=repr(exc));raise
finally:
 report['finished']=time.time();save();lock.close()
