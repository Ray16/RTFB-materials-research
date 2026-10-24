"""Launch the two predeclared small spin controls with atomic GPU reservations."""
from pathlib import Path
import json,os,socket,subprocess,sys,time
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import resources
from redox.workflows import polymer_campaign as api
BASE=ROOT/'results/validation/pcm_first_response_workspace_2026-10-23'
if (BASE/'launch.json').exists():raise RuntimeError('Inspect the existing launch before another attempt')
import pyscf,gpu4pyscf,numpy,scipy,ase,geometric
versions={m.__name__:m.__version__ for m in (pyscf,gpu4pyscf,numpy,scipy,ase,geometric)}
assert (pyscf.__version__,gpu4pyscf.__version__)==('2.14.0','1.8.1')
files=[ROOT/'scripts/pipeline/pcm_first_derivative_blocks.py',ROOT/'scripts/pipeline/smd_first_response_workspace.py',ROOT/'scripts/validation/qualify_first_response_workspace.py',Path(__file__).resolve()]
proof=dict(host=socket.gethostname(),started=time.time(),versions=versions,workspace_GiB=28,source_sha256={str(p):api.digest(p) for p in files},workers=[],purpose='Original-versus-modified full energy/coordinate/gradient/Hessian small spin controls; no candidate state')
with resources.allocation_guard():
 proof['capacity']=resources.host_capacity(4,minimum_ram_gib=16)
 slots=resources.free_slots(['lambda1'])
 suitable=[g for h,g in slots if resources.gpu_suitability(g)['memory_MiB']>=31744]
 if len(suitable)<2:raise RuntimeError('Two free suitable V100s required; no controls started')
 for mult,gpu in zip((1,2),suitable[:2]):
  cpus=resources.cpu_allocation(2)
  env=dict(os.environ,PYTHONPATH=str(ROOT/'src'))
  for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[k]='2'
  out=BASE/'controls'/str(mult)
  argv=['taskset','--cpu-list',','.join(map(str,cpus)),str(resources.GATE),'run',str(gpu),'--label',f'redox_PCM_first_blocks_mult{mult}','--',sys.executable,'-m','redox.workflows.gpu','--out',str(BASE/f'reserved_gpu{gpu}_spin{mult}'),'--probe','--',sys.executable,'-u',str(ROOT/'scripts/validation/qualify_first_response_workspace.py'),'--mult',str(mult),'--out',str(out),'--workspace-GiB','28']
  with (BASE/f'wrapper_gpu{gpu}_spin{mult}.log').open('x') as log:
   p=subprocess.Popen(argv,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
  for _ in range(100):
   if p.poll() is not None:raise RuntimeError('Reserved control wrapper stopped during launch; inspect preserved outputs')
   if os.sched_getaffinity(p.pid)==set(cpus):break
   time.sleep(.01)
  else:raise RuntimeError('Cannot confirm submitted CPU mask')
  proof['workers'].append(dict(mult=mult,gpu=gpu,device=resources.gpu_suitability(gpu),CPU_affinity=cpus,gate_pid=p.pid,process_start=Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19],argv=argv))
  atomic_write_json(BASE/'launch.json',proof)
print(json.dumps(proof))
