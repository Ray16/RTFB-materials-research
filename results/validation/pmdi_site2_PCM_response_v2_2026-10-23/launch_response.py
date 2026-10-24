"""Launch one source-bound fixed-density PCM response/resource qualification."""
from pathlib import Path
import json,os,socket,subprocess,sys,time
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import resources
from redox.workflows import polymer_campaign as api
BASE=ROOT/'results/validation/pmdi_site2_PCM_response_v2_2026-10-23'
plan=api.read_plan(BASE)
if (BASE/'launch.json').exists():raise RuntimeError('Inspect the existing attempt before submitting again')
import pyscf,gpu4pyscf,numpy,scipy,ase
versions={m.__name__:m.__version__ for m in (pyscf,gpu4pyscf,numpy,scipy,ase)}
assert (pyscf.__version__,gpu4pyscf.__version__)==('2.14.0','1.8.1')
with resources.allocation_guard():
 capacity=resources.host_capacity(2,minimum_ram_gib=plan['acceptance']['minimum_host_RAM_GiB'])
 assert plan['acceptance']['minimum_host_RAM_GiB']*2**30<.75*capacity['available_RAM_bytes']
 api.verify_sources(plan['producer_sha256']);api.verify_sources(plan['input_sha256'])
 cpus=resources.cpu_allocation(2)
 slots=resources.free_slots(['lambda1']);gpu=next((g for h,g in slots if resources.gpu_suitability(g)['memory_MiB']>=31744),None)
 if gpu is None:raise RuntimeError('No suitable free V100 now; no molecular attempt started')
 preflight=dict(host=socket.gethostname(),capacity=capacity,CPU_candidate=cpus,device=resources.gpu_suitability(gpu),versions=versions,plan_sha256=api.digest(BASE/'plan.json'),checked=time.time())
 atomic_write_json(BASE/'preflight.json',preflight)
 env=dict(os.environ,PYTHONPATH=str(ROOT/'src'))
 for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[k]='2'
 argv=['taskset','--cpu-list',','.join(map(str,cpus)),str(resources.GATE),'run',str(gpu),'--label','redox_PCM_first_response_115','--',sys.executable,'-m','redox.workflows.gpu','--out',str(BASE/f'reserved_gpu{gpu}'),'--probe','--',sys.executable,'-u',str(ROOT/'scripts/validation/qualify_pcm_first_response_v2.py'),'run','--out',str(BASE)]
 with (BASE/'wrapper.log').open('x') as log:
  p=subprocess.Popen(argv,cwd=ROOT,env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
 proof=dict(host=socket.gethostname(),pid=p.pid,process_start=Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19],argv=argv,CPU_affinity=cpus,GPU=gpu,started=time.time(),plan_sha256=api.digest(BASE/'plan.json'),launcher_sha256=api.digest(__file__))
 atomic_write_json(BASE/'launch.json',proof)
print(json.dumps(proof))

# Register the actual standalone response process for the CPU residency guard.
for _ in range(100):
 status_path=BASE/'status.json'
 if status_path.exists():
  status=json.loads(status_path.read_text())
  directory=BASE/'observation/workers';directory.mkdir(parents=True,exist_ok=True)
  atomic_write_json(directory/'response.json',dict(host=status['host'],pid=status['pid'],process_start=status['process_start'],gpu=int(status['GPU']),status='computing',scope='Fixed-density response/resource diagnostic only'))
  break
 time.sleep(.1)
else:raise RuntimeError('Response startup unconfirmed; inspect the existing launch, do not duplicate it')
