from pathlib import Path
import json,os,socket,subprocess,time,sys
R=Path('/nfs/lambda_stor_01/homes/rzhu/0_redox')
C=R/'results/validation/pmdi_site2_fine_force_2026-10-23'
from redox.workflows import resources
from redox.workflows import polymer_campaign as api
from redox.core.io import atomic_write_json
import psutil
P=json.loads((C/'plan.json').read_text())
with resources.allocation_guard():
 if (C/'launch.json').exists():raise RuntimeError('An existing launch must be inspected before another attempt')
 failure=json.loads((Path(P['source_campaign'])/'workspace_energy_status.json').read_text())
 try:
  source=psutil.Process(failure['pid'])
  start=Path(f'/proc/{source.pid}/stat').read_text().rsplit(')',1)[1].split()[19]
  if start==failure['process_start']:raise RuntimeError('Original scientific process is still live')
 except psutil.NoSuchProcess:pass
 capacity=resources.host_capacity(2,minimum_ram_gib=P['host_budget_GiB'])
 if P['host_budget_GiB']*2**30>.75*capacity['available_RAM_bytes']:raise RuntimeError('Conservative host-RAM forecast does not fit')
 cpus=resources.cpu_allocation(2);os.sched_setaffinity(0,cpus)
 import pyscf,gpu4pyscf,numpy,scipy,ase,geometric
 versions={m.__name__:m.__version__ for m in (pyscf,gpu4pyscf,numpy,scipy,ase,geometric)}
 if (pyscf.__version__,gpu4pyscf.__version__)!=('2.14.0','1.8.1'):raise RuntimeError('Native software versions differ')
 api.verify_sources(P['producer_sha256']);api.verify_sources(P['input_sha256'])
 gpu=next((g for host,g in resources.free_slots(['lambda2']) if resources.gpu_suitability(g)['memory_MiB']>=P['minimum_GPU_MiB']),None)
 if gpu is None:raise RuntimeError('No suitable free V100 available now; no molecular attempt started')
 device=resources.gpu_suitability(gpu)
 preflight=dict(host=socket.gethostname(),capacity=capacity,CPU_candidate=cpus,device=device,versions=versions,plan_sha256=api.digest(C/'plan.json'),checked=time.time(),source_worker_stopped=True)
 atomic_write_json(C/'preflight.json',preflight)
 env=dict(os.environ)
 for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):env[k]='2'
 argv=['taskset','-c',','.join(map(str,cpus)),str(resources.GATE),'run',str(gpu),'--label','redox_pmdi_site2_fine_force','--',sys.executable,'-u',str(R/'scripts/validation/polymer_gradient_profile.py'),'run','--campaign',str(C)]
 with (C/'wrapper.log').open('x') as log:
  process=subprocess.Popen(argv,cwd=R,env=env,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
 start=Path(f'/proc/{process.pid}/stat').read_text().rsplit(')',1)[1].split()[19]
 launch=dict(host=socket.gethostname(),pid=process.pid,process_start=start,argv=argv,CPU_affinity=cpus,GPU=gpu,started=time.time(),plan_sha256=api.digest(C/'plan.json'),launcher_sha256=api.digest(C/'launch_source.py'))
 atomic_write_json(C/'launch.json',launch)
 print(json.dumps(launch))
