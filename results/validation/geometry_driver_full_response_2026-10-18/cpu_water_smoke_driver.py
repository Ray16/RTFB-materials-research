import os,json,time,hashlib
from pathlib import Path
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
from redox.core.io import atomic_write_json
from redox.qm.dft import dft_smd
base=Path('/nfs/lambda_stor_01/homes/rzhu/0_redox/results/validation/geometry_driver_full_response_2026-10-18/cpu_water_smoke')
base.mkdir(exist_ok=False)
with allocation_guard():
 resources=host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
seed=base/'seed.xyz';seed.write_text('3\nwater geometry driver smoke\nO 0 0 0\nH .90 0 0\nH -.20 .90 0\n')
started=time.time()
try:
 result=dft_smd(seed,0,1,do_smd=False,do_gas=False,do_freq=False,opt_out=base/'opt.xyz',max_opt_steps=30,backend='cpu')
 assert result['optimized'] and result['geometry_optimization_converged'] and result['optimization_SCF_converged'] and result['opt_gradient_grid_response']
 report=dict(status='actual_cpu_geometry_driver_passed',result=result,resources=resources,affinity=sorted(os.sched_getaffinity(0)),started=started,finished=time.time(),dft_sha256=hashlib.sha256(Path('/nfs/lambda_stor_01/homes/rzhu/0_redox/src/redox/qm/dft.py').read_bytes()).hexdigest(),GPU_jobs_submitted=0,scientific_battery_properties_validated=False)
except BaseException as exc:
 atomic_write_json(base/'result.json',dict(status='failed',error=repr(exc),finished=time.time()));raise
atomic_write_json(base/'result.json',report)
print(report['status'])
