"""Inspect actual native pilot identity and reserved GPU/CPU resources."""
from pathlib import Path
import json,os,socket,subprocess,time
import psutil
R=Path('/nfs/lambda_stor_01/homes/rzhu/0_redox');C=R/'results/polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_2_native_xc_surface_budget'
s=json.loads((C/'thermal/workers/lambda2_gpu0_pilot.json').read_text());p=psutil.Process(s['pid']);start=Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19]
assert start==s['process_start'];e=p.environ();caps=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS');assert all(e.get(k)=='2' for k in caps);assert len(p.cpu_affinity())==2
owner=(R.parent/f'.gpu_locks/lambda2_gpu{e["CUDA_VISIBLE_DEVICES"]}/owner').read_text().split();assert int(owner[0])==p.ppid();assert owner[5]==Path(f'/proc/{p.ppid()}/stat').read_text().rsplit(')',1)[1].split()[19]
rows=[r.strip().split(', ') for r in subprocess.check_output(['nvidia-smi','--query-compute-apps=pid,gpu_uuid,used_memory','--format=csv,noheader,nounits'],text=True).splitlines() if r.strip()]
ours=[r for r in rows if int(r[0])==p.pid];assert len(ours)==1;assert len([r for r in rows if r[1]==ours[0][1]])==1
plan=json.loads((C/'plan.json').read_text());observer=json.loads((R/'results/validation/pmdi_occupancy_2026-10-23/baseline_coordinator.json').read_text());assert not set(p.cpu_affinity())&set(observer['CPU_affinity'])
proof=dict(checked=time.time(),host=socket.gethostname(),pid=p.pid,process_start=start,CPU_affinity=p.cpu_affinity(),thread_caps={k:e[k] for k in caps},GPU=e['CUDA_VISIBLE_DEVICES'],resident_GPU=ours[0],gate_parent_owner_verified=True,one_resident_process_on_assigned_GPU=True,disjoint_from_observer=True,native_host_budget_GiB=plan.get('native_host_budget_GiB'),qualified_native_G=False)
out=R/'results/validation/pmdi_site2_fine_force_2026-10-23/native_startup_verification.json';out.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))
