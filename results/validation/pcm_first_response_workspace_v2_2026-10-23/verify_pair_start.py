"""Read-only verification of the actual resumed two-graft energy worker."""
from pathlib import Path
import json,os,socket,subprocess,time
import psutil
R=Path('/nfs/lambda_stor_01/homes/rzhu/0_redox');C=R/'results/polymer_models/campaigns/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_2_3_xc_checkpoint_resume'
pth=C/'energy/workers/lambda13_gpu3_queue.json';s=json.loads(pth.read_text());p=psutil.Process(s['pid']);start=Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19];assert start==s['process_start'];assert p.cpu_affinity()==s['cpus']
e=p.environ();caps=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS');assert all(e.get(k)=='2' for k in caps)
gpu=e['CUDA_VISIBLE_DEVICES'];owner=(R.parent/f'.gpu_locks/lambda13_gpu{gpu}/owner').read_text().split();assert int(owner[0])==p.ppid();assert owner[5]==Path(f'/proc/{p.ppid()}/stat').read_text().rsplit(')',1)[1].split()[19]
rows=[r.strip().split(', ') for r in subprocess.check_output(['nvidia-smi','--query-compute-apps=pid,gpu_uuid,used_memory','--format=csv,noheader,nounits'],text=True).splitlines() if r.strip()];ours=[r for r in rows if int(r[0])==p.pid];assert len(ours)==1;assert len([r for r in rows if r[1]==ours[0][1]])==1
raw=R/'calcs/polymer_models/2026-10-23_crest_refinement/baseline__design_pmdi_ethyl__occ_2_3_xc_checkpoint_resume/energy/baseline__design_pmdi_ethyl__occ_2_3__pose0/rest/geometry_checkpoints';last=json.loads((raw/'latest.json').read_text());assert last['SCF_converged'] and not last['solution_geometry_qualified']
proof=dict(checked=time.time(),host=socket.gethostname(),pid=p.pid,process_start=start,CPU_affinity=p.cpu_affinity(),native_thread_caps={k:e[k] for k in caps},GPU=gpu,resident_GPU=ours[0],gate_parent_owner_verified=True,one_resident_process_on_assigned_GPU=True,latest_unqualified_geometry_evaluation=last,qualified_E=False,qualified_G=False)
out=R/'results/validation/pcm_first_response_workspace_v2_2026-10-23/pair_optimizer_start_verification.json';out.write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))
