from pathlib import Path
import json,os,socket,subprocess,time
import psutil
R=Path('/nfs/lambda_stor_01/homes/rzhu/0_redox');C=R/'results/validation/pmdi_site2_fine_force_2026-10-23'
s=json.loads((C/'status.json').read_text());p=psutil.Process(s['pid'])
start=Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19]
assert start==s['process_start']
assert p.cpu_affinity()==s['CPU_affinity']==[4,5]
keys=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
e=p.environ();assert all(e.get(k)=='2' for k in keys)
assert e.get('CUDA_VISIBLE_DEVICES')==s['GPU']=='2'
owner=(R.parent/'.gpu_locks/lambda2_gpu2/owner').read_text().split();assert int(owner[0])==p.ppid();assert owner[5]==Path(f'/proc/{p.ppid()}/stat').read_text().rsplit(')',1)[1].split()[19]
raw=subprocess.check_output(['nvidia-smi','--query-compute-apps=pid,gpu_uuid,used_gpu_memory','--format=csv,noheader,nounits'],text=True)
rows=[line.strip().split(', ') for line in raw.splitlines() if line.strip()]
ours=[r for r in rows if int(r[0])==p.pid];assert len(ours)==1
assert len([r for r in rows if r[1]==ours[0][1]])==1
report=json.loads((R/'calcs/validation/pmdi_site2_fine_force_2026-10-23/diagnostic/status.json').read_text())
assert report['evaluations'] and report['evaluations'][0]['name']=='baseline'
proof=dict(checked=time.time(),host=socket.gethostname(),pid=p.pid,process_start=start,CPU_affinity=p.cpu_affinity(),native_thread_caps={k:e[k] for k in keys},GPU=e['CUDA_VISIBLE_DEVICES'],resident_GPU=ours[0],one_process_on_assigned_GPU=True,gate_parent_owner_verified=True,completed_diagnostic_evaluations=len(report['evaluations']),scientific_admission_changed=False)
(C/'startup_verification.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(proof))
