from pathlib import Path
import json,socket,os,time,psutil
R=Path('/nfs/lambda_stor_01/homes/rzhu/0_redox');host=socket.gethostname().split('.')[0];base=R/'results/validation/gpu_exclusivity_2026-10-23';records=[]
for name in (host,host+'_production'):
 folder=base/name
 if not folder.exists():continue
 row=json.loads((folder/'owner.json').read_text());p=psutil.Process(row['pid']);assert Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19]==row['process_start'];assert len(p.cpu_affinity())==1
 e=p.environ();keys=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS');assert all(e.get(k)=='1' for k in keys)
 state=json.loads((folder/'state.json').read_text());assert time.time()-state['checked']<60
 records.append(dict(name=name,pid=p.pid,process_start=row['process_start'],CPU_affinity=p.cpu_affinity(),native_thread_caps={k:e[k] for k in keys},last_audit=state))
proof=dict(checked=time.time(),host=host,live_observers=records)
(base/f'verified_{host}_snapshot.json').write_text(json.dumps(proof,indent=2)+'\n');print(json.dumps(dict(host=host,live_observers=len(records))))
