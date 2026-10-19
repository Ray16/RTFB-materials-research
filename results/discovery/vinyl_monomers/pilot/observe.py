"""Read-only proof of actual private pilot processes, allocations and progress."""
import csv,json,os,socket,subprocess,time
from pathlib import Path
from redox.core.io import atomic_write_json
base=Path(__file__).resolve().parent
launch=json.loads((base/'launch_lambda13.json').read_text())
raw=subprocess.check_output(['nvidia-smi','--query-compute-apps=gpu_uuid,pid,used_gpu_memory','--format=csv,noheader,nounits'],text=True)
resident=list(csv.reader(raw.splitlines(),skipinitialspace=True))
records=[];used_cpus=set();used_uuids=set()
for entry in launch['workers']:
 host=entry['host'];gpu=entry['gpu'];path=base/'workers'/f'{host}_gpu{gpu}.json'
 if not path.exists():continue
 r=json.loads(path.read_text());pid=r['pid'];stat=Path(f'/proc/{pid}/stat')
 start=stat.read_text().rsplit(')',1)[1].split()[19]
 assert start==r['process_start']
 cpus=sorted(os.sched_getaffinity(pid));assert cpus==entry['cpu_affinity'] and not used_cpus.intersection(cpus)
 used_cpus.update(cpus)
 env=dict(s.split('=',1) for s in Path(f'/proc/{pid}/environ').read_bytes().decode().split('\0') if '=' in s)
 assert env['CUDA_VISIBLE_DEVICES']==str(gpu)
 for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS'):assert env[key]=='2'
 uuid=entry['device']['uuid'];assert uuid not in used_uuids;used_uuids.add(uuid)
 processes=[dict(pid=int(p),memory_MiB=int(m)) for u,p,m in resident if u==uuid]
 assert len(processes)==1 and processes[0]['pid']==pid,(gpu,processes,pid)
 claim=base/'claims'/r['active_key']/'owner.json';assert json.loads(claim.read_text())['pid']==pid
 status=json.loads((base/'tasks'/r['active_key']/'status.json').read_text())
 log=base/'tasks'/r['active_key']/(status['stage']+'.log')
 assert log.exists() and log.stat().st_size>0
 records.append(dict(gpu=gpu,uuid=uuid,pid=pid,process_start=start,cpu_affinity=cpus,stage=status['stage'],key=r['active_key'],resident=processes,log=str(log),log_bytes=log.stat().st_size,log_tail=log.read_text()[-600:]))
report=dict(status='actual_private_vinyl_quantum_starts_verified',host=socket.gethostname(),actual_GPU_workers=len(records),workers=records,load_1m=os.getloadavg()[0],online_cpus=os.cpu_count(),completed_state_records=len(list((base/'done').glob('*.json'))),failed_state_records=len(list((base/'failures').glob('*.json'))),expected_states=9,ground_spin_selection_complete=False,properties_published=False,broad_goal_resumed=False,checked=time.time())
assert len(records)==len(launch['workers'])
assert report['load_1m']<=report['online_cpus']
atomic_write_json(base/'start_verification.json',report)
print(json.dumps({k:v for k,v in report.items() if k!='workers'}))
