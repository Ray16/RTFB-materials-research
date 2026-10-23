import json,os,sys
from pathlib import Path
rows=[];used=set();gpu_used=set()
for name in json.load(sys.stdin):
 p=Path(name);d=json.loads(p.read_text())
 if d.get('host')!='lambda13':continue
 proc=Path('/proc')/str(d['pid']);stat=proc/'stat'
 if not stat.exists() or stat.read_text().rsplit(')',1)[1].split()[19]!=d['process_start']:continue
 affinity=set(os.sched_getaffinity(d['pid']));caps={}
 for entry in (proc/'environ').read_bytes().split(bytes([0])):
  if b'=' not in entry:continue
  k,v=entry.split(b'=',1)
  if k.decode() in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS','CUDA_VISIBLE_DEVICES'):caps[k.decode()]=v.decode()
 if d.get('stage')=='thermal':
  expected={k:'2' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')};visible=caps.pop('CUDA_VISIBLE_DEVICES');assert visible==str(d['gpu']) and d['gpu'] not in gpu_used;gpu_used.add(d['gpu']);assert len(affinity)==2
 else:
  expected=d['thread_caps'];assert len(affinity)==14
 assert caps==expected
 assert not affinity&used,(name,sorted(affinity&used));used.update(affinity)
 rows.append(dict(source=name,pid=d['pid'],process_start=d['process_start'],CPU_affinity=sorted(affinity),kind='native_thermal' if d.get('stage')=='thermal' else 'CREST',thread_caps=caps))
print(json.dumps(dict(host=os.uname().nodename,actors=rows,disjoint_CPUs=len(used),live_native_GPUs=sorted(gpu_used))))
