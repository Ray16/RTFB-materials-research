"""Start one finite CPU-only residency observer on a worker host."""
from pathlib import Path
import json,os,socket,subprocess,sys,time
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
BASE=ROOT/'results/validation/gpu_exclusivity_2026-10-23';host=socket.gethostname().split('.')[0]
record=BASE/f'launch_{host}.json'
if record.exists():raise RuntimeError('Verify the existing observer before another launch')
owner=json.loads((ROOT/'results/validation/pmdi_occupancy_2026-10-23/baseline_coordinator.json').read_text())
argv=[sys.executable,'-u',str(ROOT/'scripts/fleet/watch_gpu_exclusivity.py'),'--work-root',str(ROOT/'results/polymer_models/campaigns/2026-10-23_crest_refinement'),'--out',str(BASE/host),'--deadline',str(owner['deadline']),'--poll-seconds','10']
with (BASE/f'{host}.log').open('x') as log:
 p=subprocess.Popen(argv,cwd=ROOT,env=dict(os.environ,PYTHONPATH=str(ROOT/'src')),stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
proof=dict(host=host,pid=p.pid,process_start=Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19],argv=argv,started=time.time(),deadline=owner['deadline'],observer_type='CPU-only; no GPU reservation or CUDA context',source_sha256=__import__('hashlib').sha256((ROOT/'scripts/fleet/watch_gpu_exclusivity.py').read_bytes()).hexdigest())
atomic_write_json(record,proof);print(json.dumps(proof))
