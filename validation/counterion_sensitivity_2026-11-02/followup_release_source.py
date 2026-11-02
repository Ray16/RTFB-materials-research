"""One bounded release of the six prepared exact chloride geometry controls."""
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import load, verify
from redox.workflows.resources import allocation_guard, host_capacity, cpu_allocation, GATE

with allocation_guard():
    host_capacity(1, minimum_ram_gib=4)
    os.sched_setaffinity(0, cpu_allocation(1))
root = Path(__file__).parent
output = root / 'followup_release'
output.mkdir(exist_ok=False)
review = ROOT / 'validation/counterion_sensitivity_2026-11-02/pilot_geometry_review.json'
pilot = json.loads(review.read_text())
if (pilot['status'] != 'primary_counterion_geometry_independently_checked' or pilot['key'] != 'N0_above'
        or pilot['native_G_qualified'] or pilot['nao'] != 412
        or any(sha256_file(p) != h for p, h in pilot['source_sha256'].items())):
    raise ValueError('Exact independently checked geometry pilot required')
rows = json.loads((root/'manifest.json').read_text())['states']
capacity = load(ROOT/'scripts/pipeline/run_namespace_energy_queue.py', 'chloride_release_capacity')
storage = load(ROOT/'scripts/fleet/run_charged_polymer_geometry_panel.py', 'chloride_release_storage')
owners = load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py', 'chloride_release_owners')
peers = json.loads((ROOT/'validation/corrected_sampling_2026-11-02/current_storage_peers.json').read_text())
policy = GATE.parent.parent/'.config/gpu_reserve/config.env'
hosts = re.search(r'^ALLOWED_HOSTS="([^"]+)"', policy.read_text(), re.M).group(1).split()
slots, errors = capacity.discover_slots(hosts)
pending = [r for r in rows if r['key'] != 'N0_above']
record = dict(status='chloride_geometry_followup_bounded_release', checked=time.time(), availability_errors=errors,
              launches=[], refusals=[], native_G_qualified=False,
              source_sha256={str(p):sha256_file(p) for p in (Path(__file__), review, root/'manifest.json',
                  Path(capacity.__file__), Path(storage.__file__), policy)})
remote = '''import json,os,pathlib,subprocess,sys
from redox.workflows.screening_control import probe
from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
from redox.core.io import atomic_write_json
namespace=pathlib.Path(sys.argv[1]);placement=pathlib.Path(sys.argv[2]);gpu=sys.argv[3];review=sys.argv[4];claim=pathlib.Path(sys.argv[5])
allowed=os.sched_getaffinity(0)
with allocation_guard():host_capacity(1,minimum_ram_gib=2);os.sched_setaffinity(0,cpu_allocation(1))
probe(namespace);os.sched_setaffinity(0,allowed)
argv=[sys.executable,"-u","-m","redox.workflows.counterion_geometry","launch","--spec",str(namespace),"--gpu",gpu,"--out",str(placement),"--pilot-review",review]
with (claim/"controller.log").open("x") as log:child=subprocess.Popen(argv,stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
owner=dict(host=os.uname().nodename.split(".")[0],pid=child.pid,process_start=pathlib.Path(f"/proc/{child.pid}/stat").read_text().rsplit(")",1)[1].split()[19])
receipt=dict(status="controller_spawned_scientific_start_unverified",owner=owner,argv=argv,placement=str(placement))
atomic_write_json(claim/"try_001.json",receipt);print(json.dumps(receipt))
'''
(output/'remote_launch.py').write_text(remote)
record['source_sha256'][str(output/'remote_launch.py')] = sha256_file(output/'remote_launch.py')
for host, gpu in slots:
    if not pending:
        break
    row = pending[0]
    namespace = Path(row['namespace'])
    try:
        if sha256_file(namespace) != row['namespace_sha256']:
            raise ValueError('Prepared namespace changed')
        spec = verify(namespace)
        plan = json.loads((Path(spec['case'])/'plan.json').read_text())
        observed = capacity.capacity(host, gpu, str(ROOT.parent/'miniforge3/envs/redox/bin/python'))
        starting = sum(32*2**30 for p in record['launches'] if p['host'] == host
            and not (Path(p['namespace']).parent/'case/residency_guard/owner.json').exists())
        observed['promised'] += starting
        if not capacity.fits(plan, observed):
            continue
        space = storage.storage_headroom(rows, peers, plan['maximum_durable_output_bytes'], owners.alive, shutil.disk_usage(root).free)
        if not space['sufficient']:
            record['refusals'].append(dict(reason='aggregate_storage', **space))
            break
        claim = namespace.parent/'dispatch_claim'
        claim.mkdir()
        atomic_write_json(claim/'owner.json', dict(queue=str(output), host=os.uname().nodename,
            pid=os.getpid(), process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19]))
        pending.pop(0)
        argv = ['env', 'PYTHONPATH='+str(ROOT/'src'), 'PYTHONDONTWRITEBYTECODE=1',
            *[k+'=1' for k in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')],
            str(ROOT.parent/'miniforge3/envs/redox/bin/python'), '-c', remote, str(namespace),
            str(namespace.parent/'placement_001'), str(gpu), str(review), str(claim)]
        run = subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8',host,shlex.join(argv)],
            capture_output=True,text=True,check=True,timeout=180)
        record['launches'].append(dict(key=row['key'],namespace=str(namespace),host=host,GPU=gpu,
            capacity=observed,receipt=json.loads(run.stdout)))
    except Exception as error:
        record['refusals'].append(dict(key=row['key'],host=host,GPU=gpu,error=repr(error)))
        if (namespace.parent/'dispatch_claim').exists():
            break
    atomic_write_json(output/'dispatch.json', record)
record.update(finished=time.time(),unsubmitted=[r['key'] for r in pending])
atomic_write_json(output/'dispatch.json',record)
print(json.dumps(dict(launches=len(record['launches']),refusals=record['refusals'],unsubmitted=record['unsubmitted'])))
