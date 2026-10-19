"""Independently read completed state records without modifying calculations."""
import argparse,hashlib,importlib.util,json,os,time
from pathlib import Path
from redox.core.common import ROOT,DFT,read_manifest,free_energy,apply_thermal_policy
from redox.core.protocol import load_record
from redox.core.io import atomic_write_json
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity

parser=argparse.ArgumentParser();parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
with allocation_guard():
 resources=host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
base=Path(__file__).resolve().parent;args.out.mkdir(parents=True,exist_ok=False)
path=ROOT/'scripts/pipeline/collect_verified_candidates.py'
spec=importlib.util.spec_from_file_location('independent_cata_collection',path);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
audit=module.audit_module();plan=json.loads((base/'plan.json').read_text());records=[];checks=[]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def check(label,passed):checks.append(dict(check=label,passed=bool(passed)))
for p,h in plan['producer_sha256'].items():check('frozen producer/control '+p,sha(p)==h)
for task in plan['tasks']:
 report_path=base/'tasks'/task['key']/'status.json'
 if not report_path.exists():continue
 report=json.loads(report_path.read_text())
 if report['status']!='state_records_verified':continue
 before={p:sha(p) for p in report['output_sha256']}
 for p,h in report['output_sha256'].items():check(task['key']+' producer output '+p,before[p]==h)
 gid,st=task['id'],task['state'];row=next(r for r in read_manifest() if r['id']==gid and r['state']==st)
 G,status,inputs=module.independent_state(gid,st,int(row['charge']),module.coverage.full_key(row['smiles']),audit)
 check(task['key']+' independent physical/native admission',G is not None and status=='independently admitted native free energy')
 sd=DFT/gid/st;raw=json.loads((sd/'result.json').read_text());energy=load_record(sd,int(raw['charge']),int(raw['mult']))
 primary,quality=free_energy(apply_thermal_policy(dict(raw,e_smd_eV=energy['e_smd_eV']),sd))
 check(task['key']+' primary and independent G agreement',primary is not None and G is not None and abs(primary-G)<1e-8)
 for p,h in before.items():check(task['key']+' input stability '+p,sha(p)==h)
 records.append(dict(id=gid,state=st,independent_G_eV=G,primary_G_eV=primary,quality=quality,independent_admission=status,input_sha256=inputs,exact_graft_experimentally_validated=False,preferred_spin_not_inferred_from_this_review=True))
result=dict(status='completed_state_records_independently_checked' if all(c['passed'] for c in checks) else 'failed',states_reviewed=len(records),checks=len(checks),failed_checks=sum(not c['passed'] for c in checks),records=records,resources=resources,CPU_affinity=sorted(os.sched_getaffinity(0)),reviewed=time.time(),scope='Only completed selected-state records; adjacent potentials and alternate-spin choice remain separate')
atomic_write_json(args.out/'summary.json',result);print(json.dumps({k:v for k,v in result.items() if k not in ['records','resources']}))
assert not result['failed_checks']
