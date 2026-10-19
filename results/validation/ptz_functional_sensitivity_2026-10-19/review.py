"""Independent fixed-geometry electronic charging-gap comparison; never a fit."""
import hashlib,json,math,time
from pathlib import Path
base=Path(__file__).resolve().parent;plan=json.loads((base/'plan.json').read_text());checks=[];records={}
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def check(label,passed):checks.append(dict(check=label,passed=bool(passed)))
check('predeclared complete factorial',len(plan['tasks'])==9 and
      {(t['method']['id'],t['state']) for t in plan['tasks']}==
      {(m['id'],s) for m in plan['methods'] for s in ('neu','ox1','ox2')})
check('unique task keys',len({t['key'] for t in plan['tasks']})==9)
for p,h in plan['producer_sha256'].items():check('producer '+p,sha(p)==h)
relaxation=base.parent/'ptz_active_relaxation_2026-10-19'
for state in ('neu','ox1','ox2'):
 check(state+' frozen geometry matches completed relaxation',
       sha(base/f'input_{state}.xyz')==sha(relaxation/state/'relaxed.xyz'))
for task in plan['tasks']:
 path=base/'tasks'/task['key']/'status.json'
 if not path.exists() or json.loads(path.read_text())['status'] in ['starting','computing']:
  print('Pending point:',task['key']);raise SystemExit(2)
 report=json.loads(path.read_text());check(task['key']+' solver status',report['status']=='fixed_geometry_electronic_method_point_verified')
 if report['status']!='fixed_geometry_electronic_method_point_verified':continue
 check(task['key']+' executed declared task',report['task']==task)
 check(task['key']+' attempt after declaration',plan['created']<=report['started']<=report['finished'])
 for p,h in task['input_sha256'].items():check('geometry '+p,sha(p)==h)
 check(task['key']+' finite energy/units',math.isfinite(report['energy_Ha']) and abs(report['energy_eV']-report['energy_Ha']*27.211386245988)<1e-8)
 check(task['key']+' same stated spin model',abs(report['spin_diagnostics']['spin_contam'])<=.1 and report['SCF_converged'] is True)
 if task['method']['id']=='wb97mv':check(task['key']+' WB control',abs(report['WB_control_difference_eV'])<.001)
 records.setdefault(task['method']['id'],{})[task['state']]=report['energy_Ha']
methods=[]
for method in plan['methods']:
 e=records.get(method['id'],{});check(method['id']+' all three charges',set(e)=={'neu','ox1','ox2'})
 if len(e)==3:
  gap=math.fsum([e['neu'],e['ox2'],-2*e['ox1']])*27.211386245988
  methods.append(dict(method=method,electronic_gap_eV=gap,energies_Ha=e))
failed=sum(not c['passed'] for c in checks)
if len(methods)==3:
 control=next(r['electronic_gap_eV'] for r in methods if r['method']['id']=='wb97mv')
 for r in methods:r['difference_from_WB_eV']=r['electronic_gap_eV']-control
result=dict(status='fixed_geometry_method_sensitivity_independently_reviewed' if not failed and len(methods)==3 else 'failed',checks=len(checks),failed_checks=failed,details=checks,methods=methods,source_sha256={str(base/'plan.json'):sha(base/'plan.json'),str(base/'review.py'):sha(base/'review.py'),**{str(base/'tasks'/t['key']/'status.json'):sha(base/'tasks'/t['key']/'status.json') for t in plan['tasks']}},method_selected_for_production=False,parameters_fitted=False,canonical_values_changed=False,thermal_computed=False,experimental_free_energy_gap_is_not_electronic_validation=True,reviewed=time.time())
(base/'independent_review.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='details'}));assert not failed and len(methods)==3
