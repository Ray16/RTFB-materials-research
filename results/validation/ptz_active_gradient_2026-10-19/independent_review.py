from pathlib import Path
import json,hashlib,math,time
import numpy as np
base=Path('results/validation/ptz_active_gradient_2026-10-19').resolve()
plan=json.loads((base/'plan.json').read_text());summary=json.loads((base/'summary.json').read_text());checks=[]
def check(label,passed):checks.append(dict(check=label,passed=bool(passed)))
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
for p,h in plan['producer_sha256'].items():check('producer hash '+p,sha(p)==h)
for task in plan['tasks']:
 state=task['state'];p=base/state/'report.json';r=json.loads(p.read_text());g=np.load(base/state/'gradient_Ha_Bohr.npy').ravel()
 for source,h in task['input_sha256'].items():check(state+' input '+source,sha(source)==h)
 check(state+' terminal success',r['status']=='active_gradient_and_finite_energy_checked')
 check(state+' gradient bytes',sha(base/state/'gradient_Ha_Bohr.npy')==r['gradient_sha256'])
 check(state+' finite gradient',np.isfinite(g).all())
 check(state+' summary bytes',sha(p)==summary['states'][state]['source_report_sha256'])
 check(state+' largest component',abs(g).argmax()==r['checked_column'])
 check(state+' analytic component',abs(g[r['checked_column']]-r['analytic_component_Ha_Bohr'])<1e-13)
 check(state+' maximum force',abs(abs(g).max()-r['max_gradient_Ha_Bohr'])<1e-13)
 em,ep=r['finite_energy_components_Ha'];fd=(ep-em)/(2*plan['FD_step_Bohr'])
 check(state+' independent FD quotient',abs(fd-r['finite_energy_derivative_Ha_Bohr'])<1e-13)
 error=abs(fd-g[r['checked_column']])
 check(state+' FD threshold',error<=plan['FD_absolute_acceptance_Ha_Bohr'])
 check(state+' FD error serialization',abs(error-r['finite_energy_absolute_error_Ha_Bohr'])<1e-13)
 repro=r['energy_Ha']*27.211386245988-task['expected_eV']
 check(state+' independently reproduced energy',abs(repro)<1e-3)
 check(state+' no canonical replacement',not r['production_admission'] and not r['geometry_optimized'] and not r['thermochemistry_recomputed'])
for p in base.glob('reserved_*/execution.json'):
 d=json.loads(p.read_text());check(str(p)+' execution success',d['status']=='process_completed' and d['returncode']==0)
terminal=json.loads((base/'terminal_verification.json').read_text());check('actual workers absent',len(terminal['processes'])==6 and all(v is None for v in terminal['processes'].values()))
q=json.loads((base/'quadratic_sensitivity.json').read_text());delta={r['state']:r['estimated_electronic_change_eV'] for r in q['states']}
check('charging-gap estimate identity',abs(delta['neu']+delta['ox2']-2*delta['ox1']-q['estimated_charging_gap_change_eV'])<1e-12)
for r in q['states']:
 for p,h in r['evidence_sha256'].items():check('quadratic evidence '+p,sha(p)==h)
check('estimate not a promoted correction',not q['production_admission'] and not q['energy_correction_applied'] and not q['actual_active_minimum_found'])
result=dict(status='independent_report_and_artifact_checks_passed' if all(c['passed'] for c in checks) else 'failed',checks=len(checks),failed_checks=sum(not c['passed'] for c in checks),per_case_quantum_recalculation=False,independent_FD_quotients_and_identity_checks=True,reviewed=time.time(),details=checks)
(base/'independent_review.json').write_text(json.dumps(result,indent=2)+'\n')
print({k:v for k,v in result.items() if k!='details'})
assert not result['failed_checks']
