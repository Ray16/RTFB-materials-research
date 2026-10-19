"""Independent post-run arithmetic and saved-stencil review; no adoption."""
import hashlib,json,math,time
from pathlib import Path
base=Path(__file__).resolve().parent;plan=json.loads((base/'plan.json').read_text());reports={};checks=[]
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def check(label,passed):checks.append(dict(check=label,passed=bool(passed)))
for p,h in plan['producer_sha256'].items():check('producer '+p,sha(p)==h)
for task in plan['tasks']:
 path=base/task['state']/'status.json'
 if not path.exists() or json.loads(path.read_text())['status'] in ('starting','optimizing','checking_cold_final_gradient'):
  print('Pending state:',task['state']);raise SystemExit(2)
 r=json.loads(path.read_text());check(task['state']+' final status',r['status']=='active_level_stationary_geometry_verified')
 if r['status']!='active_level_stationary_geometry_verified':continue
 for p,h in task['input_sha256'].items():check('input '+p,sha(p)==h)
 for p,h in r['output_sha256'].items():check('output '+p,sha(p)==h)
 check(task['state']+' initial reproduction',abs(r['initial_reproduction_difference_eV'])<.001)
 check(task['state']+' optimizer/cold state',abs(r['cold_reproduction_difference_eV'])<.001 and r['geometry_optimization_converged'])
 check(task['state']+' stationary gradient',math.isfinite(r['final_max_gradient_Ha_Bohr']) and r['final_max_gradient_Ha_Bohr']<=plan['stationarity_max_Ha_Bohr'])
 import numpy as np
 g=np.load(base/task['state']/'final_gradient_Ha_Bohr.npy').ravel();minus,plus=r['finite_energy_Ha'];fd=(plus-minus)/(2*plan['FD_step_Bohr']);error=abs(fd-g[r['FD_column']])
 check(task['state']+' independent FD quotient',abs(fd-r['FD_derivative_Ha_Bohr'])<1e-12 and error<=plan['FD_absolute_acceptance_Ha_Bohr'])
 delta=(r['final_energy_Ha']-r['initial_energy_Ha'])*27.211386245988
 check(task['state']+' relaxation energy',abs(delta-r['electronic_relaxation_eV'])<1e-8 and delta<=.001)
 reports[task['state']]=r
failed=sum(not c['passed'] for c in checks)
result=dict(status='failed' if failed or len(reports)!=3 else 'electronic_geometry_sensitivity_independently_reviewed',checks=len(checks),failed_checks=failed,details=checks,production_admission=False,thermal_recomputed=False,Hessian_minima_proven=False,reviewed=time.time())
if len(reports)==3:
 HA=27.211386245988
 initial=(reports['neu']['initial_energy_Ha']+reports['ox2']['initial_energy_Ha']-2*reports['ox1']['initial_energy_Ha'])*HA
 relaxed=(reports['neu']['final_energy_Ha']+reports['ox2']['final_energy_Ha']-2*reports['ox1']['final_energy_Ha'])*HA
 result.update(initial_electronic_gap_eV=initial,relaxed_electronic_gap_eV=relaxed,actual_electronic_gap_change_eV=relaxed-initial,previous_local_estimate_eV=plan['previous_local_quadratic_gap_estimate_eV'],experimental_free_energy_gap_eV=plan['experimental_free_energy_gap_eV'],comparison_warning='Electronic gap is not a matched new free-energy prediction; old thermals must not be silently reused')
(base/'independent_review.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='details'}));assert not failed and len(reports)==3
