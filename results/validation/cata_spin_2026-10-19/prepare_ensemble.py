"""Prepare/review an explicit two-spin equilibrium; select only when writers are free."""
import hashlib,importlib.util,json,os,time
from pathlib import Path
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.validation.spin_ensemble import prepare_policy,evaluate_policy,select_policy
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
base=Path(__file__).resolve().parent;state=ROOT/'calcs/dft/hyp_cata_ndi_ethyl/red2';helper=ROOT/'calcs/spincheck/hyp_cata_ndi_ethyl/red2_m3'
policy_path=base/'prepared_policy.json'
if policy_path.exists():policy=json.loads(policy_path.read_text())
else:policy=prepare_policy(state,helper);atomic_write_json(policy_path,policy)
primary=evaluate_policy(state,policy)
spec=importlib.util.spec_from_file_location('independent_cata_spin_axes',ROOT/'scripts/validation/audit/recompute_axes.py');audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
spec=importlib.util.spec_from_file_location('independent_cata_ensemble',ROOT/'scripts/validation/audit/spin_ensemble.py');ind=importlib.util.module_from_spec(spec);spec.loader.exec_module(ind)
other=ind.read_ensemble(state,policy,root=ROOT,manifest=audit.MANIFEST,native_protocol=audit.NATIVE_THERMAL_CONTRACT,surface_reader=audit.native_single_thermal_inputs,partition_reader=audit.native_partition_inputs,spin_decision=audit.spin_decision_ok,strict=True)
checks=[dict(check=k,passed=abs(primary[k]-other[k])<1e-8,primary=primary[k],independent=other[k]) for k in ['free_energy_eV','correction_eV','alternative_population','gap_eV']]
checks += [dict(check='source '+p,passed=hashlib.sha256(Path(p).read_bytes()).hexdigest()==h) for p,h in policy['evidence_sha256'].items()]
assert all(c['passed'] for c in checks)
review=dict(status='prepared_two_spin_policy_independently_reviewed',checks=len(checks),failed_checks=0,values={k:primary[k] for k in ['free_energy_eV','correction_eV','alternative_population','gap_eV']},single_surface_kinetics_validated=False,production_policy_selected=False,details=checks,reviewed=time.time());atomic_write_json(base/'ensemble_review.json',review);print(json.dumps({k:v for k,v in review.items() if k!='details'}),flush=True)
try:
 selected=select_policy(state,policy);selected.update(selected=time.time(),review=str(base/'ensemble_review.json'),experimental_solution_rates_validated=False);atomic_write_json(base/'selection_receipt.json',selected);print('Selected source-bound equilibrium with both readers verified')
except BlockingIOError:
 atomic_write_json(base/'selection_receipt.json',dict(status='waiting_for_molecule_writers',canonical_modified=False,prepared_values=review['values'],reviewed=time.time()));print('Policy reviewed; active molecule writer prevents selection')
