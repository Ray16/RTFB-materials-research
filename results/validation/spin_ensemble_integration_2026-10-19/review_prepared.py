"""Compare production and independent equilibrium readers on preserved CN-PMDI data."""
import hashlib,importlib.util,json,os,time
from pathlib import Path
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
from redox.validation.spin_ensemble import evaluate_policy
base=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('independent_ensemble_axes',ROOT/'scripts/validation/audit/recompute_axes.py');audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
spec=importlib.util.spec_from_file_location('independent_ensemble_replay',ROOT/'scripts/validation/audit/spin_ensemble.py');ind=importlib.util.module_from_spec(spec);spec.loader.exec_module(ind)
policy=json.loads((base/'prepared_policy.json').read_text());source=ROOT/'calcs/dft/search_pmdi_1864c2c74267/red2'
production=evaluate_policy(source,policy)
review=ind.read_ensemble(source,policy,root=ROOT,manifest=audit.MANIFEST,native_protocol=audit.NATIVE_THERMAL_CONTRACT,surface_reader=audit.thermal_inputs,partition_reader=audit.native_partition_inputs,spin_decision=audit.spin_decision_ok,strict=True)
checks=[]
for key in ['free_energy_eV','correction_eV','alternative_population','gap_eV']:
 checks.append(dict(check=key,production=production[key],independent=review[key],passed=abs(production[key]-review[key])<=1e-8))
for path,h in policy['evidence_sha256'].items():checks.append(dict(check='source:'+path,passed=hashlib.sha256(Path(path).read_bytes()).hexdigest()==h))
report=dict(status='prepared_equilibrium_policy_independently_replayed',checks=len(checks),failed_checks=sum(not c['passed'] for c in checks),free_energy_eV=review['free_energy_eV'],mixing_correction_eV=review['correction_eV'],alternative_population=review['alternative_population'],single_surface_kinetics_validated=False,production_policy_selected=(source/'spin_ensemble_policy.json').exists(),canonical_records_modified=False,numerical_jobs_submitted=0,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),finished=time.time())
atomic_write_json(base/'independent_review.json',report);atomic_write_json(base/'independent_checks.json',checks)
print(report);assert report['failed_checks']==0 and not report['production_policy_selected']
