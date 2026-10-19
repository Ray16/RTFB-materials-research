"""Source-bound derived qualification; raw solution calculations stay unchanged."""
import hashlib,json,os,time
from pathlib import Path
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.qm.native_thermal import load_native_thermal,cache_path
from redox.properties.solution_reorg import summarize_relaxations
from redox.validation.solution_relaxation import verified_implementation_controls
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
base=Path(__file__).resolve().parent;association=json.loads((base/'review.json').read_text());assert association['both_selected_local_basins_associated'] and association['failed_checks']==0
pilot=ROOT/'results/validation/new_candidate_solution_relaxation_2026-10-18';raw_path=pilot/'diagnostic/result.json';raw=json.loads(raw_path.read_text());spec=json.loads((pilot/'pair.json').read_text());cache_path_review=ROOT/'results/validation/solution_contract_review_2026-10-19/pilot_cache_review.json';cache=json.loads(cache_path_review.read_text());assert cache['failed_checks']==0 and cache['records']==6
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
for path,h in {**association['evidence_sha256'],**cache['evidence_sha256']}.items():assert sha(path)==h,path
implementation=verified_implementation_controls(spec['implementation_controls'],ROOT/'src/redox/qm/nonequilibrium_smd.py',backend='gpu')
endpoints=[];evidence={str(base/'review.json'):sha(base/'review.json'),str(cache_path_review):sha(cache_path_review),str(raw_path):sha(raw_path)}
for s in spec['states']:
 sd=Path(s['state_dir']);assert sha(sd/'opt.xyz')==raw['identity']['input_sha256'][str(sd/'opt.xyz')]
 native=load_native_thermal(sd,s['charge'],s['mult']);assert native
 native_path=cache_path(sd,s['charge'],s['mult'])[0];post=Path(native['thermal_postprocessing_path'])
 endpoints.append(dict(state=s['name'],minimum_verified=True,conformer_policy='same_basin',conformer_basis='reciprocal strict optimization and proper-rotation chemical-symmetry match to selected local minimum',minimum_source=str(native_path),minimum_source_sha256=sha(native_path),active_SP_stationarity_verified=False))
 for path in [sd/'opt.xyz',native_path,post]:evidence[str(path)]=sha(path)
summary=summarize_relaxations(*cache['directional'],endpoint_evidence=endpoints,method_identity=raw['identity']['method'],implementation_validation=True)
assert summary['endpoint_model_verified'] and not summary['physical_flags']
summary.update(directional_descriptor_model_qualified=True,symmetric_approximation_diagnostic_passed=True,qualification_convention='composite electronic active-SP energies at selected geometry-level native minima',cached_arithmetic_and_densities_verified=True,independent_target_quantum_SCF_repeated=False,global_conformer_ensemble_proven=False,active_SP_stationary_endpoints_proven=False,finite_temperature_activation_validated=False,experimental_kinetic_accuracy_validated=False,polymer_constraint_proven=False)
for path,h in evidence.items():assert sha(path)==h
result=dict(status='local_pair_composite_electronic_descriptor_qualified',summary=summary,endpoint_evidence=endpoints,implementation_evidence=implementation,evidence_sha256=evidence,raw_outputs_modified=False,production_promotion=False,source_sha256=sha(Path(__file__)),qualified=time.time())
atomic_write_json(base/'descriptor_qualification.json',result);print(summary)
