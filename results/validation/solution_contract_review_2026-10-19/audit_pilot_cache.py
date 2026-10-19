"""Read-only audit of six archived PTZ SCF densities under the tightened reuse contract."""
import os,json,hashlib,importlib.util,subprocess
from pathlib import Path
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
from redox.validation.solution_relaxation import verified_density,scf_cache_matches,verified_implementation_controls
spec=importlib.util.spec_from_file_location('pair_cache_audit',ROOT/'scripts/validation/solution_reorg_pair.py');pair=importlib.util.module_from_spec(spec);spec.loader.exec_module(pair)
base=Path(__file__).resolve().parent;source=ROOT/'results/validation/new_candidate_solution_relaxation_2026-10-18';review=json.loads((source/'ptz_first_couple_audit.json').read_text());assert review['failed_checks']==0
checks=[];evidence={};records={};sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
resolved={};archive=base/'original_producers';archive.mkdir(exist_ok=True)
for path,h in review['evidence_sha256'].items():
 if sha(path)==h:resolved[path]=path;continue
 assert Path(path).suffix=='.py','Raw scientific evidence changed'
 saved=archive/(Path(path).stem+'_'+h[:12]+'.py')
 if not saved.exists():
  relative=str(Path(path).relative_to(ROOT))
  data=subprocess.check_output(['git','show','88880cc:'+relative],cwd=ROOT)
  assert hashlib.sha256(data).hexdigest()==h,'Original source cannot be resolved'
  saved.write_bytes(data)
 assert sha(saved)==h;resolved[path]=str(saved)
atomic_write_json(base/'original_evidence_resolution.json',resolved)
for path in sorted((source/'diagnostic').glob('scf_*.json')):
 record=json.loads(path.read_text());key=record['key'];state=key['state'];reference=key['reference_density_sha256'];geometry=next(p for p,h in key['identity']['input_sha256'].items() if h==key['geometry_sha256'])
 mol=pair.mol_at(state,geometry,pair.ACTIVE_SP['basis'])
 density=verified_density(record,path,evidence=review['evidence_sha256'],overlap=mol.intor('int1e_ovlp'),electrons=mol.nelectron,mult=state['mult'])
 matched=scf_cache_matches(record,state=state,geometry_sha256=key['geometry_sha256'],geometry_sha1=hashlib.sha1(Path(geometry).read_bytes()).hexdigest()[:16],method=key['identity']['method'],reference_density_sha256=reference,numerical_sources={str(p.relative_to(ROOT)):sha(p) for p in [ROOT/'src/redox/qm/nonequilibrium_smd.py',ROOT/'scripts/validation/physical_reorg.py']})
 checks.append(dict(file=str(path),passed=matched,density_shape=list(density.shape),electron_trace=record['electron_trace']));evidence[str(path)]=sha(path);evidence[str(path.with_suffix('.npz'))]=sha(path.with_suffix('.npz'))
 records[state['name'],key['geometry_sha256'],key['reference'] is not None]=record
specification=json.loads((source/'pair.json').read_text());directions=[]
for initial,final in [specification['states'],specification['states'][::-1]]:
 gi=sha(initial['geometry']);gf=sha(final['geometry'])
 ieq=records[initial['name'],gi,False];feq=records[final['name'],gf,False]
 cross=records[final['name'],gi,False];vertical=records[final['name'],gi,True]
 factor=27.211386245988
 values=dict(total_solution_relaxation_eV=(vertical['e_Ha']-feq['e_Ha'])*factor,
  equilibrium_geometry_component_eV=(cross['e_Ha']-feq['e_Ha'])*factor,
  frozen_slow_increment_eV=(vertical['e_Ha']-cross['e_Ha'])*factor,
  vertical_charge_change_eV=(vertical['e_Ha']-ieq['e_Ha'])*factor)
 original=next(d for d in json.loads((source/'diagnostic/result.json').read_text())['directional'] if d['initial']==initial['name'])
 for name,value in values.items():checks.append(dict(field=name,initial=initial['name'],passed=abs(original[name]-value)<=1e-8))
 checks.append(dict(field='same_geometry_CDS',initial=initial['name'],passed=abs(cross['cds_Ha']-vertical['cds_Ha'])<=1e-10))
 directions.append(dict(initial=initial['name'],final=final['name'],**values))
controls=json.loads((source/'pair.json').read_text())['implementation_controls'];implementation=verified_implementation_controls(controls,ROOT/'src/redox/qm/nonequilibrium_smd.py',backend='gpu')
result=dict(status='archived_PTZ_cache_revalidated_read_only',records=len(records),checks=len(checks),failed_checks=sum(not c['passed'] for c in checks),directional=directions,symmetric_mean_eV=sum(d['total_solution_relaxation_eV'] for d in directions)/2,implementation_reference_verified=implementation['implementation_reference_verified'],individual_quantum_SCF_not_repeated=True,legacy_reuse_requires_explicit_review=True,raw_files_modified=False,quantum_jobs_submitted=0,evidence_sha256=evidence,source_sha256=sha(Path(__file__)))
atomic_write_json(base/'pilot_cache_review.json',result);atomic_write_json(base/'pilot_cache_checks.json',checks);print({k:v for k,v in result.items() if k!='evidence_sha256'});assert len(records)==6 and not result['failed_checks']
