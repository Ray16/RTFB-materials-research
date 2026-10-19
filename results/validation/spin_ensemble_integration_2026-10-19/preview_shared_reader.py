"""Read-only CPU preview of the staged hook; no canonical ensemble policy exists."""
import importlib.util,json,hashlib,os,time
from pathlib import Path
from unittest.mock import patch
from redox.core import common
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
base=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('staged_equilibrium_common',base/'common_candidate.py');staged=importlib.util.module_from_spec(spec);spec.loader.exec_module(staged)
for key in ['ROOT','CONFIG_DIR','LIBRARY','CALCS','DFT','UMA','RESULTS']:setattr(staged,key,getattr(common,key))
spec=importlib.util.spec_from_file_location('independent_equilibrium_preview',ROOT/'scripts/validation/audit/recompute_axes.py');audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
spec=importlib.util.spec_from_file_location('collector_equilibrium_preview',ROOT/'scripts/pipeline/collect_verified_candidates.py');collector=importlib.util.module_from_spec(spec);spec.loader.exec_module(collector)
gid='search_pmdi_1864c2c74267';state='red2';sd=ROOT/'calcs/dft'/gid/state;selection=sd/'spin_ensemble_policy.json';assert not selection.exists()
policy_bytes=(base/'prepared_policy.json').read_bytes();old_exists,old_file,old_read,old_bytes=Path.exists,Path.is_file,Path.read_text,Path.read_bytes
checks=[]
def check(name,v):checks.append(dict(check=name,passed=bool(v)))
# Only this missing selection file is virtualized, inside this CPU process.
# Every component geometry, energy, raw Hessian and partition is read unchanged.
with patch.object(Path,'exists',lambda p:True if p==selection else old_exists(p)),patch.object(Path,'is_file',lambda p:True if p==selection else old_file(p)),patch.object(Path,'read_text',lambda p,*a,**kw:policy_bytes.decode() if p==selection else old_read(p,*a,**kw)),patch.object(Path,'read_bytes',lambda p:policy_bytes if p==selection else old_bytes(p)),patch.object(common,'apply_thermal_policy',staged.apply_thermal_policy):
 selected=common.read_result(gid,state);mixed_G,qc=common.free_energy(selected)
 raw=json.loads((sd/'result.json').read_text());g,ni,definition=audit.thermal_inputs(sd,raw)
 from redox.core.protocol import load_record
 independent_G=load_record(sd,-2,1)['e_smd_eV']+g
 check('staged_primary_applies_selection',selected['spin_ensemble_status']=='applied' and qc=='ok')
 check('primary_and_independent_G',abs(mixed_G-independent_G)<=1e-8 and ni==0)
 spin=collector.state_spin_evidence(gid,state)
 check('candidate_primary_gate',spin['spin_primary_admissible'] and spin['spin_confidence']=='DFT_native_equilibrium_ensemble')
 from redox.screening.scorecard import _min_spin_gap_eV
 score_gap,score_source=_min_spin_gap_eV(gid)
 check('scorecard_primary_gate',score_gap is not None and score_source.startswith('equilibrium_spin_ensemble:'))
 row=next(r for r in common.read_manifest() if r['id']==gid and r['state']==state)
 state_G,_,bound=collector.independent_state(gid,state,-2,collector.coverage.full_key(row['smiles']),audit)
 check('candidate_independent_state',state_G is not None and abs(state_G-mixed_G)<1e-8 and str(selection) in bound)
 independent_spin=audit.my_spin_evidence(gid)
 check('independent_spin_gate',independent_spin[2]!='DFT_check_required')
 red1=common.free_energy(common.read_result(gid,'red1'))[0]
 fc=common.free_energy(common.read_result('ferrocene','ox'))[0]-common.free_energy(common.read_result('ferrocene','neu'))[0]
 E2=red1-mixed_G-fc
 check('kinetics_remain_unvalidated',selected['single_surface_kinetics_validated'] is False)
check('no_production_policy_installed',not selection.exists())
report=dict(status='isolated_shared_reader_preview',checks=len(checks),failed_checks=sum(not c['passed'] for c in checks),preview_E2_vs_Fc_V=E2,preview_G_eV=mixed_G,independent_G_eV=independent_G,alternative_population=selected['spin_ensemble_alternative_population'],production_values_released=False,canonical_records_modified=False,source_frozen_common_unchanged=hashlib.sha256((ROOT/'src/redox/core/common.py').read_bytes()).hexdigest()==json.loads((base/'activation.json').read_text())['expected_source_sha256'],staged_common_sha256=hashlib.sha256((base/'common_candidate.py').read_bytes()).hexdigest(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),finished=time.time())
atomic_write_json(base/'shared_reader_preview.json',report);atomic_write_json(base/'shared_reader_preview_checks.json',checks);print(report)
assert report['failed_checks']==0 and report['source_frozen_common_unchanged']
