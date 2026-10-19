"""Prepare two already-computed minima for independent physical review on CPU."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,time
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.core.protocol import load_record
from redox.qm.native_thermal import cache_path,NATIVE_THERMAL,load_native_thermal,select_policy
from redox.qm.hessian_memory import derivative_certificate
from redox.qm.projected_thermal import ensure_projected_thermal
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
a=argparse.ArgumentParser();a.add_argument('--name',choices=['minus','plus'],required=True);arg=a.parse_args()
base=Path(__file__).resolve().parent/arg.name;base.mkdir(exist_ok=False)
original=ROOT/'results/validation/ndi_methyl_radical_minimum_2026-10-19'
branch=original/arg.name;source=json.loads((original/'source_plus.json').read_text());branch_status=json.loads((branch/'status.json').read_text());assert branch_status['status']=='private_minimum_and_energy_verified_pending_certificate'
case=branch/'qualification/candidate';qual=json.loads((case/'status.json').read_text());assert qual['status']=='isolated_target_verified';q,m=source['charge'],source['mult'];canonical=Path(source['source_state_dir'])
for path,h in source['inputs_sha256'].items():assert sha(path)==h,path
proof={str(p):sha(p) for p in [case/'status.json',branch/'status.json',branch/'geometric_convergence.json',branch/'optimization_progress.json',original/'source_plus.json',branch/'opt.xyz',branch/'result.json']}
for path,h in qual['evidence_sha256'].items():assert sha(path)==h,path;proof[path]=h
producers=base/'producers';producers.mkdir();resolutions=[]
for path,h in source['producer_sha256'].items():
 p=Path(path)
 if sha(p)!=h:
  p=ROOT/'results/validation/spin_ensemble_integration_2026-10-19/common_before.py'
  assert Path(path).name=='common.py' and sha(p)==h,'Only the recorded analysis hook changed after source guards'
 assert p.is_file() and sha(p)==h,path
 dst=producers/(Path(path).stem+'_'+h[:12]+'.py');shutil.copy2(p,dst);proof[str(dst)]=h;resolutions.append(dict(declared=path,resolved=str(p),sha256=h))
atomic_write_json(base/'producer_resolution.json',resolutions);proof[str(base/'producer_resolution.json')]=sha(base/'producer_resolution.json')
prepared=base/'prepared_state';shutil.copytree(case/'work',prepared)
# Geometry is shared by the independent derivative and active-energy calculations.
assert sha(branch/'opt.xyz')==sha(prepared/'opt.xyz')
shutil.copytree(branch/'sp',prepared/'sp',dirs_exist_ok=True)
active=load_record(prepared,q,m);assert active
raw=json.loads((branch/'result.json').read_text());raw.update(opt_out=str(canonical/'opt.xyz'),minimum_repair_source=str(base/'superseded_canonical/opt.xyz'),minimum_repair_source_sha256=sha(canonical/'opt.xyz'),minimum_repair_review=str(base/'independent_review.json'),active_energy_source_record=branch_status['active_SP'],active_energy_source_sha256=branch_status['active_SP_sha256'])
atomic_write_json(prepared/'result.json',raw)
path,_=cache_path(prepared,q,m);native=json.loads(path.read_text());assert native['status']=='isolated_CDS_1e6_qualification_pending'
native.update(status='derivative_checked_solution_thermal',execution_derivatives=dict(cds_step_Bohr=1e-6,qualification_sha256=derivative_certificate(NATIVE_THERMAL,1e-6),qualification_only=False,scientific_model_changed=False,source_geometry_sha256=sha(prepared/'opt.xyz'),worst_contact_indices=qual['worst_contact_indices'],worst_contact_checks=qual['worst_contact_checks'],target_evidence_sha256={str(Path(k).relative_to(ROOT)):v for k,v in proof.items()}),production_admission_review=str(base/'independent_review.json'))
atomic_write_json(path,native)
kind,partition=ensure_projected_thermal(path,prepared/'opt.xyz');select_policy(prepared);assert load_native_thermal(prepared,q,m)
workers=json.loads((original/'launch.json').read_text())['workers'];owner=next(w for w in workers if w['id']=='ndi_methyl_'+arg.name)
receipt=dict(status='prepared_pending_independent_review',id=source['id'],state=source['state'],charge=q,mult=m,smiles=source['smiles'],canonical=str(canonical),prepared=str(prepared),native_record=str(path),qualified_case=str(case),branch=str(branch),active_SP=branch_status['active_SP'],canonical_inputs_sha256={str(canonical/n):sha(canonical/n) for n in ['opt.xyz','result.json']},evidence_sha256=proof,producer_controller=owner,minimum_frequency_cm=qual['minimum_frequency_cm'],native_gradient_max_Ha_Bohr=native['gradient_max'],g_thermal_total_eV=load_native_thermal(prepared,q,m)['g_thermal_total_eV'],source_optimization_claim='All five strict geomeTRIC criteria and final coordinate identity are independently rechecked.',canonical_modified=False,created=time.time())
atomic_write_json(base/'preparation.json',receipt)
# Retain the earlier independently implemented review/adoption algorithms, adding
# strict optimizer evidence checks in the separate reviewer before use.
for name in ['review.py','adopt.py']:shutil.copy2(ROOT/'results/validation/ndi_radical_adoptions_2026-10-19/ndi_methoxyethyl'/name,base/name)
print(arg.name,receipt['minimum_frequency_cm'],receipt['g_thermal_total_eV'])
