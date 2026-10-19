"""CPU-only independent equilibrium spin-population diagnostic; no primary admission."""
from pathlib import Path
import hashlib,json,math,time,importlib.util,os
from redox.core.common import ROOT
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
from redox.core.io import atomic_write_json
from redox.core.protocol import load_record
with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
spec=importlib.util.spec_from_file_location('independent_CN_partition',ROOT/'scripts/validation/audit/recompute_axes.py');audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
b=Path(__file__).resolve().parent;gid='search_pmdi_1864c2c74267';states=[];evidence={}
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
for path,m in [(ROOT/'calcs/dft'/gid/'red2',1),(ROOT/'calcs/spincheck'/gid/'red2_m3',3)]:
 raw=json.loads((path/'result.json').read_text());energy=load_record(path,-2,m);g,imag,definition=audit.thermal_inputs(path,raw)
 assert energy and g is not None and imag==0 and definition=='native_solution_dft_qrrho_spin'
 states.append(dict(mult=m,electronic_eV=energy['e_smd_eV'],thermal_eV=g,G_eV=energy['e_smd_eV']+g))
 for p in path.rglob('*'):
  if p.is_file() and p.suffix in ('.json','.xyz','.npy','.py'):evidence[str(p)]=sha(p)
kt=8.617333262145e-5*298.15;gap=states[1]['G_eV']-states[0]['G_eV'];weight=math.exp(-gap/kt);population=weight/(1+weight);correction=-kt*math.log1p(weight)
assert .01<population<.99
for p,h in evidence.items():assert sha(p)==h
report=dict(status='independent_equilibrium_spin_diagnostic',states=states,thermochemical_gap_eV=gap,alternative_population=population,ensemble_correction_from_selected_eV=correction,evidence_sha256=evidence,source_sha256=sha(Path(__file__)),primary_ensemble_integration_complete=False,kinetics_validated=False,production_records_modified=False,quantum_jobs_submitted=0,finished=time.time(),next_action='Integrate a source-bound two-spin equilibrium free energy into both primary readers and property arithmetic; retain separate kinetic admission.')
atomic_write_json(b/'summary.json',report);print({k:v for k,v in report.items() if k not in ['states','evidence_sha256']})
