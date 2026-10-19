#!/usr/bin/env python
"""Collect original39 and additional registered candidate results with independent native-input arithmetic.

Timestamped immutable CSV/XLSX snapshots and an atomic latest.json feed reports/slides.
--watch only reads calculations/queue markers; stops at23:55Chicago or a terminal queue.
No quantum jobs, parent-property substitutions or production-table writes.
"""
from __future__ import annotations
import argparse,csv,fcntl,hashlib,importlib.util,json,math,os,signal,shutil,socket,sys,time
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
sys.path.insert(0,str(ROOT/'scripts/mining'))
from rdkit import Chem
from redox.core.common import load_config,read_manifest,DFT,UMA,read_result,free_energy,apply_thermal_policy
from redox.core.protocol import load_record,record_path
from redox.build.contracts import validate_group
from redox.validation.provenance import referenced_files, calculation_gid
from redox.screening.comparisons import paired_potential_shifts
from redox.validation.spin_evidence import qualify_DFT_spin_gap, helper_spin_report_admission, native_partition_dependencies_ok
import collect_proposal_coverage as coverage


def audit_module():
 spec=importlib.util.spec_from_file_location('candidate_raw_independent_audit',ROOT/'scripts/validation/audit/recompute_axes.py')
 mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
 return mod


def sha(path):
 return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_atomic(path,value):
 tmp=path.with_name(path.name+f'.{os.getpid()}.tmp')
 tmp.write_text(json.dumps(value,indent=2)+'\n');tmp.replace(path)


def finite(value):
 return type(value) in (float,int) and math.isfinite(value)


def exploration_metadata(group,row):
 """Preserve search axes and route uncertainty in every export's common schema."""
 axes=group.get('search_axes') or {}
 row.update(search_core=axes.get('core'),search_anchor=axes.get('anchor'),
            search_cap=axes.get('cap'),search_ring=axes.get('ring'),
            search_contrast=group.get('comparison_contrast'),comparison_baseline=group.get('comparison_baseline'))
 if not row.get('chemistry_notes'):
  notes=[group.get('hypothesis'),group.get('provenance',{}).get('route_status')]
  row['chemistry_notes']=' '.join(n for n in notes if n) or None
 return row


def comparison_control_groups(rows,registered):
 """Resolve every declared control, including controls outside custom_candidates."""
 required={r['comparison_baseline'] for r in rows if r.get('comparison_baseline')}
 present={r.get('registry_id') for r in rows}
 definitions={g['id']:g for g in registered.values()}
 missing=required-set(definitions)
 if missing:raise ValueError('Unregistered comparison controls: '+', '.join(sorted(missing)))
 return [definitions[gid] for gid in sorted(required-present)]


def decision_input_paths(gids):
 """Include spin evidence whose mutation can withdraw an otherwise clean surface."""
 from redox.validation.provenance import spin_decision_paths
 return spin_decision_paths(ROOT,gids)


def validation_source_paths():
 paths={Path(__file__),ROOT/'scripts/mining/collect_proposal_coverage.py',
        ROOT/'scripts/validation/audit/recompute_axes.py',ROOT/'library/manifest.csv',
        ROOT/'results/discovery/parallel_expansion_2026-10-10/proposals.json',
        ROOT/'results/discovery/parallel_expansion_2026-10-10/sources.json'}
 paths.update((ROOT/'config').glob('*.py'));paths.update((ROOT/'config').glob('*.json'))
 paths.update((ROOT/'src/redox').rglob('*.py'))
 paths.update((ROOT/'scripts/validation/audit').glob('*.py'))
 for p in (ROOT/'config').glob('native*.json'):
  paths.update(referenced_files(coverage.safe_json(p) or {},ROOT))
 model_review=ROOT/'config/model_validation_limits.json'
 if model_review.is_file():paths.update(referenced_files(coverage.safe_json(model_review) or {},ROOT))
 return {p for p in paths if p.is_file()}


def independent_state(gid,st,q,key,audit):
 """Admit raw exact-graph state: addressed SP plus independent native thermal reader."""
 sd=DFT/gid/st;raw=coverage.safe_json(sd/'result.json')
 if not raw:return None,'pending raw geometry/result',{}
 if raw.get('charge')!=q:return None,'raw charge mismatch',{}
 rows=[r for r in read_manifest() if r['id']==gid and r['state']==st]
 if len(rows)!=1 or coverage.full_key(rows[0]['smiles'])!=key:return None,'manifest exact-structure mismatch',{}
 try:
  energy=load_record(sd,q,int(raw['mult']))
  correction,nimag,definition=audit.thermal_inputs(sd,raw)
  inputs={str(sd/'result.json'):sha(sd/'result.json'),str(sd/'opt.xyz'):sha(sd/'opt.xyz')}
  if energy:
   inputs[str(record_path(sd,q,raw['mult']))]=sha(record_path(sd,q,raw['mult']))
  if (not energy or definition!='native_solution_dft_qrrho_spin' or nimag!=0 or not finite(correction)
      or not audit.spin_decision_ok(gid,st,raw)):
   return None,'pending admitted active-SP/native thermal',inputs
  policy=sd/'thermal_policy.json';inputs[str(policy)]=sha(policy)
  identity=dict(geometry_sha256=sha(sd/'opt.xyz'),charge=q,mult=raw['mult'],protocol=audit.NATIVE_THERMAL_CONTRACT)
  address=hashlib.sha256(json.dumps(identity,sort_keys=True).encode()).hexdigest()[:16]
  thermal_path=sd/'native_thermal'/f'{address}.json';inputs[str(thermal_path)]=sha(thermal_path)
  for p in referenced_files(coverage.safe_json(thermal_path) or {},ROOT):inputs[str(p)]=sha(p)
  for name in ('hessian_cartesian.npy','hessian_diagnostics.json'):
   artifact=thermal_path.with_suffix('')/name
   if artifact.exists():inputs[str(artifact)]=sha(artifact)
  partition=audit.native_partition_inputs(sd,raw)
  if partition is None:return None,'corrected internal partition pending',inputs
  inputs.update(partition['evidence_sha256'])
  ensemble=audit.spin_ensemble_inputs(sd,raw)
  if ensemble is not None:inputs.update(ensemble['evidence_sha256'])
  total=energy['e_smd_eV']+correction
  return total,'independently admitted native free energy',inputs
 except (OSError,ValueError,KeyError,TypeError):
  return None,'current record unreadable or changing; pending',{}


def queue_progress(queue):
 data=coverage.safe_json(queue/'queue.json') or {};ids=data.get('ids',[])
 if data.get('jobs'):
  states={j['name']:coverage.safe_json(queue/'job_status'/f"{j['name']}.json") or {} for j in data['jobs']}
  terminal=all(s.get('status') in ('verified','failed_preserved','dependency_failed') for s in states.values())
  result={gid:'preserved earlier records; field-specific QC separate' for gid in ids}
  for job in data['jobs']:
   gid=job.get('gid');status=states[job['name']].get('status','queued')
   if gid in result and status!='verified':result[gid]=status+'; field-specific QC separate'
  return result,terminal
 result={gid:'pending' for gid in ids}
 for gid in ids:
  if (queue/'done'/f'{gid}.json').exists():result[gid]='worker done marker; record QC separate'
  elif (queue/'claims'/gid/'owner.json').exists():result[gid]='claimed/running'
  else:
   fail=coverage.safe_json(queue/'failures'/f'{gid}.json') or {}
   if fail.get('attempts',0)>=2:result[gid]='terminal failure; partial data preserved'
   elif fail:result[gid]='retry pending; partial data preserved'
 terminal=bool(ids) and all(v in ('worker done marker; record QC separate','terminal failure; partial data preserved') for v in result.values())
 return result,terminal


# Compatibility export; the reusable state reader is owned by validation.
from redox.validation.spin_evidence import state_spin_evidence


def apply_spin_gates(rows,states,couples):
 by_state={}
 for s in states:
  spin=state_spin_evidence(s['registry_id'],s['state']) if s['registry_id'] else dict(
     spin_confidence='spin_evidence_pending',spin_gap_eV=None,spin_gap_source='pending',spin_evidence_path=None,spin_primary_admissible=False)
  s.update(spin);s['native_surface_G_eV']=s['native_solution_G_eV']
  if not spin['spin_primary_admissible']:
   s['native_solution_G_eV']=None
   if s['native_surface_G_eV'] is not None:s['native_solution_status']='surface calculation verified; spin validation pending'
  by_state[s['proposal_id'],s['state']]=s
 for c in couples:
  o,r=c['event'].split('->');allowed=all(by_state[c['proposal_id'],st]['spin_primary_admissible'] for st in (o,r))
  c['spin_qc']='admissible under stated spin evidence' if allowed else 'DFT spin validation required/wrong ground state'
  c['E_native_surface_estimate_vs_Fc_V']=c['E_native_vs_Fc_V']
  if not allowed:c['E_native_vs_Fc_V']=None;c['potential_status']='pending DFT spin validation'
 for row in rows:
  ss=[s for s in states if s['proposal_id']==row['id']];bad=[s for s in ss if not s['spin_primary_admissible']]
  row['native_surface_calculation_states']=sum(s['native_surface_G_eV'] is not None for s in ss)
  row['native_solution_ready_states']=sum(s['native_solution_G_eV'] is not None for s in ss)
  row['spin_confidence']=('wrong_ground_state' if any(s['spin_confidence']=='wrong_ground_state' for s in bad) else
     ('spin_ensemble_required' if any(s['spin_confidence']=='spin_ensemble_required' for s in bad) else
      ('matching_native_spin_thermals_required' if any(s['spin_confidence']=='matching_native_spin_thermals_required' for s in bad) else
       ('DFT_check_required' if bad else 'admissible stated spin model'))))
  row['spin_state_details']='; '.join(s['state']+':'+s['spin_confidence'] for s in ss)
  # Match the declared direction; a clean first couple may be shown independently.
  declared=json.loads(row['states']);rest=next(q for st,q,m,ne in declared if ne==0);direction=1 if max(ne for st,q,m,ne in declared)>0 else -1
  cc={c['q_ox']:c for c in couples if c['proposal_id']==row['id']}
  ordered=[];q=rest
  while q+(0 if direction<0 else 1) in cc:
   ordered.append(cc[q+(0 if direction<0 else 1)]);q+=direction
  row['E1_vs_Fc_V']=ordered[0]['E_native_vs_Fc_V'] if ordered else None
  row['E2_vs_Fc_V']=ordered[1]['E_native_vs_Fc_V'] if len(ordered)>1 else None
  row['E1_surface_estimate_V']=ordered[0]['E_native_surface_estimate_vs_Fc_V'] if ordered else None
  row['E2_surface_estimate_V']=ordered[1]['E_native_surface_estimate_vs_Fc_V'] if len(ordered)>1 else None
  # The separation is reference-free and signed: G(end)+G(rest)-2G(middle)
  # for either direction of a two-step ladder. Do not use abs(E1-E2), which
  # would hide potential inversion, or require a complete Fc calculation.
  declared_charge={state:charge for state,charge,mult,ne in declared}
  by_charge={state.get('requested_charge',declared_charge.get(state['state'])):state for state in ss}
  values=[by_charge.get(rest+i*direction,{}).get('native_solution_G_eV') for i in range(3)]
  row['wave_spacing_V']=values[0]+values[2]-2*values[1] if all(finite(v) for v in values) else None
  if bad:
   row['solution_potentials_status']='partial/pending DFT spin validation'
   row['dG_disproportionation_kJ_mol']=None;row['wave_spacing_V']=None;row['accessible_capacity_mAh_g']=None
   row['disproportionation_status']='pending DFT spin validation' if len(ss)>2 else row['disproportionation_status']
   row['accessible_capacity_status']='pending DFT spin validation'
  if bad and row['native_surface_calculation_states']:
   row['result_stage']=('native surfaces complete; spin validation pending' if
      row['native_surface_calculation_states']==row['declared_states'] else
      'partial native surfaces; remaining states and spin validation pending')


def apply_thermal_publication_hold(rows, states, couples, reference, path=None):
 """Withdraw affected primary quantities while their thermal method is corrected.

 Stored geometry/electronic evidence and immutable historical snapshots remain
 available. Missing/unrecognized review metadata fails closed.
 """
 path=Path(path or ROOT/'config/thermal_publication_review.json')
 review=coverage.safe_json(path) or {}
 if review.get('status')=='corrected_protocol_independently_admitted':return False
 reason='pending rigid-mode-projected thermochemistry'
 for row in rows:
  for field in ('E1_vs_Fc_V','E2_vs_Fc_V','E1_surface_estimate_V','E2_surface_estimate_V',
                'wave_spacing_V','dG_disproportionation_kJ_mol','accessible_capacity_mAh_g'):
   if field in row:row[field]=None
  row.update(solution_potentials_status=reason,disproportionation_status=reason,
             accessible_capacity_status=reason,native_solution_ready_states=0,
             native_thermal_QC_states=0,result_stage='stored calculations; thermal projection review pending')
 for state in states:
  for field in ('native_solution_G_eV','native_surface_G_eV'):
   if field in state:state[field]=None
  state['native_solution_status']=reason
 for couple in couples:
  for field in ('E_native_vs_Fc_V','E_native_surface_estimate_vs_Fc_V'):
   if field in couple:couple[field]=None
  couple['potential_status']=reason
 reference.update(reference_native_Fc_absolute_V=None,thermal_publication_status=reason,
                  thermal_publication_review=str(path),thermal_publication_review_sha256=sha(path) if path.is_file() else None)
 return True


def apply_model_selection_gates(rows,review_path=None):
 """Retain calculated surfaces but withhold unqualified decision-level capacity."""
 from redox.validation.model_limits import selection_admission
 for row in rows:
  declared=json.loads(row['states'])
  role='catholyte' if any(state[3]>0 for state in declared) else 'anolyte'
  review=selection_admission(row['exact_grafted_smiles'],[state[1] for state in declared],role,
                             review_path=review_path)
  row.update(review)
  row['model_accessible_capacity_estimate_mAh_g']=row.get('accessible_capacity_mAh_g')
  if not review['model_primary_admissible']:
   row['accessible_capacity_mAh_g']=None
   row['accessible_capacity_status']='provisional model window/path estimate; '+review['model_validation_status']


def raw_audit_and_quarantine(rows,states,couples,audit,reference=None,frozen_inputs=None,input_inventory=None,
                            input_owners=None,global_inputs=None):
 evidence=dict(frozen_inputs or {});checks=[];manual={};bad_ids=set()
 owners={str(p):set(gids) for p,gids in (input_owners or {}).items()}
 for s in states:
  if not s['registry_id']:continue
  G,reason,files=independent_state(s['registry_id'],s['state'],s['requested_charge'],s['exact_grafted_inchikey'],audit)
  for path in files:owners.setdefault(path,set()).add(s['registry_id'])
  evidence.update(files);manual[s['registry_id'],s['state']]=G
  published=s.get('native_surface_G_eV',s['native_solution_G_eV'])
  if published is not None:
   ok=G is not None and abs(G-published)<1e-8
   checks.append(dict(id=s['proposal_id'],field='native_solution_G_eV',state=s['state'],passed=ok,production_value=published,independent_value=G,reason=reason))
   if not ok:bad_ids.add(s['proposal_id']);s['native_solution_G_eV']=None;s['native_solution_status']='quarantined independent admission mismatch'
 fc=[]
 for st,q in [('ox',1),('neu',0)]:
  mr=next((r for r in read_manifest() if r['id']=='ferrocene' and r['state']==st),None)
  G,reason,files=independent_state('ferrocene',st,q,coverage.full_key(mr['smiles']) if mr else '',audit)
  for path in files:owners.setdefault(path,set()).add('ferrocene')
  evidence.update(files);fc.append(G)
  if reference is not None and G is not None:
   for name in files:
    path=Path(name)
    if path.parent==DFT/'ferrocene'/st/'native_thermal' and path.suffix=='.json':
     native=coverage.safe_json(path) or {};freq=native.get('frequencies_cm',[])
     if freq and all(finite(v) for v in freq):
      reference.setdefault('reference_native_Fc_minimum_frequencies_cm',{})[st]=min(freq)
 ref=fc[0]-fc[1] if all(v is not None for v in fc) else None
 if reference is not None:
  reference['reference_native_Fc_absolute_V']=ref
 for c in couples:
  o,r=c['event'].split('->');go=manual.get((c['registry_id'],o));gr=manual.get((c['registry_id'],r))
  surface=c.get('E_native_surface_estimate_vs_Fc_V',c['E_native_vs_Fc_V'])
  if surface is not None:
   expected=-(gr-go)-ref if all(v is not None for v in (go,gr,ref)) else None
   ok=expected is not None and abs(expected-surface)<1e-8
   checks.append(dict(id=c['proposal_id'],field='E_native_vs_Fc_V',state=c['event'],passed=ok,
      production_value=surface,independent_value=expected,reason='raw admitted candidate/Fc arithmetic; spin qualification separate'))
   if not ok:bad_ids.add(c['proposal_id']);c['E_native_vs_Fc_V']=None;c['potential_status']='quarantined independent arithmetic mismatch'
  # GAS Nelsen remains a model-labelled diagnostic. Only fully QC-clean values
  # appear in the verified collector; flagged/missing legs stay blank.
  if c['gas_Nelsen_status']!='QC-clean gas diagnostic':c['gas_Nelsen_sum_eV']=None
  if c['gas_Nelsen_sum_eV'] is not None:
   try:
    ro=coverage.safe_json(DFT/c['registry_id']/o/'result.json');rr=coverage.safe_json(DFT/c['registry_id']/r/'result.json')
    eo=load_record(DFT/c['registry_id']/o,ro['charge'],ro['mult']);er=load_record(DFT/c['registry_id']/r,rr['charge'],rr['mult'])
    xo=audit.cross_path(c['registry_id'],o,r,ro['charge'],ro['mult']);xr=audit.cross_path(c['registry_id'],r,o,rr['charge'],rr['mult'])
    co=coverage.safe_json(xo);cr=coverage.safe_json(xr)
    admissible=bool(eo and er and co and cr and not audit.cross_valid(co,ro['charge'],ro['mult'],audit.geom_sha1(DFT/c['registry_id']/r/'opt.xyz'))
       and not audit.cross_valid(cr,rr['charge'],rr['mult'],audit.geom_sha1(DFT/c['registry_id']/o/'opt.xyz')))
    expected=(co['e_gas_eV']-eo['e_gas_eV'])+(cr['e_gas_eV']-er['e_gas_eV']) if admissible else None
    ok=expected is not None and abs(expected-c['gas_Nelsen_sum_eV'])<.000051
    checks.append(dict(id=c['proposal_id'],field='gas_Nelsen_sum_eV',state=c['event'],passed=ok,
       production_value=c['gas_Nelsen_sum_eV'],independent_value=expected,reason='raw addressed four-point GAS diagnostic; solution validity not inferred'))
    if ok:evidence.update({str(xo):sha(xo),str(xr):sha(xr)})
    else:c['gas_Nelsen_sum_eV']=None;c['gas_Nelsen_status']='quarantined independent raw four-point mismatch'
   except (OSError,TypeError,ValueError,KeyError):
    c['gas_Nelsen_sum_eV']=None;c['gas_Nelsen_status']='pending current cross-point records changing'
 for row in rows:
  cohort=[s for s in states if s['proposal_id']==row['id']]
  ladder=sorted(cohort,key=lambda s:-s['requested_charge'])
  if 'wave_spacing_V' in row:
   definition=json.loads(row['states']);origin=next(charge for state,charge,mult,ne in definition if ne==0)
   step=1 if max(ne for state,charge,mult,ne in definition)>0 else -1
   states_by_charge={s['requested_charge']:s for s in cohort};members=[states_by_charge.get(origin+i*step) for i in range(3)]
   values=[manual.get((row['registry_id'],s['state'])) if s else None for s in members]
   admitted=all(s is not None and s.get('spin_primary_admissible',True) for s in members)
   expected=values[0]+values[2]-2*values[1] if admitted and all(v is not None for v in values) else None
   published=row['wave_spacing_V'];ok=(published is None and expected is None) or (finite(published) and expected is not None and abs(published-expected)<1e-8)
   checks.append(dict(id=row['id'],field='wave_spacing_V',state='two-step declared path',passed=ok,
       production_value=published,independent_value=expected,reason='raw charge-ordered native G; signed and reference-free; completeness checked'))
   if not ok:bad_ids.add(row['id'])
  if row.get('dG_disproportionation_kJ_mol') is not None:
   values=[manual.get((row['registry_id'],s['state'])) for s in ladder]
   expected=(values[0]+values[2]-2*values[1])*audit.FARADAY/1000 if len(values)==3 and all(v is not None for v in values) else None
   ok=expected is not None and abs(expected-row['dG_disproportionation_kJ_mol'])<1e-7
   checks.append(dict(id=row['id'],field='dG_disproportionation_kJ_mol',state='three-state native ladder',
       passed=ok,production_value=row['dG_disproportionation_kJ_mol'],independent_value=expected,reason='raw native G(O)+G(R)-2G(I)'))
   if not ok:bad_ids.add(row['id'])
  if row.get('accessible_capacity_mAh_g') is not None:
   definition=json.loads(row['states']);rest=next(q for st,q,m,ne in definition if ne==0)
   direction=1 if max(ne for st,q,m,ne in definition)>0 else -1
   cc={c['q_ox']:c for c in couples if c['proposal_id']==row['id']}
   lo,hi=audit.ELE.WINDOW_V_VS_FC;divider=audit.ELE.ANOLYTE_CATHOLYTE_DIVIDER_V;q=rest;count=0
   while q+(0 if direction<0 else 1) in cc:
    c=cc[q+(0 if direction<0 else 1)];E=c['E_native_vs_Fc_V']
    if E is None or c['integrity_status']!='intact_bound' or not lo<=E<=hi or not (E<divider if direction<0 else E>=divider):break
    count+=1;q+=direction
   mass=sum(a.GetMass() for a in Chem.AddHs(Chem.MolFromSmiles(row['exact_grafted_smiles'])).GetAtoms())
   expected=count*audit.FARADAY/(3.6*mass);ok=abs(expected-row['accessible_capacity_mAh_g'])<.0001
   checks.append(dict(id=row['id'],field='accessible_capacity_mAh_g',state='contiguous native window/integrity path',
       passed=ok,production_value=row['accessible_capacity_mAh_g'],independent_value=expected,reason='independent path traversal plus atomic mass'))
   if not ok:bad_ids.add(row['id'])
  if row['id'] in bad_ids:
   for k in ('E1_vs_Fc_V','E2_vs_Fc_V','wave_spacing_V','dG_disproportionation_kJ_mol','accessible_capacity_mAh_g'):
    row[k]=None
   row.update(solution_potentials_status='quarantined independent raw-input mismatch',
      result_stage='quarantined independent raw-input mismatch',native_solution_ready_states=0,
      accessible_capacity_status='quarantined')
  mol=Chem.AddHs(Chem.MolFromSmiles(row['exact_grafted_smiles']));mw=sum(a.GetMass() for a in mol.GetAtoms())
  nominal=row['n_e_nominal']*audit.FARADAY/(3.6*mw)
  checks.append(dict(id=row['id'],field='capacity_nominal_mAh_g',state='resting graph',
      passed=abs(nominal-row['capacity_nominal_mAh_g'])<=.051,production_value=row['capacity_nominal_mAh_g'],
      independent_value=nominal,reason='independent atomic-mass sum; declared electrons only'))
 # A file changed during collection: quarantine its dependent scientific numbers
 # instead of treating an in-flight geometry repair as an arithmetic discrepancy.
 changing=[]
 for path,digest in evidence.items():
  try:
   if sha(path)!=digest:changing.append(path)
  except OSError:changing.append(path)
 if input_inventory is not None:
  before,read_paths=input_inventory
  changing.extend(str(p) for p in set(before)^set(read_paths()))
 changing=sorted(set(changing));scopes={};global_change=False;affected=set()
 global_inputs={str(p) for p in (global_inputs or [])}
 for path in changing:
  dependencies=set(owners.get(path,[]));gid=calculation_gid(path,ROOT)
  if gid is not None:dependencies.add(gid)
  scopes[path]=sorted(dependencies)
  if path in global_inputs or not dependencies or 'ferrocene' in dependencies:global_change=True
  affected.update(dependencies)
 if global_change:affected.update(r['registry_id'] for r in rows)
 if reference is not None:
  reference.update(global_input_changed=global_change,input_change_scope=scopes,
      changed_input_gids=sorted(gid for gid in affected if gid),
      quarantined_input_gids=sorted(affected & {r['registry_id'] for r in rows if r['registry_id']}))
 if changing:
  for row in rows:
   if row['registry_id'] not in affected:continue
   for k in ('E1_vs_Fc_V','E2_vs_Fc_V','E1_surface_estimate_V','E2_surface_estimate_V','wave_spacing_V','dG_disproportionation_kJ_mol','accessible_capacity_mAh_g'):row[k]=None
   row['solution_potentials_status']='pending input changed during snapshot'
   row['native_solution_ready_states']=0;row['native_thermal_QC_states']=0
   if 'native_surface_calculation_states' in row:row['native_surface_calculation_states']=0
   row['result_stage']='input changing; physical numbers pending'
  for s in states:
   if s['registry_id'] not in affected:continue
   for k in ('native_solution_G_eV','native_surface_G_eV','native_g_thermal_eV','active_smd_energy_eV'):s[k]=None
   s['native_solution_status']='pending input changed during snapshot'
  for c in couples:
   if c['registry_id'] not in affected:continue
   c['E_native_vs_Fc_V']=None;c['E_native_surface_estimate_vs_Fc_V']=None
   c['potential_status']='pending input changed during snapshot'
   c['gas_Nelsen_sum_eV']=None;c['gas_Nelsen_status']='pending input changed during snapshot'
  if not global_change:
   for path in list(evidence):
    dependencies=set(owners.get(path,[]));gid=calculation_gid(path,ROOT)
    if gid is not None:dependencies.add(gid)
    if dependencies and dependencies<=affected and path not in global_inputs:del evidence[path]
 return checks,evidence,changing


def snapshot(out,queue):
 source=ROOT/'results/discovery/parallel_expansion_2026-10-10';groups=json.loads((source/'proposals.json').read_text());assert len(groups)==39
 sources=json.loads((source/'sources.json').read_text());manifest=read_manifest();registered={}
 for name in coverage.CONFIGS:
  cfg=load_config(name)
  for g in cfg.GROUPS:registered[Chem.MolToInchiKey(validate_group(cfg.SCAFFOLD,g))]=g
 faraday=load_config('electrolyte').FARADAY;rows=[];states=[];couples=[]
 # Freeze code/configuration and spin decisions before reading any property.
 gids={g['id'] for g in groups}|{g['id'] for g in registered.values()}|{'ferrocene'}
 def decision_paths():return validation_source_paths()|decision_input_paths(gids)
 initial_source_paths=validation_source_paths()
 spin_by_gid={gid:decision_input_paths({gid}) for gid in gids}
 initial_spin_paths=set().union(*spin_by_gid.values())
 input_owners={}
 for gid,paths in spin_by_gid.items():
  for p in paths:input_owners.setdefault(str(p),set()).add(gid)
 initial_paths=initial_source_paths|initial_spin_paths;frozen={str(p):sha(p) for p in initial_paths}
 for g in groups:
  r,s,c=coverage.collect(g,manifest,registered,sources,faraday);exploration_metadata(g,r);rows.append(r);states+=s;couples+=c
 keys={r['exact_grafted_inchikey'] for r in rows};extras=[]
 for g in load_config('custom_candidates').GROUPS:
  key=Chem.MolToInchiKey(validate_group('Cc1ccc(C[*:1])cc1',g))
  if key not in keys:
   r,s,c=coverage.collect(g,manifest,registered,sources,faraday);exploration_metadata(g,r);extras.append(r);states+=s;couples+=c
 controls=[]
 for g in comparison_control_groups(rows+extras,registered):
  r,s,c=coverage.collect(g,manifest,registered,sources,faraday);exploration_metadata(g,r);controls.append(r);states+=s;couples+=c
 allrows=rows+extras+controls;apply_spin_gates(allrows,states,couples)
 audit=audit_module();reference={};checks,evidence,changing=raw_audit_and_quarantine(
  allrows,states,couples,audit,reference,frozen,(initial_paths,decision_paths),input_owners,initial_source_paths)
 retained_gids=gids-set(reference.get('changed_input_gids',[]))
 apply_thermal_publication_hold(allrows,states,couples,reference)
 apply_model_selection_gates(allrows)
 comparisons=paired_potential_shifts(allrows)
 progress,terminal=queue_progress(queue)
 for r in allrows:
  r['fleet_progress']=progress.get(r['registry_id'],'outside current queue')
  r['experimentally_validated_for_exact_graft']=False
 now=datetime.now(ZoneInfo('America/Chicago'));dest=out/'snapshots'/now.strftime('%Y%m%d_%H%M%S_%f');dest.mkdir(parents=True)
 coverage.write_csv(dest/'original39.csv',rows);coverage.write_csv(dest/'states.csv',states);coverage.write_csv(dest/'couples.csv',couples)
 if extras:coverage.write_csv(dest/'additional_candidates.csv',extras)
 if controls:coverage.write_csv(dest/'comparison_controls.csv',controls)
 if comparisons:coverage.write_csv(dest/'matched_comparisons.csv',comparisons)
 coverage.write_csv(dest/'independent_arithmetic_audit.csv',checks)
 notes=dict(snapshot_Chicago=now.isoformat(),phase='native MeCN free energies vs matched native Fc; gas Nelsen is separate diagnostic',
  capacity_nominal_mAh_g='Declared electrons and exact active-species mass; no counterion/backbone correction.',
  accessible_capacity_mAh_g='Native potential window and integrity-contiguous path bookkeeping, withheld where a scoped model limitation applies; kinetics assessed separately.',
  model_accessible_capacity_estimate_mAh_g='Original model window/path estimate retained for diagnosis; not qualified capacity when model_primary_admissible is false.',
  model_max_observed_error_V='Largest reviewed free-parent residual for the scoped model limitation; never a measurement or error bound for an exact graft.',
  gas_Nelsen_sum_eV='Only QC-clean current-protocol four-point GAS diagnostic; never a validated solution rate or full solution reorganization energy.',
  empty_numeric_cells='Pending, flagged or not independently admissible. Never substituted parent/zero values.',
  experimental_accuracy='No exact graft has an experimentally validated property in this dataset.',
  paired_shifts='Candidate minus declared exact-graft control, same charge/electron ladder and snapshot. Each event requires two admitted values; no surface/parent substitution or improvement claim.',
  spin_qc='UMA gaps<0.5eV require DFT alternatives; DFT electronic gaps within0.2eV require matching native spin thermals; native gaps within0.1eV require spin-ensemble evaluation. Raw surfaces remain diagnostics; far positive electronic gaps retain explicit thermal uncertainty.',
  native_g_thermal_eV='Corrected native partition of the selected single spin surface, before any equilibrium mixing correction.',
  effective_g_thermal_eV='Native term added to the selected active SMD energy; includes an admitted equilibrium spin correction when selected.',
  spin_ensemble_correction_eV='Equilibrium G minus selected-surface G. Component spin degeneracy is included once, not multiplied again.',
  spin_ensemble_alternative_population='Alternative-spin equilibrium population at 298.15 K under the selected two-spin model; this does not validate kinetics.',
  workbook='Immutable snapshot of independently audited external pipeline records, no spreadsheet formulas.')
 sheets={'Original39':rows,'States':states,'Couples':couples}
 if extras:sheets['Additional']=extras
 if controls:sheets['Controls']=controls
 if comparisons:sheets['Comparisons']=comparisons
 cells=coverage.workbook(dest/'verified_candidates.xlsx',sheets,notes)
 summary=dict(as_of_Chicago=now.isoformat(),original39_count=39,additional_candidate_count=len(extras),comparison_control_count=len(controls),
  paired_comparison_count=len(comparisons),paired_complete_potential_count=sum(r['status']=='both potential shifts admitted' for r in comparisons),
  registered_original39=sum(r['registered'] for r in rows),
  original39_native_potential_molecules=sum(r['solution_potentials_status']=='verified native calculation' for r in rows),
  original39_native_ready_states=sum(s['native_solution_status']=='verified calculation' for s in states if s['proposal_id'] in {r['id'] for r in rows}),
  original39_DFT_geometry_states=sum(s['dft_geometry_present'] for s in states if s['proposal_id'] in {r['id'] for r in rows}),
  arithmetic_checks=len(checks),arithmetic_discrepancies=sum(not r['passed'] for r in checks),
  workbook_cells_verified=cells,inputs_changed_during_snapshot=changing,queue_terminal=terminal,
  additional_candidates_with_model_selection_hold=sum(r.get('model_primary_admissible') is False for r in extras),
  spin_input_gids=sorted(retained_gids),
  spin_input_paths=sorted({str(p) for gid in retained_gids for p in spin_by_gid[gid]}),
  spin_input_paths_by_gid={gid:sorted(str(p) for p in spin_by_gid[gid]) for gid in retained_gids},
  phase_conventions=notes,**reference,source_sha256={str(p):frozen[str(p)] for p in initial_source_paths},
  admitted_input_sha256=evidence)
 (dest/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
 payload=dict(snapshot_path=str(dest),summary_path=str(dest/'summary.json'),workbook_path=str(dest/'verified_candidates.xlsx'),
     feed_schema=2,original39_csv=str(dest/'original39.csv'),additional_csv=str(dest/'additional_candidates.csv') if extras else None,
     as_of_Chicago=now.isoformat(),summary=summary,original39=rows,additional_candidates=extras,comparison_controls=controls,matched_comparisons=comparisons)
 # Slides consume only the latest independently checked snapshot. Partial states
 # remain visible as pending; an arithmetic discrepancy blocks publication.
 if summary['arithmetic_discrepancies']==0:
  save_atomic(out/'latest.json',payload)
  exports=[(dest/'verified_candidates.xlsx','LATEST_verified_candidates.xlsx'),
           (dest/'original39.csv','LATEST_original39.csv')]
  if extras:exports.append((dest/'additional_candidates.csv','LATEST_additional_candidates.csv'))
  if comparisons:exports.append((dest/'matched_comparisons.csv','LATEST_matched_comparisons.csv'))
  for source,name in exports:
   target=out/name;staged=out/(f'.{name}.{os.getpid()}.tmp')
   shutil.copy2(source,staged);staged.replace(target)
 return summary,dest


def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',required=True,type=Path)
 ap.add_argument('--queue',type=Path,default=ROOT/'results/validation/candidate_fleet_today_2026-10-11')
 ap.add_argument('--watch',action='store_true');ap.add_argument('--interval-seconds',type=int,default=120)
 ap.add_argument('--max-hours',type=float,help='optional watch duration spanning midnight')
 a=ap.parse_args();out=a.out.resolve();out.mkdir(parents=True,exist_ok=True)
 if a.interval_seconds<30:ap.error('interval must be at least30seconds')
 lock=(out/'collector.lock').open('a+');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
 from redox.workflows.publication import watch_deadline
 now=datetime.now(ZoneInfo('America/Chicago'))
 try:deadline=datetime.fromtimestamp(watch_deadline(a.max_hours,now=now),ZoneInfo('America/Chicago'))
 except ValueError as exc:ap.error(str(exc))
 status=dict(status='collecting',host=socket.gethostname(),pid=os.getpid(),
  process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
  affinity=sorted(os.sched_getaffinity(0)),threads=1,deadline_Chicago=deadline.isoformat(),snapshots=0,quantum_jobs_submitted=0)
 def stop(signum,frame):raise KeyboardInterrupt(f'signal{signum}')
 signal.signal(signal.SIGTERM,stop)
 try:
  while True:
   summary,dest=snapshot(out,a.queue.resolve());status.update(snapshots=status['snapshots']+1,last_snapshot=str(dest),
       last_summary=summary,status='watching' if a.watch else 'complete');save_atomic(out/'watch_status.json',status)
   print(json.dumps({k:summary[k] for k in ('as_of_Chicago','original39_native_potential_molecules','arithmetic_discrepancies','queue_terminal')}),flush=True)
   if not a.watch or summary['queue_terminal'] or datetime.now(ZoneInfo('America/Chicago'))>=deadline:break
   time.sleep(min(a.interval_seconds,max(0,(deadline-datetime.now(ZoneInfo('America/Chicago'))).total_seconds())))
  status['status']='complete';save_atomic(out/'watch_status.json',status)
 except (Exception,KeyboardInterrupt) as exc:
  status.update(status='stopped',error=repr(exc));save_atomic(out/'watch_status.json',status);raise
 finally:lock.close()
if __name__=='__main__':main()
