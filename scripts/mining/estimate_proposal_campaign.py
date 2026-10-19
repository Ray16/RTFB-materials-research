#!/usr/bin/env python
"""Estimate candidate fleet costs from observed UMA intervals and native Hessian walls.

Forecasts are scheduling estimates. No predicted redox/kinetics values are generated.
"""
import argparse,csv,json,statistics,sys
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/'src'))
from rdkit import Chem
from rdkit.Chem import Descriptors,rdMolDescriptors
from redox.core.common import load_config
from redox.build.contracts import validate_group
from redox.qm.uma import candidate_mults
from redox.build.structures import n_conformers


def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',required=True,type=Path);a=ap.parse_args()
 if a.out.exists():ap.error('choose a fresh estimate directory')
 cfg=load_config('custom_candidates');groups=cfg.GROUPS
 originals={g['id'] for g in json.loads((ROOT/'results/discovery/parallel_expansion_2026-10-10/proposals.json').read_text())}
 # Completion intervals from one continuously reserved UMA pilot, excluding its
 # first state (whose elapsed launch time included waiting and model startup).
 specs=[('exp_viologen_methoxyethyl','ox2','ox1'),('exp_viologen_methoxyethyl','ox1','neu'),
        ('exp_viologen_methoxyethyl','neu','exp_phthalimide_control/neu'),('exp_phthalimide_control','neu','red1')]
 observations=[]
 for gid,first,last in specs:
  p1=ROOT/'calcs/uma'/gid/first/'result.json'
  p2=(ROOT/'calcs/uma'/last/'result.json') if '/' in last else ROOT/'calcs/uma'/gid/last/'result.json'
  r=json.loads(p2.read_text());m=Chem.AddHs(Chem.MolFromSmiles(r['smiles']))
  dt=p2.stat().st_mtime-p1.stat().st_mtime;scans=r['n_conf']*len(r['mults_scanned'])
  observations.append(dict(source=str(p2),seconds_completion_interval=dt,atoms=m.GetNumAtoms(),relaxations=scans,
        seconds_per_relaxation=dt/scans,scaled_coefficient_s_per_atom3=dt/(scans*m.GetNumAtoms()**3),
        caveat='File-completion interval is an approximate warm-start wall; contains small IO/write overhead. No wait-inclusive first-state interval used.'))
 uma_coef=statistics.median(o['scaled_coefficient_s_per_atom3'] for o in observations)
 thermal_sources=[]
 for gid,st in [('aq_benzylamino','neu'),('bm_1_4_benzoquinone','red1')]:
  for p in (ROOT/'calcs/dft'/gid/st/'native_thermal').glob('*.json'):
   r=json.loads(p.read_text())
   if r.get('wall_s') is not None:
    atoms=int((p.parent.parent/'opt.xyz').read_text().splitlines()[0])
    thermal_sources.append(dict(source=str(p),wall_s=r['wall_s'],atoms=atoms))
 thermal_sources.sort(key=lambda x:x['atoms']);small,large=thermal_sources[0],thermal_sources[-1]
 slope=(large['wall_s']-small['wall_s'])/(large['atoms']**3-small['atoms']**3)
 intercept=small['wall_s']-slope*small['atoms']**3
 rows=[]
 for g in groups:
  m=validate_group(cfg.SCAFFOLD,g);mh=Chem.AddHs(m);n=mh.GetNumAtoms();nuclei=sum(a.GetAtomicNum() for a in mh.GetAtoms())
  confs=list((ROOT/'library'/g['id']/'conformers').glob('conf_*.xyz'))
  nseed=len(confs) if confs else min(10,n_conformers(m))
  scans=sum(len(candidate_mults(nuclei-q,mult)) for st,q,mult,ne in g['states'])*nseed
  u=uma_coef*n**3*scans;t=max(0,intercept+slope*n**3)*len(g['states'])
  rows.append(dict(id=g['id'],original39=g['id'] in originals,family=g['family'],formula=rdMolDescriptors.CalcMolFormula(m),
    inchikey=Chem.MolToInchiKey(m),total_atoms=n,heavy_atoms=m.GetNumHeavyAtoms(),rotatable_bonds=Descriptors.NumRotatableBonds(m),
    declared_states=len(g['states']),canonical_conformers=len(confs),uma_relaxation_count_estimate=scans,
    solution_DFT_seed_optimizations_nominal=3*len(g['states']),
    uma_worker_hours_estimate=round(u/3600,3),native_hessian_worker_hours_estimate=round(t/3600,3),
    covered_stage_lower_bound_worker_hours=round((u+t)/3600,3),
    size_order_proxy=round(len(g['states'])*n**3/1e6,4),
    uncalibrated_stages='3-seed DFT minima; active SP; alternative-spin/SCF audits; solution reorganization; possible minimum repairs',
    scientific_properties_predicted=False))
 rows.sort(key=lambda r:(r['size_order_proxy'],r['rotatable_bonds'],r['id']))
 for order,r in enumerate(rows,1):r['cost_order']=order
 a.out.mkdir(parents=True)
 with (a.out/'all40_cost_order.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 (a.out/'all40_cost_order_ids.txt').write_text('\n'.join(r['id'] for r in rows)+'\n')
 now=datetime.now(ZoneInfo('America/Chicago'));deadline=now.replace(hour=23,minute=59,second=59,microsecond=0)
 remaining=(deadline-now).total_seconds()/3600
 aggregates={}
 for label,cohort in [('original39',[r for r in rows if r['original39']]),('first4',[r for r in rows if r['id'] in ('exp_phthalimide_control','exp_viologen_methoxyethyl','exp_aq_1oxy','exp_ndi_hexyl_pentyloxy')]),('all40',rows)]:
  hours=sum(r['covered_stage_lower_bound_worker_hours'] for r in cohort)
  active_workers=min(12,len(cohort));longest=max(r['covered_stage_lower_bound_worker_hours'] for r in cohort)
  aggregates[label]=dict(molecules=len(cohort),states=sum(r['declared_states'] for r in cohort),
     uma_plus_native_hessian_worker_hours_estimate=round(hours,2),
     usable_independent_molecule_workers=active_workers,
     ideal_parallel_hours_for_only_calibrated_stages=round(max(hours/active_workers,longest),2),
     sum_atom3_proxy=sum(r['size_order_proxy'] for r in cohort))
 summary=dict(as_of_Chicago=now.isoformat(),deadline_Chicago=deadline.isoformat(),hours_remaining=round(remaining,2),
     planned_parallel_workers=12,uma_observations=observations,uma_atom3_coefficient=uma_coef,
     native_hessian_observations=thermal_sources,thermal_atom3_fit=dict(intercept_seconds=intercept,slope=slope),
     aggregates=aggregates,
     resource_envelope=dict(cpu_threads_at12workers=24,V100_32GiB_workers_VRAM_not_pooled=True,
                            observed_dependencies='Root preflight passes on lambda1/2/4/11; lambda13 driver/NVML mismatch excluded'),
     limitations=['No complete new-candidate end-to-end timing exists; whole39 completion today cannot be guaranteed.',
       'Timing anchors are from lambda1; other GPUs/platforms and radicals may be slower.',
       'Fitted thermal curve has only two organic examples and excludes metal Fc. Extrapolation is coarse.',
       'Estimates intentionally omit uncalibrated DFT optimization/SP/reorganization/minimum-repair costs; they are lower bounds, not completion forecasts.',
       'Descriptor coverage today is complete; quantum/property coverage must be refreshed from actual records.'])
 (a.out/'throughput_estimate.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:summary[k] for k in ('hours_remaining','aggregates','limitations')},indent=2))
if __name__=='__main__':main()
