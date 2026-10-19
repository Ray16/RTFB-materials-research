#!/usr/bin/env python
"""Read-only matched-parent benchmark with independent G and explicit spin gates.

    python scripts/validation/benchmark/review_matched_parents.py --only ID,ID --out NEW_DIR

No quantum calculations, fitted corrections or graft-property inheritance.
"""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import time

ROOT=Path(__file__).resolve().parents[3]


def load(path,name):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--only',required=True,help='comma-separated exact benchmark parent IDs')
    ap.add_argument('--out',type=Path,required=True)
    a=ap.parse_args()
    for name in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS'):
        os.environ[name]='2'
    from redox.core.common import load_config,read_manifest,read_result,free_energy,DFT
    from redox.core.io import atomic_write_json
    from redox.properties.potentials import fc_live
    from redox.validation.provenance import spin_decision_paths
    from redox.validation.spin_evidence import spin_comparison_running
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    with allocation_guard():
        resources=host_capacity(2);cpus=cpu_allocation(2);os.sched_setaffinity(0,cpus)
    selected=set(a.only.split(','))
    target_inputs={str(p):hashlib.sha256(p.read_bytes()).hexdigest()
                   for p in (ROOT/'data/raw/validation/two_wave_mecn').glob('*.csv')}
    entries=[b for b in load_config('benchmark_mecn').BENCHMARK if b['id'] in selected]
    if {b['id'] for b in entries}!=selected:raise ValueError('Unknown experimental parent ID')
    collector=load(ROOT/'scripts/pipeline/collect_verified_candidates.py','parent_benchmark_collector')
    audit=collector.audit_module()
    scorer=load(Path(__file__).with_name('score_benchmark_mecn.py'),'parent_benchmark_statistics')
    manifest={(r['id'],r['state']):r for r in read_manifest()}
    inputs=dict(target_inputs);checks=[];states={};rows=[]
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    def remember(path):
        path=Path(path)
        if path.is_file():inputs[str(path)]=sha(path)
    for path in collector.validation_source_paths()|spin_decision_paths(ROOT,selected|{'ferrocene'}):
        remember(path)
    remember(__file__);remember(Path(__file__).with_name('score_benchmark_mecn.py'))
    reference=fc_live();independent_reference=audit.fc_abs_live()
    checks.append(dict(check='independent_matched_Fc',ok=reference is not None
                       and independent_reference is not None and abs(reference-independent_reference)<1e-9))
    def read_state(gid,state,require_spin=True):
        key=(gid,state)
        if key in states:return states[key]
        row=manifest[key];q=int(row['charge']);molecule_key=collector.coverage.full_key(row['smiles'])
        total,quality,evidence=collector.independent_state(gid,state,q,molecule_key,audit)
        inputs.update(evidence)
        raw=read_result(gid,state);primary,pquality=free_energy(raw)
        matched=total is not None and primary is not None and pquality=='ok' and abs(total-primary)<1e-9
        if total is not None:checks.append(dict(check='independent_native_G',id=gid,state=state,ok=matched))
        spin=(collector.state_spin_evidence(gid,state) if require_spin else
              dict(spin_primary_admissible=True,spin_confidence='qualified_reference_model'))
        ready=matched and spin['spin_primary_admissible']
        result=dict(G_eV=total if ready else None,native_G_eV=total,
                    status='ready' if ready else quality+'; '+spin['spin_confidence'],
                    spin_confidence=spin['spin_confidence'],spin_gap_eV=spin.get('spin_gap_eV'),
                    thermal_definition=raw.get('thermal_definition') if raw else None)
        states[key]=result;return result
    for state in ('ox','neu'):read_state('ferrocene',state,False)
    for b in entries:
        event_rows={}
        busy=spin_comparison_running(b['id'])
        for e in b['events']:
            if e['kind']!='E':continue
            row=dict(id=b['id'],name=b['name'],family=b['bench_family'],event=e['event'],
                     kind='E',tier=e['tier'],exp_V=e['exp_V_vs_Fc'],calc_V=None,
                     calc_unrounded_V=None,error_V=None,doi=e['doi'],source=e['source'],where=e['where'],
                     status=e['tier_reason'],spin_qc=None,exact_graft_accuracy_validated=False)
            if e['tier']=='A' and not b.get('excluded') and not busy:
                ox,red=e['event'].split('->');O=read_state(b['id'],ox);R=read_state(b['id'],red)
                if O['G_eV'] is not None and R['G_eV'] is not None and checks[0]['ok']:
                    definitions={O['thermal_definition'],R['thermal_definition'],states['ferrocene','ox']['thermal_definition']}
                    if definitions=={'native_solution_dft_qrrho_spin'}:
                        value=O['G_eV']-R['G_eV']-independent_reference
                        primary=(free_energy(read_result(b['id'],ox))[0]-free_energy(read_result(b['id'],red))[0]-reference)
                        checks.append(dict(check='independent_potential',id=b['id'],event=e['event'],ok=abs(value-primary)<1e-9))
                        row.update(calc_V=round(value,3),calc_unrounded_V=value,
                                   error_V=round(value,3)-row['exp_V'],status='checked_native_parent',
                                   spin_qc=O['spin_confidence']+'; '+R['spin_confidence'])
                else:row['status']='pending endpoint/spin/reference qualification'
            elif b.get('excluded'):row['status']='excluded molecule: '+b['excluded']
            elif busy and e['tier']=='A':row['status']='pending active molecule spin writer'
            rows.append(row);event_rows[e['source'],e['event']]=row
        for e in b['events']:
            if e['kind']!='spacing':continue
            endpoints=[r for (source,event),r in event_rows.items() if source==e['source']]
            value=None
            if e['tier']=='A' and len(endpoints)==2 and all(r['calc_V'] is not None for r in endpoints):
                value=endpoints[0]['calc_V']-endpoints[1]['calc_V']
            rows.append(dict(id=b['id'],name=b['name'],family=b['bench_family'],event='dE12',kind='spacing',
                tier=e['tier'],exp_V=e['exp_V'],calc_V=value,calc_unrounded_V=None,
                error_V=value-e['exp_V'] if value is not None else None,doi=e['doi'],source=e['source'],where=e['where'],
                status='checked_native_parent' if value is not None else 'pending or source-excluded spacing',
                spin_qc=None,exact_graft_accuracy_validated=False))
    changed=[p for p,digest in inputs.items() if not Path(p).is_file() or sha(p)!=digest]
    failed=sum(not c['ok'] for c in checks)
    if changed or failed:
        for row in rows:
            row.update(calc_V=None,calc_unrounded_V=None,error_V=None,status='withdrawn: changed input or failed independent check')
    scored=[r for r in rows if r['kind']=='E' and r['tier']=='A' and r['calc_V'] is not None]
    requested=sum(e['kind']=='E' and e['tier']=='A' for b in entries if not b.get('excluded') for e in b['events'])
    summary=dict(status=('checked' if len(scored)==requested else 'partial_checked') if not changed and not failed else 'failed',
        scope='Selected exact free-parent calculations; not global tables or graft validation',selected_ids=sorted(selected),
        native_Fc_absolute_V=reference,arithmetic_checks=len(checks),arithmetic_discrepancies=failed,changed_inputs=changed,
        tier_A_potential_events=len(scored),tier_A_potential_molecules=len({r['id'] for r in scored}),
        requested_tier_A_potentials=requested,pending_tier_A_potentials=requested-len(scored),
        tier_A_original_studies=len({r['doi'] for r in scored}),tier_A_potential_stats=scorer.stats([r['error_V'] for r in scored]),
        potential_serialization_decimals=3,quantum_SCF_repeated=False,fitted_correction=False,
        exact_graft_accuracy_validated=False,scientific_accuracy_fully_validated=False,
        spin_model_limits='State-specific electronic/native/UMA confidence is retained in rows; electronic-only alternatives retain thermal uncertainty',
        resources=resources,cpu_affinity=cpus,completed=time.time())
    atomic_write_json(out/'input_inventory.json',inputs);atomic_write_json(out/'checks.json',checks)
    atomic_write_json(out/'states.json',{g+':'+s:r for (g,s),r in states.items()})
    with (out/'comparisons.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]) if rows else []);writer.writeheader();writer.writerows(rows)
    atomic_write_json(out/'summary.json',summary);print(json.dumps(summary,indent=2))
    if changed or failed:raise SystemExit(1)


if __name__=='__main__':main()
