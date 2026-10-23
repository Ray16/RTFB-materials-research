#!/usr/bin/env python
"""Publish an all-axes decision ledger from exact identities and saved evidence.

This is a successor to molecular descriptor/route evaluation. It submits no jobs,
never replaces missing material properties with molecular proxies, and retains
all pending and excluded cases. Explicit future axis receipts can be supplied
only with matching graph/scope identities and unchanged backing files.
"""
import argparse,csv,hashlib,importlib.util,json,os,time
from collections import Counter
from pathlib import Path

from rdkit import Chem
from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.screening.expansion_review import attach_evaluations,identity
from redox.screening.feasibility import assess_candidate,feasible_ids,validate_policy,pareto_after_feasibility,finite
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def collect(policy_path,out,receipts_path=None,feed_path=None):
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    out.mkdir(parents=True,exist_ok=False)
    paths=dict(policy=policy_path,pool=ROOT/'results/discovery/systematic_families_2026-10-18/structures.json',
        feed=feed_path or ROOT/'results/candidates/current.json',routes=ROOT/'results/grafting/route_evidence.json',
        assessment=ROOT/'results/grafting/selection_assessment.json',bases=ROOT/'results/grafting/base_assisted/current.json',
        polymer=ROOT/'results/polymer_models/interaction_checked.json',
        spacing=ROOT/'results/reviews/2026-10-22/spacing_native_complete/validation.json')
    if receipts_path:paths['axis_receipts']=receipts_path
    blobs={k:p.read_bytes() for k,p in paths.items()};data={k:json.loads(v) for k,v in blobs.items()}
    snapshots=out/'source_snapshot';snapshots.mkdir()
    for key,blob in blobs.items():(snapshots/(key+'.json')).write_bytes(blob)
    source={str(path.resolve()):hashlib.sha256(blobs[key]).hexdigest() for key,path in paths.items()}
    for path in (Path(__file__),ROOT/'src/redox/screening/feasibility.py',ROOT/'src/redox/screening/expansion_review.py',
        ROOT/'src/redox/build/attachment_identity.py',ROOT/'scripts/pipeline/screen_grafting_routes.py'):
        source[str(path)]=sha(path)
        captured=snapshots/'producers'/path.relative_to(ROOT);captured.parent.mkdir(parents=True,exist_ok=True);captured.write_bytes(path.read_bytes())
    policy=data['policy'];validate_policy(policy)
    spec=importlib.util.spec_from_file_location('decision_route_gate',ROOT/'scripts/pipeline/screen_grafting_routes.py')
    gate=importlib.util.module_from_spec(spec);spec.loader.exec_module(gate)
    axes={a['id']:a for a in policy['required_axes']};rows=[];overrides={};candidates=[]
    if receipts_path:
        extra=data['axis_receipts']
        if extra['profile_id']!=policy['profile_id']:raise ValueError('Axis receipts belong to a different policy')
        for receipt in extra['records']:
            key=receipt['product_key']
            if key in overrides:raise ValueError('Duplicate explicit axis receipt')
            for axis in receipt['axes'].values():
                for item in axis.get('sources',[]):
                    path=Path(item['path']);path=path if path.is_absolute() else ROOT/path
                    if sha(path)!=item['sha256']:raise ValueError('Explicit axis backing evidence changed')
                    source[str(path.resolve())]=item['sha256']
            overrides[key]=receipt
    bases={r['id']:r for r in data['bases']['candidates']}
    for record in attach_evaluations(data['pool'],data['feed']):
        row=record['evaluation'];identifier=(row.get('registry_id') or row['id']) if row else record['pool_ids'][0]
        key=record['inchikey'];evidence={}
        def entry(axis,reason,status='missing',observations=None,source_key=None,available=False):
            value=dict(product_key=key,scope=axes[axis]['scope'],status=status,reason=reason,observations=observations or [],supporting_data_available=available)
            if source_key:value['sources']=[dict(path=str(paths[source_key].resolve()),sha256=source[str(paths[source_key].resolve())])]
            evidence[axis]=value;return value
        for axis in axes:entry(axis,'Required exact-case evidence has not been established')
        group=row or dict(id=identifier,frag=record['fragment'])
        route=gate.screen(group,data['routes'].get('candidates',{}),data['assessment'])
        reviewed=data['routes'].get('candidates',{}).get(identifier)
        if reviewed:
            from redox.build.attachment_identity import exact_attachment
            from redox.build.grafting import precursor
            pre,_=precursor(record['fragment']);states=record['states']
            states=json.loads(states) if isinstance(states,str) else states
            task=dict(identity=dict(attachment_fragment=record['fragment'],scaffold_smiles='Cc1ccc(C[*:1])cc1',
                product_smiles=record['grafted_smiles'],precursor_smiles=Chem.MolToSmiles(pre)),
                payload=dict(attachment_fragment=record['fragment'],product_smiles=record['grafted_smiles'],
                    precursor_smiles=Chem.MolToSmiles(pre),family=group.get('family',record['core']),states=states))
            binding=exact_attachment(task)['route_binding']
            if reviewed.get('route_binding')!=binding:raise ValueError('Route evidence is not bound to this exact new bond')
        if route['route_status']=='proceed' and not reviewed:raise ValueError('Route pass has no exact evidence')
        route_record=entry('attachment_route','Exact precursor/selectivity/conditions/core-retention review incomplete',
            status='qualified' if reviewed else 'partial',observations=[dict(route_status=route['route_status'],missing=route.get('missing_route_evidence',[]))],source_key='routes',available=bool(reviewed))
        route_record['value']=route['route_status']
        if row:
            observations=[dict(E1_vs_Fc_V=row.get('E1_vs_Fc_V'),E2_vs_Fc_V=row.get('E2_vs_Fc_V'),QC=row.get('solution_potentials_status'),
                model_primary_admissible=row.get('model_primary_admissible'))]
            entry('redox_window','Calculated capped-graft potentials do not supply qualified candidate-specific uncertainty bounds',
                'partial',observations,'feed',available=any(finite(row.get(k)) for k in ('E1_vs_Fc_V','E2_vs_Fc_V')))
            entry('capacity','Only molecular nominal/path estimates exist; material loading, counterion and backbone mass scenarios are not qualified',
                'partial',[dict(scope='capped_molecule_only',nominal_mAh_g=row.get('capacity_nominal_mAh_g'),molecular_path_mAh_g=row.get('accessible_capacity_mAh_g'))],'feed',available=finite(row.get('capacity_nominal_mAh_g')))
            entry('state_stability','State QC/disproportionation does not establish all required chemical-stability criteria',
                'partial',[dict(native_states=row.get('native_solution_ready_states'),disproportionation_kJ_mol=row.get('dG_disproportionation_kJ_mol'),
                    disproportionation_scope=row.get('disproportionation_status'))],'feed',available=bool(row.get('native_solution_ready_states')) or finite(row.get('dG_disproportionation_kJ_mol')))
            entry('charge_transfer','Qualified solution-state charge-transfer evidence and an application criterion are incomplete',
                'partial',[dict(lambda_i_eV=row.get('lambda_i_eV'),status=row.get('reorganization_status'))],'feed',available=finite(row.get('lambda_i_eV')))
        base=bases.get(identifier)
        if base:
            if identity(base['reaction']['product_smiles'])[0]!=record['grafted_smiles']:
                raise ValueError('Grafting terms belong to another exact product')
            for path,h in base['source_sha256'].items():
                if sha(path)!=h:raise ValueError('Grafting evidence changed')
                source[path]=h
            entry('grafting_thermodynamics','Base/ion representations and external capture-term checks do not establish full-route uncertainty or conformer convergence',
                'partial',[dict(reference_HCl_kJ_mol=base['reference_deltaG_kJ_mol'],base_models=base['base_models'],
                    uncertainty_interval=False)],'bases',available=finite(base['reference_deltaG_kJ_mol']))
        matching=[m for m in data['polymer']['models'] if m.get('candidate_id')==identifier and m.get('stage')=='stationary_DFT_electronic_structure' and not m.get('structural_review',{}).get('issues')]
        if matching:
            entry('motif_association','Existing full-segment/dense-oligomer electronic results are not the required motif-interaction assessment',
                'partial',[dict(checked_model_keys=[m['key'] for m in matching],scope='electronic_segment_or_oligomer_controls_only')],'polymer',available=True)
        if key in overrides:
            receipt=overrides.pop(key)
            if identity(receipt['grafted_smiles'])[0]!=record['grafted_smiles']:raise ValueError('Explicit receipt graph differs')
            evidence.update(receipt['axes'])
        states=json.loads(record['states']) if isinstance(record['states'],str) else record['states']
        pools=(["anolyte"] if any(s[3]<0 for s in states) else [])+(["catholyte"] if any(s[3]>0 for s in states) else [])
        candidate=dict(id=identifier,product_key=key,role=record['evaluation_role'] or 'candidate',axes=evidence,electrode_pools=pools)
        candidates.append(candidate)
        decision=assess_candidate(candidate,policy)
        rows.append(dict(**decision,electrode_pools=pools,grafted_smiles=record['grafted_smiles'],core=record['core'],
            pool_aliases=record['pool_ids'],role=candidate['role'],evidence=evidence,
            full_quantum_campaign_submitted=False))
    if overrides:raise ValueError('Explicit evidence refers to a product outside this pool')
    if any(sha(p)!=h for p,h in source.items()):raise ValueError('Decision source changed during collection')
    eligible=sorted(feasible_ids(rows,policy=policy));counts=Counter(r['status'] for r in rows)
    summary=dict(status='all_axes_decision_ledger_checked',profile_id=policy['profile_id'],
        input_candidate_as_of=data['feed'].get('as_of_Chicago'),source_sha256=source,
        policy_fully_numeric_targets_declared=all(a['kind']!='numeric' or a.get('minimum') is not None or a.get('maximum') is not None for a in axes.values()),
        unset_numeric_cutoffs=[a['id'] for a in axes.values() if a['kind']=='numeric' and a.get('minimum') is None and a.get('maximum') is None],
        exact_products=len(rows),decisions=dict(counts),eligible_ids=eligible,
        axis_counts={axis:dict(Counter(next(a['status'] for a in r['axes'] if a['axis']==axis) for r in rows)) for axis in axes},
        rule=policy['rule'],quantum_jobs_submitted=0,created=time.time(),
        interpretation='Pending denotes incomplete criteria/evidence, not demonstrated poor chemistry. No full decision front exists unless all required axes pass.')
    fronts={}
    for pool in ('anolyte','catholyte'):
        objectives=[dict(axis='redox_window',direction='minimize' if pool=='anolyte' else 'maximize'),
            dict(axis='capacity',direction='maximize'),dict(axis='charge_transfer',direction='minimize'),
            dict(axis='motif_association',direction='minimize_magnitude'),dict(axis='neighbor_coupling',direction='minimize_magnitude')]
        front=pareto_after_feasibility(candidates,policy,objectives,pool=pool)
        fronts[pool]={k:v for k,v in front.items() if k!='decisions'}
    atomic_write_json(out/'decision_fronts.json',fronts)
    atomic_write_json(out/'decisions.json',dict(profile_id=policy['profile_id'],candidates=rows))
    atomic_write_json(out/'summary.json',summary)
    with (out/'decisions.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=['id','core','status','eligible_for_Pareto',*axes])
        writer.writeheader()
        for row in rows:writer.writerow(dict(id=row['id'],core=row['core'],status=row['status'],eligible_for_Pareto=row['eligible_for_Pareto'],
            **{axis['axis']:axis['status'] for axis in row['axes']}))
    print(json.dumps({k:v for k,v in summary.items() if k not in ('source_sha256','axis_counts')}))
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--policy',type=Path,default=ROOT/'config/candidate_feasibility.json')
    parser.add_argument('--out',type=Path,required=True);parser.add_argument('--axis-receipts',type=Path);parser.add_argument('--feed',type=Path);args=parser.parse_args()
    collect(args.policy,args.out,args.axis_receipts,args.feed)


if __name__=='__main__':main()
