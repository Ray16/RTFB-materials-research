"""Freeze the complete declared 20-start, twelve-cell polymer sampling roster.

Only starting geometries and their identities are selected. No old electronic or
thermal values are imported into the corrected cohort. The original results
remain available as a separate historical collection.
"""
from collections import Counter
import json
from pathlib import Path

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.validation.archived_sources import ArchivedSources
from redox.validation.namespace_sources import NamespaceSources


def inventory(archive,census_path,energy_roster,reference_roots=()):
    reader=ArchivedSources(ROOT,archive,relative_to_prior=True)
    census=json.loads(census_path.read_text());expected={};choices={};sources={str(census_path):sha256_file(census_path),str(energy_roster):sha256_file(energy_roster)}
    for group in census['groups']:
        label='_'.join(map(str,group['occupied_sites'])) if group['occupied_sites'] else 'none'
        model_key=group['configuration']+'__occ_'+label
        path=archive/'results/polymer_models/campaigns/2026-10-25_diverse_conformer_inputs'/model_key/'plan.json'
        data=json.loads(path.read_text());sources[str(path)]=sha256_file(path)
        if len(data['tasks'])!=20 or data['experimental_occupancy_weights']:
            raise ValueError('Exactly twenty declared starts without invented occupancy weights required')
        for task in data['tasks']:
            if (task['key'] in expected or task['model_key']!=model_key or task['charge']!=0 or task['mult']!=1
                    or task['candidate_id']!='exp_phthalimide_control' or 'method_control' not in task['role']):
                raise ValueError('Distinct exact neutral method-control tasks required')
            expected[task['key']]=task
    if len(expected)!=240:raise ValueError('The complete twelve-by-twenty roster is required')
    def select(key,geometry,digest,kind,receipt,receipt_hash,priority):
        if key not in expected:return
        geometry=Path(geometry)
        if not geometry.is_file() or sha256_file(geometry)!=digest:
            raise ValueError('Available starting geometry must match its exact source hash: '+str(geometry))
        old=choices.get(key)
        row=dict(key=key,task=expected[key],geometry=str(geometry),geometry_sha256=digest,
                 source_kind=kind,source_receipt=str(receipt),source_receipt_sha256=receipt_hash)
        if old and priority==old['priority'] and old['geometry_sha256']!=digest:
            # Later physically replayed records deliberately represent the
            # documented continuation of that original input; keep both refs.
            if kind!='physically_replayed_native_geometry':
                raise ValueError('Ambiguous same-priority starting geometries for '+key)
            row['earlier_geometry']=dict(path=old['geometry'],sha256=old['geometry_sha256'],receipt=old['source_receipt'])
        if old is None or priority>=old['priority']:choices[key]={**row,'priority':priority}
    for ns in map(Path,json.loads(energy_roster.read_text())['specs']):
        spec=json.loads(ns.read_text());plan_path=Path(spec['attempt'])/'plan.json';plan=json.loads(plan_path.read_text());task=plan['tasks'][0]
        if task['key'] not in expected:raise ValueError('Electronic source outside the declared sampling roster')
        if task['smiles']!=expected[task['key']]['smiles'] or task['graph_sha256']!=expected[task['key']]['graph_sha256']:
            raise ValueError('Electronic task differs from the exact prospective graph')
        select(task['key'],Path(spec['attempt'])/'input.xyz',plan['scored_geometry_sha256'],
               'prepared_scored_geometry',ns,sha256_file(ns),1)
        sources[str(ns)]=sha256_file(ns);sources[str(plan_path)]=sha256_file(plan_path)
    for root in reference_roots:
        for plan_path in sorted(Path(root).glob('*/campaign/plan.json')):
            plan=json.loads(plan_path.read_text());task=plan['tasks'][0];key=task['key']
            receipt=plan_path.with_name('independent_energy_review.json')
            if key not in expected or not receipt.exists():continue
            report=json.loads(receipt.read_text())
            if (report['status']!='alternate_electronic_geometry_and_resource_checked_native_G_pending'
                    or report['portable_method_id']!=census['method_id']
                    or task['smiles']!=expected[key]['smiles'] or report['task']!=task):
                raise ValueError('Alternate-model reference task or review differs')
            select(key,plan_path.with_name('input.xyz'),plan['scored_geometry_sha256'],
                   'checked_alternate_electronic_reference_geometry',receipt,sha256_file(receipt),2)
            sources[str(receipt)]=sha256_file(receipt);sources[str(plan_path)]=sha256_file(plan_path)
    checked=archive/'results/reviews/2026-10-26/geometry_queue_resumption/independent_observer/checked'
    for path in sorted(checked.glob('*.json')):
        report=json.loads(path.read_text());key=Path(report.get('campaign','')).name
        if (key not in expected or report.get('status')!='alternate_smooth_geometry_and_resource_pilot_independently_checked'
                or report.get('portable_method_id')!=census['method_id']):continue
        original=Path(report['state_dir'])/'opt.xyz';digest=report['source_sha256'].get(str(original))
        if digest is None:raise ValueError('Geometry review must bind the actual optimized coordinates')
        geometry=reader.verify(original,digest)
        select(key,geometry,digest,'independently_checked_geometry',path,sha256_file(path),2)
        sources[str(path)]=sha256_file(path)
    for replay in census['source_replays']:
        path=Path(replay['path'])
        if sha256_file(path)!=replay['sha256']:raise ValueError('Historical census replay changed')
        data=json.loads(path.read_text());sources[str(path)]=sha256_file(path)
        bindings={r['review_sha256']:r for r in data.get('bindings',[])}
        for row in data['records']:
            key=row['key']
            if (key not in expected or row['graph_sha256']!=expected[key]['graph_sha256']
                    or row['method_id']!=census['method_id']):raise ValueError('Native source graph/model differs from prospective roster')
            review_path=reader.verify(row['review'],row['review_sha256']);report=json.loads(review_path.read_text())
            binding=bindings.get(row['review_sha256'])
            state_reader=reader if binding is None else NamespaceSources(json.loads(Path(binding['namespace']).read_text()))
            plan_path=state_reader.path(Path(report['campaign'])/'plan.json');plan=json.loads(plan_path.read_text())
            original=Path(plan['raw_root'])/'geometry.xyz';digest=report['source_sha256'].get(str(original))
            if digest is None:raise ValueError('Native review must bind its actual geometry')
            geometry=state_reader.verify(original,digest)
            select(key,geometry,digest,'physically_replayed_native_geometry',review_path,row['review_sha256'],3)
            sources[str(review_path)]=row['review_sha256'];sources[str(plan_path)]=sha256_file(plan_path)
    rows=[]
    for key in sorted(choices):
        row=choices[key];row.pop('priority');rows.append(row)
    missing=sorted(set(expected)-set(choices))
    return dict(status='complete_corrected_sampling_start_inventory' if not missing else 'sampling_start_inventory_missing_sources',
        declared_inputs=240,selected_inputs=len(rows),missing=missing,states=rows,
        source_selection='Original start IDs; replayed final geometry preferred, then checked geometry, then prepared scored geometry',
        historical_values_transferred=False,new_quantum_jobs_submitted=0,new_native_G_qualified=False,
        conformer_convergence_established=False,experimental_occupancy_weights=False,primary_candidate_admission=False,
        source_kind_counts=dict(Counter(r['source_kind'] for r in rows)),source_sha256=sources,
        producer_sha256=sha256_file(__file__))
