#!/usr/bin/env python
"""Freeze and independently check the small evidence set used in the discussion deck.

Read-only with respect to calculation records, registries and candidate decisions.
No quantum jobs are submitted and no missing property is imputed.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import time

from rdkit import Chem

ROOT = Path(__file__).resolve().parents[2]
EPOCH = ROOT / 'manuscript/evidence/2026-10-26_native_replay'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def canonical(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError('Invalid molecular graph')
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def metrics(errors):
    return dict(n=len(errors), MAE_V=math.fsum(map(abs, errors))/len(errors),
                RMSE_V=math.sqrt(math.fsum(x*x for x in errors)/len(errors)),
                bias_V=math.fsum(errors)/len(errors))


def collect(out):
    if out.exists():
        raise FileExistsError('Preserve the existing presentation snapshot: '+str(out))
    sources = {}

    def read(path):
        p = ROOT / path
        sources[str(p.relative_to(ROOT))] = digest(p)
        return json.loads(p.read_text())

    pool = read('results/discovery/systematic_families_2026-10-18/structures.json')
    enumeration = read('results/discovery/systematic_families_2026-10-18/summary.json')
    feed = read('results/candidates/current.json')
    policy_pointer = read('results/selection/current.json')
    policy = read(policy_pointer['summary'])
    plan = read('config/targeted_expansion_2026-10-23.json')
    benchmark = read(str((EPOCH/'benchmark.json').relative_to(ROOT)))
    baseline = read('manuscript/evidence/2026-10-21/benchmark.json')
    baseline_polymer = read('manuscript/evidence/2026-10-21/polymer_validation.json')
    bias = read('results/validation/reference_bias_2026-10-22/summary.json')
    basis = read('results/reviews/2026-10-22/neutral_basis_complete/summary.json')
    base = read('results/grafting/base_assisted/current.json')
    neighbor = read('results/reviews/2026-10-25/pmdi_neighbor_cycle_checked/summary.json')
    archive = read('results/releases/neighbor_cycle/current.json')
    pareto = read('results/figures/selection/descriptor_pareto.json')
    # This is the exact frozen feed underlying the descriptor figure, not a live overlay.
    pareto_feed_path = next(p for p in pareto['source_sha256'] if p.endswith('/feed.json'))
    pareto_feed = read(str(Path(pareto_feed_path).relative_to(ROOT)))
    nitro = read('results/validation/uma_nitro_2026-10-23/summary.json')
    v100 = read('results/validation/uma_v100_2026-10-24/summary.json')
    families = Counter(r['core'] for r in pool)
    recipes = sum(len(r['recipes']) for r in pool)
    pool_graphs = {canonical(r['grafted_smiles']) for r in pool}
    assert len(pool_graphs) == len(pool) == enumeration['exact_grafts'] == 837
    assert recipes == enumeration['recipes'] == 1593
    assert dict(families) == {k:v['total'] for k,v in enumeration['by_family'].items()}
    candidates = feed['original39'] + feed['additional_candidates']
    by_id = {r['id']:r for r in candidates + feed['comparison_controls']}
    candidate_graphs = {canonical(r['exact_grafted_smiles']) for r in candidates}
    assert len(candidates) == 75
    control_graphs = {canonical(r['exact_grafted_smiles']) for r in feed['comparison_controls']}
    assert len(candidate_graphs | pool_graphs | control_graphs) == policy['exact_products'] == 879
    assert len((candidate_graphs | pool_graphs) - control_graphs) == 878
    assert policy['decisions'] == {'pending':878,'comparison_control':1}
    matched = []
    for row in feed['matched_comparisons']:
        a, b = by_id[row['id']], by_id[row['baseline']]
        for wave in (1,2):
            field = f'E{wave}_vs_Fc_V'
            delta = a[field] - b[field]
            assert math.isclose(delta, row[f'delta_E{wave}_V'], abs_tol=1e-11)
        matched.append({**row, 'product_smiles':a['exact_grafted_smiles'],
                        'precursor_smiles':a['precursor_smiles']})
    assert len(matched) == 26
    systematic = [r for r in matched if r['core'] in families]
    assert len(systematic) == 18
    endpoints = benchmark['outliers']
    errors = [r['calc_V']-r['exp_V'] for r in endpoints]
    checked_metrics = metrics(errors)
    for k,v in checked_metrics.items():
        assert math.isclose(v, benchmark['metrics'][k], abs_tol=1e-12)
    cv = []
    for row in bias['cross_validation']:
        if row['holdout'] != 'id':
            continue
        residuals = []
        for endpoint in endpoints:
            training = [r for r in endpoints if r['id'] != endpoint['id']
                        and (row['model'] == 'constant_offset' or r['wave'] == endpoint['wave'])]
            offset = math.fsum(r['calc_V']-r['exp_V'] for r in training)/len(training)
            residuals.append(endpoint['calc_V']-endpoint['exp_V']-offset)
        checked_cv = metrics(residuals)
        for k,v in checked_cv.items():
            assert math.isclose(v,row['metrics'][k],abs_tol=1e-12)
        cv.append(dict(holdout='molecule', model=row['model'], metrics=checked_cv))
    assert plan['exact_new_molecules_enumerated'] == 0
    assert len(plan['slots']) == plan['design_slots'] == 34
    cutoff_path = ROOT/'results/figures/validation/thermal_cutoff_sensitivity_v2_benchmark.csv'
    sources[str(cutoff_path.relative_to(ROOT))] = digest(cutoff_path)
    cutoff = list(csv.DictReader(cutoff_path.open()))
    assert len(cutoff) == 39
    # Identity-checked finite PMDI cycle: independent sum of stored state G values.
    g = neighbor['by_cutoff']['50']['state_G_eV']
    ddg = math.fsum([g['pair'],g['empty'],-g['site2'],-g['site3']])*96.48533212331002
    assert math.isclose(ddg,neighbor['deltaDeltaG_kJ_mol'],abs_tol=1e-8)
    sources[str(Path(__file__).relative_to(ROOT))] = digest(__file__)
    snapshot = dict(status='discussion_evidence_arithmetic_and_identity_checked', created=time.time(),
        enumeration=dict(exact_designs=len(pool),recipes=recipes,families=dict(families),
                         evaluated_overlap=len(candidate_graphs & pool_graphs),
                         quantum_complete_overlap=sum(canonical(r['exact_grafted_smiles']) in pool_graphs
                             and isinstance(r.get('E1_vs_Fc_V'),(int,float))
                             and isinstance(r.get('E2_vs_Fc_V'),(int,float)) for r in candidates)),
        evaluation=dict(candidate_count=75,comparison_controls=len(feed['comparison_controls']),
            native_potential_candidates=sum(r['solution_potentials_status']=='verified native calculation' for r in candidates),
            two_potential_candidates=sum(isinstance(r.get('E1_vs_Fc_V'),(int,float))
                and isinstance(r.get('E2_vs_Fc_V'),(int,float)) for r in candidates),
            matched_comparisons=len(matched),systematic_matched_comparisons=len(systematic),
            as_of=feed['as_of_Chicago'],source_arithmetic_checks=feed['summary']['arithmetic_checks'],
            source_discrepancies=feed['summary']['arithmetic_discrepancies']),
        systematic_matches=systematic,matched_comparisons=matched,
        exact_examples=[by_id[k] for k in ['design_pmdi_ethyl','search_pmdi_e5caac5870e1','search_pmdi_1864c2c74267']],
        expansion_plan=plan,benchmark=dict(metrics=checked_metrics,points=endpoints,
            by_family=benchmark['by_family'],by_wave=benchmark['by_wave'],
            molecule_bootstrap=bias['uncertainty'][0],molecule_held_out=cv,
            held_out_external=False,exact_graft_accuracy_validated=False),
        last_week=dict(benchmark_endpoints=baseline['admitted_events'],
            spacing_selected_native_G=baseline_polymer['qualified_native_G'],
            original_oligomers=baseline_polymer['original_checked']),
        sensitivity=dict(basis_rows=basis['rows'],rotor_cutoff_rows=cutoff,
            maximum_benchmark_cutoff_span_mV=max(float(r['range_mV']) for r in cutoff)),
        base_assisted=base,PMDI_neighbor=neighbor,neighbor_archive=archive,
        selection=dict(summary=policy,descriptor_front=pareto,
            descriptor_candidates=pareto_feed['original39']+pareto_feed['additional_candidates']),
        UMA=dict(nitro=nitro,v100=v100), source_sha256=sources,
        limits=['Enumeration is not a completed quantum screen.',
            'Matched computed shifts are not experimentally validated improvements.',
            'Bootstrap and model settings spans are not candidate prediction intervals.',
            'Finite native cycles are not converged conformer ensembles.',
            'No candidate currently passes every required material-selection axis.'])
    for name,expected in sources.items():
        if digest(ROOT/name) != expected:
            raise RuntimeError('Source changed while preparing discussion: '+name)
    out.mkdir(parents=True)
    (out/'evidence.json').write_text(json.dumps(snapshot,indent=2)+'\n')
    print(json.dumps(dict(output=str(out), enumeration=snapshot['enumeration'],evaluation=snapshot['evaluation'],
                         benchmark=checked_metrics,source_files=len(sources)),indent=2))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    args=parser.parse_args()
    collect(args.out.resolve())
