#!/usr/bin/env python
"""Join exact enumerated PMDI/NDI evidence, retaining incomplete dimensions."""
import argparse,hashlib,importlib.util,json,os
from pathlib import Path
from redox.core.common import ROOT,read_manifest
from redox.core.io import atomic_write_json
from redox.screening.route_stage import screen
from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,required=True);args=parser.parse_args()
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    args.out.mkdir(parents=True,exist_ok=False)
    identity=ROOT/'results/reviews/2026-10-22/end_to_end_selection/identity.json'
    paths=[identity,ROOT/'results/grafting/route_evidence.json',ROOT/'results/grafting/base_assisted/current.json',ROOT/'results/candidates/current.json',Path(__file__)]
    sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
    sources={str(p):sha(p) for p in paths};cases=json.loads(identity.read_text())['cases']
    evidence=json.loads(paths[1].read_text())['candidates'];grafting=json.loads(paths[2].read_text())['candidates'];feed=json.loads(paths[3].read_text())
    groups={r['id']:r for r in feed['original39']+feed['additional_candidates']}
    manifest={(r['id'],r['state']):r for r in read_manifest()}
    spec=importlib.util.spec_from_file_location('exact_case_collector',ROOT/'scripts/pipeline/collect_verified_candidates.py')
    collector=importlib.util.module_from_spec(spec);spec.loader.exec_module(collector);audit=collector.audit_module()
    def state(key,label,reference=False):
        row=manifest[key,label]
        value,qc,backing=collector.independent_state(key,label,int(row['charge']),collector.coverage.full_key(row['smiles']),audit)
        if value is None:raise ValueError('Exact state does not independently reconstruct: '+key+':'+label+' '+qc)
        sources.update(backing)
        if not reference and not collector.state_spin_evidence(key,label)['spin_primary_admissible']:
            raise ValueError('Exact candidate spin selection remains unqualified')
        return value
    fc=state('ferrocene','ox',True)-state('ferrocene','neu',True);rows=[]
    for case in cases:
        key=case['candidate_id'];group=groups[key];binding=case['exact_attachment']
        route=screen(dict(payload=dict(id=case['enumerated_id'],fragment=group['frag'],reviewed=evidence[key],assessment={},attachment=binding)))
        if not route['exact_route_binding_verified'] or route['proceed_to_full_redox']:raise ValueError('Incomplete exact route unexpectedly proceeds')
        free={label:state(key,label) for label in ('neu','red1','red2')}
        reaction=next(r for r in grafting if r['id']==key)
        if reaction['reaction']['product_smiles']!=case['product_smiles'] or reaction['reaction']['precursor_smiles']!=case['precursor_smiles']:
            raise ValueError('Grafting product/precursor differs from exact selected case')
        for path,h in reaction['source_sha256'].items():
            if sha(path)!=h:raise ValueError('Grafting backing evidence changed')
            sources[path]=h
        rows.append(dict(candidate_id=key,enumerated_id=case['enumerated_id'],exact_design_key=binding['exact_design_key'],
            state_G_eV=free,E1_V_vs_Fc=free['neu']-free['red1']-fc,E2_V_vs_Fc=free['red1']-free['red2']-fc,
            redox_QC='exact native states and spin selection independently checked; external graft accuracy unvalidated',
            reference_grafting_kJ_mol=reaction['reference_deltaG_kJ_mol'],base_models=reaction['base_models'],
            route=route,conditional_occupancy_deltaDeltaG_kJ_mol=None,occupancy_QC='matched exact-motif occupancy cycle and conformer sampling pending',
            synthesis_decision='hold',reason='Favorable modeled grafting does not override missing exact precursor/selectivity/redox-retention evidence; conditional polymer assessment pending',
            end_to_end_complete=False,quantum_jobs_submitted=0))
    if any(sha(p)!=h for p,h in sources.items()):raise ValueError('Evidence changed while joining cases')
    atomic_write_json(args.out/'summary.json',dict(status='two_exact_cases_joined_with_explicit_route_and_occupancy_holds',cases=rows,
        Fc_absolute_V=fc,source_sha256=sources,end_to_end_complete=False))
    print(json.dumps(dict(cases=len(rows),decisions=[r['synthesis_decision'] for r in rows],end_to_end_complete=False)))


if __name__=='__main__':main()
