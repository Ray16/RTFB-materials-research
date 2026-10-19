#!/usr/bin/env python
"""Attach current minimum/conformer evidence to immutable solution diagnostics.

No raw calculation record, production ranking or geometry is replaced. A rigid
heavy graph permits a same-basin classification only; minimum qualification
still requires addressed, derivative-checked native solution Hessian evidence.
"""
import argparse
import hashlib
import json
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Lipinski

from redox.core.common import ROOT,DFT,read_manifest
from redox.core.protocol import geom_sha1
from redox.properties.solution_reorg import summarize_relaxations
from redox.qm.native_thermal import load_native_thermal,cache_path
from redox.validation.solution_relaxation import verified_implementation_controls
from redox.core.io import atomic_write_json

OUT=ROOT/'results/validation/physical_reorg_2026-10-11'


def main():
    p=argparse.ArgumentParser();p.add_argument('--limit',type=int,default=4)
    p.add_argument('--out',type=Path,help='new derived review path; existing outputs are preserved')
    a=p.parse_args()
    target=a.out or OUT/'native_smd_qualification_review.json'
    if target.exists():p.error('choose a new --out; existing qualification reviews are preserved')
    controls=[OUT/'native_smd_independent_reference.json',OUT/'native_smd_gpu_reference.json']
    implementation=verified_implementation_controls({str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in controls},
        ROOT/'src/redox/qm/nonequilibrium_smd.py',backend='gpu')
    queue=json.loads((OUT/'native_smd_queue.json').read_text())
    manifest={r['id']:r for r in read_manifest()}
    cases=[]
    for c in queue['cases'][:a.limit]:
        source=OUT/f'native_smd_{c["id"]}_{c["oxidized"]}_{c["reduced"]}.json'
        if not source.exists():
            cases.append(dict(id=c['id'],status='pending'));continue
        record=json.loads(source.read_text());stale=[];evidence=[]
        if record.get('code_sha256',{}).get('src/redox/qm/nonequilibrium_smd.py')!=implementation['module_sha256']:
            cases.append(dict(id=c['id'],status='implementation_source_stale'));continue
        mol=Chem.MolFromSmiles(manifest[c['id']]['smiles'])
        rotors=mol.GetSubstructMatches(Lipinski.RotatableBondSmarts)
        for i in c['states']:
            sd=DFT/i['id']/i['state']
            state=dict(state=i['state'],conformer_policy='same_basin' if not rotors else 'different_minima',
                conformer_basis='rigid heavy-atom graph' if not rotors else 'flexible heavy torsions; paired/ensemble evidence pending',
                heavy_rotatable_bonds=len(rotors),minimum_verified=False)
            if geom_sha1(sd/'opt.xyz')!=i['geom_sha1']:
                stale.append(i['state']);state['source_geometry_stale']=True
            else:
                thermal=load_native_thermal(sd,i['charge'],i['mult'])
                if thermal:
                    path,_=cache_path(sd,i['charge'],i['mult'])
                    state.update(minimum_verified=True,minimum_source=str(path),
                        minimum_source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                        minimum_energy_surface=thermal['identity']['protocol'],
                        active_SP_stationarity_verified=False,
                        minimum_frequency_cm=min(thermal['frequencies_cm']),
                        gradient_max_Ha_Bohr=thermal['gradient_max'])
                else:state['minimum_evidence_status']='pending_derivative_checked_solution_Hessian'
            evidence.append(state)
        summary=summarize_relaxations(*record['directional'],endpoint_evidence=evidence,
                                     method_identity=record['method'],implementation_validation=True)
        flags={f for e in record['electronic_state_evidence'] for f in e['flags']}
        if flags:summary['physical_flags']+=sorted(flags);summary['endpoint_model_verified']=False
        if stale:summary['physical_flags'].append('source_geometry_stale');summary['endpoint_model_verified']=False
        # Unequal directional curvatures do not invalidate their individually
        # calculated relaxation energies; their average is a separate assumption.
        summary['directional_descriptor_model_qualified']=summary['endpoint_model_verified']
        summary['qualification_convention']='composite active-SP energies at geometry-level solution minima'
        summary['active_SP_stationary_endpoints_proven']=False
        summary['symmetric_Marcus_approximation_qualified']=summary['endpoint_model_verified'] and (
            'equal_curvature_approximation_unverified' not in summary['physical_flags'])
        item=dict(id=c['id'],couple=f'{c["oxidized"]}->{c["reduced"]}',source=str(source),
            source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
            endpoint_evidence=evidence,summary=summary,
            raw_record_modified=False,production_promotion=False)
        if c['kind']=='matched_experimental_benchmark':
            mean=summary['symmetric_single_body_approximation_eV']
            item.update(experimental_interval_eV=[1.26,1.28],benchmark_difference_to_midpoint_eV=mean-1.27,
                benchmark_comparison_status='reported_with_CV_fit_and_double_layer_limitations')
        cases.append(item)
    result=dict(cases=cases,expected_cases=min(a.limit,len(queue['cases'])),
        endpoint_models_qualified=sum(c.get('summary',{}).get('endpoint_model_verified',False) for c in cases),
        symmetric_models_qualified=sum(c.get('summary',{}).get('symmetric_Marcus_approximation_qualified',False) for c in cases),
        production_promotion=False,implementation_evidence=implementation,
        per_case_independent_replay_not_inferred_from_controls=True)
    atomic_write_json(target,result)
    print(json.dumps(result))


if __name__=='__main__':main()
