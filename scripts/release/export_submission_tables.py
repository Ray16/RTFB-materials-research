#!/usr/bin/env python3
"""Export exact design identities and the frozen retrospective benchmark as CSV."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path

from rdkit import Chem

from redox.core.hashing import sha256_file


def export(registry_path, benchmark_path, state_path, protocol_path, output):
    if output.exists():raise FileExistsError('Preserve earlier submission tables')
    registry=json.loads(registry_path.read_text());benchmark=json.loads(benchmark_path.read_text())
    states=json.loads(state_path.read_text());protocol=json.loads(protocol_path.read_text())['package']
    aliases={}
    for row in states['states']:
        for alias in row['aliases']:
            value=dict(charge=row['charge'],multiplicity=row['mult'])
            if alias in aliases and aliases[alias]!=value:raise ValueError('Ambiguous benchmark state')
            aliases[alias]=value
    context=dict(id='frozen_parent_benchmark_primary_model',quantity='one_electron_couple_potential',unit='V',
        solvent='acetonitrile',temperature_K=protocol['native_protocol']['temperature_K'],
        electronic_protocol=protocol['active_protocol'],native_protocol=protocol['native_protocol'],
        thermal_definition='Projected unscaled native solution-DFT modes, xTB qRRHO at 50 cm-1; ground spin multiplicity once',
        reference=dict(electrode='Fc/Fc+',neutral=aliases['benchmark:ferrocene:neu'],oxidized=aliases['benchmark:ferrocene:ox']),
        standard_state='Same source ideal-translation convention; balanced one-electron Fc cycle cancels common 1-M conversion',
        numerical_scope='Frozen reported benchmark values; later electronic implementation repeats do not replace nuclear terms',
        QC_scope='Export rechecks identity and residual arithmetic only; source validation is separate',
        graph_scope='Parent SMILES is a connectivity representation; model charge/spin are explicit columns, not localization assignments',
        experimental_exact_graft_validation=False)
    designs=[]
    for row in registry['candidates']:
        mol=Chem.MolFromSmiles(row['product_smiles'])
        if mol is None or any(a.GetAtomicNum()==0 for a in mol.GetAtoms()):
            raise ValueError('A complete exact product graph is required')
        smiles=Chem.MolToSmiles(mol,isomericSmiles=True)
        if (hashlib.sha256(smiles.encode()).hexdigest()!=row['graph_key']
                or Chem.MolToInchiKey(mol)!=row['inchikey']):
            raise ValueError('Registry identity differs from its SMILES')
        designs.append(dict(graph_key=row['graph_key'],product_SMILES=smiles,
            product_InChIKey=row['inchikey'],aliases_json=json.dumps(row['aliases']),
            source_roles_json=json.dumps(row['roles']),route_status=row['route_status'],
            reported_evaluation_count=len(row['evaluations']),material_decision=row['material_decision'],
            scope='capped local graft design; enumeration is not experimental or material validation'))
    endpoints=[]
    for row in benchmark['outliers']:
        if Chem.MolFromSmiles(row['parent_smiles']) is None:raise ValueError('Benchmark parent SMILES invalid')
        if not all(math.isfinite(row[k]) for k in ('exp_V','calc_V','error_V')):
            raise ValueError('Finite benchmark values required')
        if not math.isclose(row['calc_V']-row['exp_V'],row['error_V'],abs_tol=1e-12):
            raise ValueError('Benchmark residual does not close')
        start,end=row['event'].split('->')
        first=aliases['benchmark:'+row['id']+':'+start];second=aliases['benchmark:'+row['id']+':'+end]
        endpoints.append(dict(id=row['id'],name=row['name'],parent_SMILES=row['parent_smiles'],
            family=row['family'],event=row['event'],wave=row['wave'],
            model_initial_charge=first['charge'],model_initial_multiplicity=first['multiplicity'],
            model_final_charge=second['charge'],model_final_multiplicity=second['multiplicity'],
            model_solvent=context['solvent'],model_temperature_K=context['temperature_K'],
            context_file='benchmark_context.json',context_id=context['id'],
            experimental_V_vs_Fc=row['exp_V'],calculated_V_vs_Fc=row['calc_V'],
            signed_error_V=row['error_V'],source_DOI=row['doi'],
            scope='retrospective parent benchmark; not exact-graft experimental validation'))
    if len({r['graph_key'] for r in designs})!=len(designs):raise ValueError('Duplicate design graph')
    if len({(r['id'],r['event'],r['source_DOI']) for r in endpoints})!=len(endpoints):
        raise ValueError('Duplicate benchmark endpoint')
    output.mkdir(parents=True)
    (output/'benchmark_context.json').write_text(json.dumps(context,indent=2)+'\n')
    for name,rows in [('design_SMILES.csv',designs),('benchmark_endpoints.csv',endpoints)]:
        with (output/name).open('w',newline='') as stream:
            writer=csv.DictWriter(stream,fieldnames=list(rows[0]),lineterminator='\n');writer.writeheader();writer.writerows(rows)
    result=dict(status='submission_tables_identity_and_arithmetic_checked',designs=len(designs),
        benchmark_endpoints=len(endpoints),benchmark_molecules=len({r['id'] for r in endpoints}),
        source_sha256={str(p):sha256_file(p) for p in (registry_path,benchmark_path,state_path,protocol_path,Path(__file__))},
        files_sha256={p.name:sha256_file(p) for p in output.iterdir()},
        original_quantum_QC_repeated=False,material_admission_inferred=False)
    (output/'manifest.json').write_text(json.dumps(result,indent=2)+'\n');return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('registry','benchmark','states','protocol','out'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    r=export(a.registry.resolve(),a.benchmark.resolve(),a.states.resolve(),a.protocol.resolve(),a.out.resolve())
    print(json.dumps({k:r[k] for k in ('status','designs','benchmark_endpoints','benchmark_molecules')}))
