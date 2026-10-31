"""Exact-identity warm starts for corrected empty-backbone reference calculations."""
import argparse
import json
import os
from pathlib import Path
import shutil

from rdkit import Chem

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.charged_polymer_pilot import prepare_optimization
from redox.workflows.molecular_native_namespace import load


def select_seed(analysis,case,solvent_model):
    data=json.loads(analysis.read_text());choices=[]
    prefix=case['configuration']+'__occ_none__'
    canonical=lambda value:Chem.MolToSmiles(Chem.MolFromSmiles(value),isomericSmiles=True)
    if case['occupied_sites'] or case['charge']!=0:raise ValueError('A matched unoccupied graph is required')
    for key,record in data['inputs'].items():
        if not key.startswith(prefix):continue
        review_path=Path(record['review']);geometry=Path(record['geometry'])
        if sha256_file(review_path)!=record['review_sha256'] or sha256_file(geometry)!=record['geometry_sha256']:
            raise ValueError('Changed saved empty-reference evidence')
        review=json.loads(review_path.read_text());task=review['task']
        if (review['status']!='alternate_native_state_independently_reconstructed'
                or review['experimental_model']!=solvent_model or task['charge']!=0 or task['mult']!=1
                or task['occupied_sites'] or task['graph_sha256']!=case['graph_sha256']
                or canonical(task['smiles'])!=canonical(case['smiles'])
                or task['repeat_units']!=case['repeat_units'] or task['reactive_sites']!=case['reactive_sites']):
            raise ValueError('Saved reference differs in graph/state/model identity')
        value=review['G_by_cutoff_eV']['50']
        import math
        if not math.isfinite(value):raise ValueError('Finite saved seed-selection value required')
        choices.append((value,key,record))
    if not choices:raise ValueError('No exact matched empty-reference warm start')
    value,key,record=min(choices,key=lambda row:(row[0],row[1]))
    return dict(key=key,geometry=record['geometry'],geometry_sha256=record['geometry_sha256'],
        review=record['review'],review_sha256=record['review_sha256'],selection='lowest stored G at 50 cm-1 within the declared old finite batch',
        old_G_used_for_seed_selection_only=True,old_E_or_G_admitted_to_new_cycle=False,candidates_checked=len(choices),
        analysis=str(analysis),analysis_sha256=sha256_file(analysis))


def prepare(archive,preparation,agreement,analysis,output,key):
    result=prepare_optimization(archive,preparation,agreement,output,key,allow_empty_reference=True)
    case=output/'case';sources=output/'sources';plan=json.loads((case/'plan.json').read_text())
    if not plan['matched_empty_polymer_reference']:raise ValueError('The selected model is not an explicit empty reference')
    seed=select_seed(analysis,plan['prepared_case'],plan['solvent_model']);shutil.copyfile(seed['geometry'],sources/'input.xyz')
    shutil.copyfile(analysis,sources/'warm_start_analysis.json')
    shutil.copyfile(seed['review'],sources/'warm_start_review.json')
    seed['snapshot_files']=dict(analysis='warm_start_analysis.json',review='warm_start_review.json',geometry='input.xyz')
    checker=load(ROOT/'scripts/validation/review_charged_polymer_optimization.py','matched_empty_geometry_checks')
    atomic_write_json(sources/'reference_geometry_review.json',checker.graph_checks(plan['prepared_case']['smiles'],sources/'input.xyz',sources/'input.xyz'))
    atomic_write_json(sources/'warm_start_lineage.json',seed)
    plan['source_sha256'].update({str(p):sha256_file(p) for p in sources.glob('*')})
    plan.update(selected_initial_geometry_sha256=sha256_file(sources/'input.xyz'),warm_start=seed,
        purpose='Fresh corrected matched empty reference; old thermochemical values select geometry only')
    atomic_write_json(case/'plan.json',plan);spec=json.loads((output/'namespace.json').read_text());spec['plan_sha256']=sha256_file(case/'plan.json');atomic_write_json(output/'namespace.json',spec)
    result['status']='matched_empty_reference_fresh_optimization_prepared';return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('archive','preparation','agreement','analysis','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--key',required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    print(json.dumps(prepare(a.archive.resolve(),a.preparation.resolve(),a.agreement.resolve(),a.analysis.resolve(),a.out.resolve(),a.key)))
