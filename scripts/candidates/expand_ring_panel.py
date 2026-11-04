#!/usr/bin/env python
"""Enumerate a predeclared one-substitution panel retaining every site recipe."""
import argparse
import hashlib
import json
from pathlib import Path
from rdkit import Chem
from redox.build.contracts import validate_group
from redox.build.grafting import precursor
from redox.build.attachment_identity import exact_attachment
from redox.core.io import atomic_write_json
from redox.screening.family_search import ring_variants


def prepare(spec):
    baseline=spec['baseline']; records={}
    for label, site, fragment in ring_variants(baseline['frag'],spec['substituents']):
        row=dict(baseline,frag=fragment)
        molecule=validate_group('Cc1ccc(C[*:1])cc1',row)
        product=Chem.MolToSmiles(molecule,isomericSmiles=True)
        key=hashlib.sha256(product.encode()).hexdigest()
        row.update(id=baseline['id'] if label=='H' else 'proposal_ptz_'+key[:16],
                   name=baseline['name']+' / '+label, product_smiles=product,
                   precursor_smiles=Chem.MolToSmiles(precursor(fragment)[0]),
                   baseline_id=baseline['id'],hypothesis=spec['hypothesis'],
                   properties_computed=False,route_scope=spec['route_scope'],recipes=[])
        identity=dict(scaffold_smiles='Cc1ccc(C[*:1])cc1',attachment_fragment=fragment,
                      product_smiles=product,precursor_smiles=row['precursor_smiles'])
        payload=dict(attachment_fragment=fragment,product_smiles=product,
                     precursor_smiles=row['precursor_smiles'],family=row['family'],states=row['states'])
        attachment=exact_attachment(dict(identity=identity,payload=payload))
        row.update(route_binding=attachment['route_binding'],exact_design_key=attachment['exact_design_key'],
                   reference_reaction=attachment['balanced_reaction'])
        recipe=dict(axis='single_ring_substitution',substituent=label,template_atom=site,
                    baseline_id=baseline['id'])
        if key not in records: records[key]=row
        records[key]['recipes'].append(recipe)
    return list(records.values())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--spec',default='config/candidate_ring_expansion.json');p.add_argument('--out',required=True);a=p.parse_args()
    import os
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    if Path(a.out).exists(): raise ValueError('Preserve existing proposal panel')
    rows=prepare(json.loads(Path(a.spec).read_text()));atomic_write_json(a.out,rows)
    print(json.dumps(dict(exact_graphs=len(rows),recipes=sum(len(r['recipes']) for r in rows),quantum_jobs=0)))
