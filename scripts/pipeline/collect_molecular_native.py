#!/usr/bin/env python3
"""Join explicitly selected checked namespaces without relabeling older failures."""
import argparse
import json
from pathlib import Path

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.molecular_native import MolecularSources
from redox.workflows.molecular_native_namespace import verify


def collect(prepared,collections,output):
    if output.exists():raise FileExistsError('Preserve earlier molecular collections')
    wanted=json.loads(prepared.read_text())['namespaces']
    if len(wanted)<2 or len(wanted)!=len(set(wanted)):raise ValueError('Distinct explicitly selected final namespaces required')
    sources={str(prepared):sha256_file(prepared)};candidates={};history=[]
    for path in collections:
        sources[str(path)]=sha256_file(path)
        for row in json.loads(path.read_text())['entries']:
            history.append(dict(collection=str(path),namespace=row['namespace'],state=row['state'],status=row['status']))
            if row['namespace'] not in wanted or row['status']!='exact_molecular_native_G_independently_checked':continue
            if row['namespace'] in candidates:raise ValueError('Ambiguous duplicate checked namespace')
            candidates[row['namespace']]=row
    if set(candidates)!=set(wanted):raise ValueError('A selected native state is still missing or unqualified')
    rows=[];identities=set();charges=set()
    for namespace in wanted:
        row=candidates[namespace];path=Path(namespace);spec=verify(path)
        review_path=Path(row['review'])
        if sha256_file(review_path)!=row['review_sha256']:raise ValueError('Checked review changed')
        review=json.loads(review_path.read_text())
        if (review['status']!='exact_molecular_native_G_independently_checked'
                or review['namespace']!=namespace or review['namespace_sha256']!=sha256_file(path)
                or review['state']!=spec['state'] or review['G_by_cutoff_eV']!=row['G_by_cutoff_eV']):
            raise ValueError('State identity or review arithmetic differs')
        reader=MolecularSources(spec);reader.verify_all(review['source_sha256']);reader.recheck()
        identities.add((review['candidate_alias'],review['product_inchikey']))
        if review['charge'] in charges:raise ValueError('Declare one selected spin/conformer model per charge')
        charges.add(review['charge']);sources[str(review_path)]=sha256_file(review_path);rows.append(row)
    if len(identities)!=1:raise ValueError('The collection must refer to one exact candidate graph')
    if any(sha256_file(p)!=h for p,h in sources.items()):raise ValueError('Collection inputs changed during assembly')
    result=dict(status='explicit_recovered_native_state_collection_checked',entries=rows,
        expected=len(rows),checked=len(rows),terminal=len(rows),source_sha256=sources,
        original_outcome_history=history,older_failures_relabelled=False,material_admission=False)
    atomic_write_json(output,result);return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--prepared',type=Path,required=True)
    p.add_argument('--collections',type=Path,nargs='+',required=True);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    import os
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    result=collect(a.prepared.resolve(),[p.resolve() for p in a.collections],a.out.resolve())
    print(json.dumps(dict(status=result['status'],checked=result['checked'])))
