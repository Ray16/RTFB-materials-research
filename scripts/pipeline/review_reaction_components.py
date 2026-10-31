#!/usr/bin/env python3
"""Review finished reagent conformers without accepting process exit as QC."""
import argparse
import json
import os
from pathlib import Path

from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.molecular_native import review
from redox.workflows.molecular_native_namespace import verify


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--prepared',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    a.out.mkdir(parents=True,exist_ok=False);rows=[]
    paths=[Path(p) for p in json.loads(a.prepared.read_text())['namespaces']]
    for path in paths:
        spec=verify(path);case=Path(spec['case']);state=case/'status.json'
        row=dict(namespace=str(path),namespace_sha256=sha256_file(path),status='not_started')
        if state.exists():
            data=json.loads(state.read_text());row['status']=data['status']
            # The closing guard files are mandatory in review(); partial work
            # is never used in a cycle or promoted to a candidate quantity.
            if data['status']=='reaction_component_native_finished_independent_review_pending':
                output=a.out/(path.parent.name+'_independent_review.json')
                result=review(path,output,reaction_component=True)
                row.update(status=result['status'],review=str(output),review_sha256=sha256_file(output),G_by_cutoff_eV=result['G_by_cutoff_eV'])
        rows.append(row)
    result=dict(status='reaction_component_review_inventory',entries=rows,expected=len(paths),
        checked=sum(r['status']=='reaction_component_native_G_independently_checked' for r in rows),
        material_admission=False,quantum_calculations_run=False)
    atomic_write_json(a.out/'summary.json',result);print(json.dumps(result))


if __name__=='__main__':main()
