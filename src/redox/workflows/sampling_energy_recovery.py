"""Replace one rejected electronic slot with an exact-state fresh recovery.

The original namespace and its interrupted/exclusivity evidence remain intact.
Only its checked same-frame density initializes a new solve; every normal
convergence, residency, resource and independent review gate remains required.
"""
import argparse
import copy
import json
import os
from pathlib import Path
import shutil

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.molecular_native_namespace import verify
from redox.workflows.sampling_energy import density_restart_source


def prepare(manifest_path,key,output,successor):
    manifest=json.loads(manifest_path.read_text())
    if manifest['target_states']!=240 or len(manifest['states'])!=240 or successor.exists():
        raise ValueError('Complete cohort and fresh immutable successor manifest required')
    row=next(r for r in manifest['states'] if r['key']==key);namespace=Path(row['namespace']);spec=verify(namespace)
    old_case=Path(spec['case']);old_sources=Path(spec['sources']);old=json.loads((old_case/'plan.json').read_text())
    proof=density_restart_source(namespace,old['tasks'][0],old['experimental_method_id'],old['scored_geometry_sha256'],reason='foreign_residency')
    if output.exists() or any(output.is_relative_to(Path(spec[k])) for k in ('prior_root','archive','case','sources')):
        raise FileExistsError('Fresh recovery namespace outside every scientific source required')
    output.mkdir(parents=True);case=output/'case';sources=output/'sources';case.mkdir()
    shutil.copytree(old_sources,sources,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    shutil.copyfile(old_case/'input.xyz',case/'input.xyz');shutil.copyfile(proof['checkpoint'],sources/'restart_occupied.h5')
    atomic_write_json(sources/'rejected_source.json',proof);shutil.copyfile(__file__,sources/'recovery_preparer.py')
    plan=copy.deepcopy(old)
    plan.update(raw_root=str(case/'raw'),density_restart=dict(solvent='acetonitrile',
        checkpoint=str(sources/'restart_occupied.h5'),sha256=sha256_file(sources/'restart_occupied.h5'),
        source_receipt=str(sources/'rejected_source.json'),source_receipt_sha256=sha256_file(sources/'rejected_source.json')),
        source_rejected_namespace=str(namespace),source_qualification_transferred=False)
    hashes={**old['source_sha256'],**{str(p):sha256_file(p) for p in [*sources.glob('*'),case/'input.xyz',manifest_path] if p.is_file()}}
    plan.update(source_sha256=hashes,producer_sha256=hashes,input_sha256={});atomic_write_json(case/'plan.json',plan)
    new=dict(spec,case=str(case),sources=str(sources),plan_sha256=sha256_file(case/'plan.json'))
    atomic_write_json(output/'namespace.json',new);verify(output/'namespace.json')
    result=copy.deepcopy(manifest)
    replacement=next(r for r in result['states'] if r['key']==key)
    replacement.update(namespace=str(output/'namespace.json'),namespace_sha256=sha256_file(output/'namespace.json'))
    result.update(status='complete_corrected_cohort_with_source_bound_rejected_slot_recovery',
        predecessor_manifest=str(manifest_path),predecessor_manifest_sha256=sha256_file(manifest_path),
        recovered_slot=dict(key=key,rejected_namespace=str(namespace),replacement_namespace=str(output/'namespace.json'),
                            source_qualification_transferred=False))
    atomic_write_json(successor,result);return result['recovered_slot']


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('manifest','out','successor'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--key',required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    print(json.dumps(prepare(a.manifest.resolve(),a.key,a.out.resolve(),a.successor.resolve())))
