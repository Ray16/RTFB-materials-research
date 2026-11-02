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
from redox.workflows.molecular_native_namespace import verify,load
from redox.workflows.sampling_energy import density_restart_source


def guard_preflight_source(namespace):
    """Prove an exited initial guard never entered its scientific program."""
    spec=verify(namespace);case=Path(spec['case']);plan=json.loads((case/'plan.json').read_text())
    placements=list(case.parent.glob('placement*/completion.json'))
    if len(placements)!=1:raise ValueError('One exact closed initial-guard placement required')
    completion_path=placements[0];launch_path=completion_path.with_name('launch.json');owner_path=case/'residency_guard/owner.json'
    completion=json.loads(completion_path.read_text());launch=json.loads(launch_path.read_text());owner=json.loads(owner_path.read_text())
    alive=load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py','guard_recovery_exact_owners').alive
    identities=[owner,dict(host=launch['host'],pid=launch['pid'],process_start=launch['process_start']),launch['budget_observer_owner']]
    if (completion['returncode'] in (None,0) or completion['budget_stop'] is not None
            or any(alive(identity) is not False for identity in identities)
            or any((case/name).exists() for name in ('status.json','host_guard','tensor_runtime'))
            or Path(plan['raw_root']).exists() or (case/'residency_guard/interruption.json').exists()):
        raise ValueError('Closed guard failure before any scientific worker/SCF is required')
    log=completion_path.with_name('command.log');text=log.read_text()
    if 'TimeoutExpired' not in text or 'nvidia-smi' not in text:
        raise ValueError('The preserved initial GPU-query failure must be explicit')
    paths=[namespace,case/'plan.json',completion_path,launch_path,owner_path,log]
    return dict(status='closed_initial_guard_query_failure_science_never_started',checkpoint=None,
        namespace=str(namespace),source_qualification_transferred=False,
        source_sha256={str(p):sha256_file(p) for p in paths})


def prepare(manifest_path,key,output,successor,reason='foreign_residency'):
    manifest=json.loads(manifest_path.read_text())
    if (manifest['target_states']!=240 or len(manifest['states'])!=240
            or len({r['key'] for r in manifest['states']})!=240 or successor.exists()):
        raise ValueError('Complete cohort and fresh immutable successor manifest required')
    row=next(r for r in manifest['states'] if r['key']==key);namespace=Path(row['namespace']);spec=verify(namespace)
    old_case=Path(spec['case']);old_sources=Path(spec['sources']);old=json.loads((old_case/'plan.json').read_text())
    if reason=='guard_preflight':proof=guard_preflight_source(namespace)
    elif reason=='foreign_residency':proof=density_restart_source(namespace,old['tasks'][0],old['experimental_method_id'],old['scored_geometry_sha256'],reason=reason)
    else:raise ValueError('An explicit supported closed-failure recovery reason is required')
    if output.exists() or any(output.is_relative_to(Path(spec[k])) for k in ('prior_root','archive','case','sources')):
        raise FileExistsError('Fresh recovery namespace outside every scientific source required')
    output.mkdir(parents=True);case=output/'case';sources=output/'sources';case.mkdir()
    shutil.copytree(old_sources,sources,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    shutil.copyfile(old_case/'input.xyz',case/'input.xyz')
    if proof['checkpoint'] is not None:shutil.copyfile(proof['checkpoint'],sources/'restart_occupied.h5')
    atomic_write_json(sources/'rejected_source.json',proof);shutil.copyfile(__file__,sources/'recovery_preparer.py')
    plan=copy.deepcopy(old)
    restart=(dict(solvent='acetonitrile',
        checkpoint=str(sources/'restart_occupied.h5'),sha256=sha256_file(sources/'restart_occupied.h5'),
        source_receipt=str(sources/'rejected_source.json'),source_receipt_sha256=sha256_file(sources/'rejected_source.json')) if proof['checkpoint'] is not None else None)
    plan.update(raw_root=str(case/'raw'),density_restart=restart,
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
    p.add_argument('--key',required=True);p.add_argument('--reason',choices=('foreign_residency','guard_preflight'),default='foreign_residency');a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    print(json.dumps(prepare(a.manifest.resolve(),a.key,a.out.resolve(),a.successor.resolve(),a.reason)))
