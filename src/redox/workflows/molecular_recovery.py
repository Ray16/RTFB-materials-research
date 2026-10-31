"""Reuse verified minimum-repair geometries after a post-QM namespace failure.

Old attempts remain unchanged. Fresh native checks run in a new guarded case;
successful optimization and active-SP work is retained instead of recomputed.
"""
import argparse
import json
from pathlib import Path
import shutil

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.core.protocol import load_record,record_path
from redox.validation.protonation import hydrogen_topology_status
from redox.workflows.molecular_native_namespace import prepare,verify


def recovered_metadata(raw,original_state,source_record,target):
    """Relocate provenance only; do not claim the old optimizer produced a new SP."""
    result=dict(raw)
    result.pop('energy_protocol',None);result.pop('energy_geom_sha1',None)
    result.update(opt_out=str(target/'opt.xyz'),minimum_repair_source=str(original_state/'opt.xyz'),
        legacy_thermal_archive=str(original_state/'result.json'),
        recovery_input_record=str(source_record),recovery_input_record_sha256=sha256_file(source_record),
        recovery_metadata_scope='provenance-path relocation; chemistry, geometry and addressed energy unchanged')
    return result


def select_verified_repair(spec):
    case=Path(spec['case']);plan=json.loads((case/'plan.json').read_text())
    state=json.loads((case/'status.json').read_text())
    if (state['status']!='molecular_native_failed_preserved'
            or state.get('error')!="OSError(18, 'Invalid cross-device link')" or 'finished' not in state):
        raise ValueError('Only the diagnosed closed post-QM cross-mount promotion failure is supported')
    guard=json.loads((case/'residency_guard/owner.json').read_text())
    if not guard.get('guard_stopped') or (case/'residency_guard/interruption.json').exists():
        raise ValueError('Old scientific ownership must be closed without interruption')
    statuses=list((case/'repairs'/spec['source_id']/spec['state']).glob('*/status.json'))
    if len(statuses)!=1:raise ValueError('One explicit completed two-branch repair required')
    path=statuses[0];repair=json.loads(path.read_text())
    if (len(repair['branches'])!=2 or {r['sign'] for r in repair['branches']}!={-1,1}
            or any(r['status'] not in ('verified','failed') for r in repair['branches'])
            or not any(r['status']=='verified' for r in repair['branches'])):
        raise ValueError('Both original signed-mode branches must finish before selection')
    candidates=[]
    prefix=Path(spec['prior_root'])/'results/validation/physical_thermal_spin_2026-10-10/native_minimum_repairs'
    for row in repair['branches']:
        if row['status']!='verified':continue
        geometry=Path(row['geometry'])
        if not geometry.is_relative_to(prefix):raise ValueError('Unexpected repaired geometry namespace')
        directory=case/'repairs'/geometry.parent.relative_to(prefix)
        record=load_record(directory,plan['charge'],plan['mult'])
        if (record is None or abs(record['e_smd_eV']-row['e_smd_eV'])>1e-8
                or hydrogen_topology_status(plan['product_smiles'],directory/'opt.xyz')!='ok'
                or not row['geometry_optimization_converged'] or not row['optimization_SCF_converged']):
            raise ValueError('Saved repair energy/graph/optimization does not reconstruct')
        candidates.append((record['e_smd_eV'],directory))
    selected=min(candidates,key=lambda r:r[0])[1]
    return plan,state,path,selected


def prepare_recovery(old_namespace,registry,routes,assessment,output,alive):
    old_namespace=Path(old_namespace).resolve();spec=verify(old_namespace)
    plan,state,repair_status,selected=select_verified_repair(spec)
    if alive(state['owner']) is not False:raise ValueError('Old scientific owner is live or unverified')
    source=[old_namespace,Path(spec['case'])/'status.json',repair_status,
        selected/'opt.xyz',selected/'result.json',record_path(selected,plan['charge'],plan['mult'])]
    before={str(p):sha256_file(p) for p in source}
    prepared=prepare(Path(spec['archive']),registry,routes,assessment,plan['candidate_alias'],plan['source_id'],[plan['state']],output)
    namespace=Path(prepared['namespaces'][0]);new=json.loads(namespace.read_text());case=Path(new['case'])
    target=case/'dft'/plan['source_id']/plan['state'];backup=case/'initial_source_before_recovery';backup.mkdir()
    for name in ('opt.xyz','result.json'):
        shutil.copyfile(target/name,backup/name);shutil.copyfile(selected/name,target/name)
    shutil.copyfile(selected/'result.json',Path(new['sources'])/'recovered_optimizer_record.json')
    preparer=Path(new['sources'])/'recovery_preparer.py';shutil.copyfile(Path(__file__).resolve(),preparer)
    metadata=recovered_metadata(json.loads((selected/'result.json').read_text()),
        Path(spec['case'])/'dft'/plan['source_id']/plan['state'],selected/'result.json',target)
    atomic_write_json(target/'result.json',metadata)
    energy=record_path(selected,plan['charge'],plan['mult']);shutil.copyfile(energy,target/'sp'/energy.name)
    updated=json.loads((case/'plan.json').read_text())
    updated.update(source_geometry_sha256=sha256_file(selected/'opt.xyz'),
        original_archive_geometry_sha256=plan['source_geometry_sha256'],
        recovered_from=dict(namespace=str(old_namespace),namespace_sha256=before[str(old_namespace)],
            selected_branch=str(selected),selection='lowest active SP energy among verified branches after both signed-mode attempts terminate',
            failure='post-QM cross-mount canonical promotion',native_cache_reused=False,
            optimization_and_active_SP_reused=True,fresh_native_derivatives_required=True),
        source_sha256={**updated['source_sha256'],**before,str(preparer):sha256_file(preparer),
            str(Path(new['sources'])/'recovered_optimizer_record.json'):sha256_file(Path(new['sources'])/'recovered_optimizer_record.json')})
    if any(sha256_file(p)!=h for p,h in before.items()):raise ValueError('Closed repair inputs changed during recovery preparation')
    atomic_write_json(case/'plan.json',updated);new['plan_sha256']=sha256_file(case/'plan.json');atomic_write_json(namespace,new)
    prepared.update(status='verified_repair_geometry_prepared_for_fresh_native_checks',source_namespace=str(old_namespace),
        source_namespace_sha256=before[str(old_namespace)],optimization_reused=True,new_quantum_jobs_submitted=0)
    atomic_write_json(Path(output)/'prepared.json',prepared);return prepared


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('source','registry','routes','assessment','out'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    import importlib.util,os
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    definition=importlib.util.spec_from_file_location('repaired_molecular_owner',ROOT/'scripts/pipeline/observe_archived_native_namespace.py')
    owners=importlib.util.module_from_spec(definition);definition.loader.exec_module(owners)
    print(json.dumps(prepare_recovery(a.source,a.registry,a.routes,a.assessment,a.out,owners.alive)))
