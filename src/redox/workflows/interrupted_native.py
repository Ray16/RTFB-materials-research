"""Continue an exact native state interrupted by a proven foreign GPU resident."""
import argparse
import json
import os
from pathlib import Path
import re
import shutil
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.namespace_sources import NamespaceSources
from redox.workflows.archived_native_namespace import verify_spec
from redox.workflows.submission_refusal import producer_inventory


def prepare(source_namespace,destination,alive,hosts):
    source=verify_spec(source_namespace);reader=NamespaceSources(source);case=Path(source['attempt'])
    if destination.exists() or any(destination.resolve().is_relative_to(Path(source[k]).resolve()) for k in ('prior_root','archive','attempt','raw_parent')):
        raise ValueError('Use a fresh destination outside every source tree')
    stage_path=case/'experimental_native/status.json';guard_path=case/'residency_guard/owner.json';interruption_path=case/'residency_guard/interruption.json'
    stage=json.loads(stage_path.read_text());guard=json.loads(guard_path.read_text());interruption=json.loads(interruption_path.read_text())
    if (stage['status']!='experimental_native_failed_preserved' or 'finished' not in stage
            or not guard.get('guard_stopped') or not interruption.get('foreign_resident_PIDs')
            or interruption.get('foreign_processes_signaled') is not False
            or interruption.get('signal_target')!='own reserved process only'
            or guard['pid']!=stage['owner']['pid'] or guard['process_start']!=stage['owner']['process_start']
            or interruption['pid']!=guard['pid'] or 'KeyboardInterrupt' not in stage.get('error','')):
        raise ValueError('A closed owned interruption from a foreign GPU resident is required')
    if alive(stage['owner']) is not False:raise ValueError('Original scientific ownership is live or unverified')
    inventory=producer_inventory(source,source_namespace,hosts)
    if set(inventory)!=set(hosts) or any(inventory.values()):raise ValueError('A matching producer is live or fleet ownership is incomplete')
    plan=json.loads((case/'plan.json').read_text());reader.verify_all(plan['producer_sha256']);reader.verify_all(plan['input_sha256'])
    if (len(plan['tasks'])!=1 or (plan['tasks'][0]['charge'],plan['tasks'][0]['mult'])!=(0,1)
            or source['maximum_host_RSS_GiB']!=80 or plan['maximum_host_RSS_GiB']!=80
            or plan['maximum_peak_GPU_MiB']!=31744 or plan['experimental_model']['solvent']!='acetonitrile'
            or plan['execution_driver']!=source['prior_root']+'/scripts/pipeline/smooth_polymer_native.py'
            or sha256_file(case/'input.xyz')!=plan['scored_geometry_sha256']):raise ValueError('Exact declared native geometry/state/resource class required')
    if alive(stage['owner']) is not False:raise ValueError('Original owner changed during continuation preparation')
    destination.mkdir(parents=True);parent=destination/'campaign_parent';parent.mkdir();newcase=parent/source['key'];newcase.mkdir();raw=destination/'raw';raw.mkdir()
    for name in ('plan.json','input.xyz'):shutil.copyfile(case/name,newcase/name)
    result=dict(source,status='native_continuation_after_foreign_residency',attempt=str(newcase),raw_parent=str(raw),
        interrupted_namespace=str(source_namespace),interrupted_namespace_sha256=sha256_file(source_namespace),
        interrupted_owner=stage['owner'],scientific_settings_changed=False,original_partial_outputs_preserved=True)
    proof=dict(status='foreign_residency_interruption_prepared_for_exact_native_continuation',source_namespace=str(source_namespace),
        source_namespace_sha256=sha256_file(source_namespace),producer_inventory=inventory,stopped_owner=stage['owner'],
        source_sha256={str(p):sha256_file(p) for p in (stage_path,guard_path,interruption_path,case/'plan.json',case/'input.xyz')},
        raw_partial_outputs_preserved=str(Path(source['raw_parent'])/source['key']),new_native_G_qualified=False,checked=time.time())
    atomic_write_json(destination/'continuation.json',proof);shutil.copyfile(__file__,destination/'preparer.py')
    atomic_write_json(destination/'namespace.json',result);reader.recheck();verify_spec(destination/'namespace.json');return result


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--recipe',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    helper=ROOT/'scripts/pipeline/observe_archived_native_namespace.py'
    import importlib.util
    definition=importlib.util.spec_from_file_location('continuation_owner_check',helper);owners=importlib.util.module_from_spec(definition);definition.loader.exec_module(owners)
    config=json.loads(a.recipe.read_text());policy=(ROOT.parent/'.config/gpu_reserve/config.env').read_text()
    hosts=re.search(r'^ALLOWED_HOSTS="([^"]+)"',policy,re.M)[1].split()
    result=prepare(Path(config['source_namespace']),Path(config['destination']),owners.alive,hosts)
    print(json.dumps(dict(status=result['status'],namespace=str(Path(config['destination'])/'namespace.json'))))
