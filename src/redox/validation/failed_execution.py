"""Closed physical failures, distinct from resource or residency interruption."""
from pathlib import Path
from redox.validation.execution_review import host_usage_peak


def physical_failure(case,read,alive,*,status,error,host_GiB):
    case=Path(case);placements=list(case.parent.glob('placement*/launch.json'))
    if len(placements)!=1:raise ValueError('One actual failed scientific placement required')
    stage=read(case/'status.json');launch=read(placements[0]);completion=read(placements[0].with_name('completion.json'))
    guard=read(case/'residency_guard/owner.json');resident=read(case/'residency_guard/state.json')
    host=read(case/'host_guard/status.json');runtime=read(case/'tensor_runtime/status.json')
    owners=[stage['owner'],dict(host=launch['host'],pid=launch['pid'],process_start=launch['process_start']),launch['budget_observer_owner']]
    if (stage['status']!=status or stage.get('error')!=error or stage.get('native_G_qualified') is not False
            or any(alive(owner) is not False for owner in owners) or completion['returncode'] in (None,0)
            or completion['budget_stop'] is not None or not guard['guard_stopped']
            or resident['foreign_resident_PIDs'] or not resident['own_residency_observed']
            or (case/'residency_guard/external_interruption.json').exists()
            or host['maximum_RSS_GiB']!=host_GiB or host['minimum_free_GiB']!=32
            or host_usage_peak(host)>host_GiB*2**30 or host['free_RAM_bytes']<32*2**30
            or runtime['status']!='runtime_or_task_failed_preserved' or runtime['program']!=str(case.parent/'sources/worker.py')):
        raise ValueError('Closed physical failure without resource or residency interference required')
    owner=stage['owner']
    if (not owner['pid']==guard['pid']==host['pid']==runtime['pid']
            or not owner['process_start']==guard['process_start']==host['process_start']
            or not owner['host']==guard['host']==host['host']==launch['host']
            or not guard['CPU_affinity']==host['CPU_affinity']==runtime['CPU_affinity']==launch['CPU_affinity']
            or len(set(guard['CPU_affinity']))!=2):
        raise ValueError('Exact failed-worker ownership and two-CPU chain required')
    return stage
