"""Common closure, ownership and resource checks for independent QM reviews."""
import math
from pathlib import Path
import hashlib

SOURCE_BYTES=Path(__file__).read_bytes()


def host_usage_peak(host, expected_profile=None):
    """Review the explicitly declared metric, retaining native/RSS admission."""
    profile=host.get('resource_profile')
    if expected_profile is not None and profile!=expected_profile:
        raise ValueError('Host-memory profile differs from the predeclared calculation')
    names=['peak_tree_RSS_bytes','self_VmHWM_bytes']
    if profile=='corrected_hybrid_energy_pss':
        if (expected_profile!=profile or host.get('maximum_RSS_GiB') not in (160,240)
                or host.get('memory_accounting')!='max_self_HWM_tree_PSS'):
            raise ValueError('Explicit supported proportional-memory energy profile required')
        names+=['peak_tree_PSS_bytes','peak_guard_usage_bytes']
    elif profile not in (None,'native','corrected_hybrid_energy'):
        raise ValueError('Unknown host-memory accounting profile')
    if any(type(host.get(k)) not in (int,float) or not math.isfinite(host[k]) or host[k]<0 for k in names):
        raise ValueError('Finite nonnegative host-memory observations required')
    if profile=='corrected_hybrid_energy_pss':
        peak=max(host['self_VmHWM_bytes'],host['peak_tree_PSS_bytes'])
        if host['peak_tree_PSS_bytes']<=0 or host['peak_guard_usage_bytes']!=peak:
            raise ValueError('Proportional guard peak does not reconstruct')
        return peak
    return max(host[k] for k in names)


def freeze_source(case, source, snapshot=None):
    from redox.core.hashing import sha256_file
    source=Path(source);content=source.read_bytes() if snapshot is None else snapshot
    digest=hashlib.sha256(content).hexdigest()
    target=Path(case)/'review_sources'/digest/source.name
    if not target.exists():target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(content)
    if sha256_file(target)!=digest:raise ValueError('Frozen review source differs')
    return target


def validate_closed(stage, launch, completion, guard, residency, host, runtime, *,
                    owner_alive, expected_status, host_GiB, GPU_MiB, worker, expected_memory_profile=None):
    peak=stage.get('observed_peak_GPU_MiB')
    if (stage.get('status')!=expected_status or stage.get('native_G_qualified') is not False
            or owner_alive is not False or completion.get('returncode')!=0
            or completion.get('budget_stop') is not None or not guard.get('guard_stopped')
            or not residency.get('own_residency_observed') or residency.get('foreign_resident_PIDs')
            or host.get('status')!='effective_memory_guard_finished_profile_review_pending'
            or runtime.get('status')!='task_returned_under_checked_cutensor_runtime'
            or host.get('maximum_RSS_GiB')!=host_GiB or host.get('minimum_free_GiB')!=32
            or type(peak) not in (int,float) or not math.isfinite(peak) or not 0<peak<=GPU_MiB
            or runtime.get('program')!=str(worker)):
        raise ValueError('Normal closed execution and declared resource limits required')
    quantities=[host.get('free_RAM_bytes')]
    if (any(type(v) not in (int,float) or not math.isfinite(v) or v<0 for v in quantities)
            or host_usage_peak(host,expected_memory_profile)>host_GiB*2**30 or quantities[0]<32*2**30):
        raise ValueError('Finite host-memory resource checks required')
    owner=stage['owner']
    if (not owner['pid']==guard['pid']==host['pid']==runtime['pid']
            or not owner['process_start']==guard['process_start']==host['process_start']
            or not owner['host']==guard['host']==host['host']==launch['host']
            or not guard['CPU_affinity']==host['CPU_affinity']==runtime['CPU_affinity']==launch['CPU_affinity']
            or len(set(host['CPU_affinity']))!=2 or len(host['CPU_affinity'])!=2):
        raise ValueError('Exact ownership and disjoint two-CPU execution chain required')


def review_closure(case, read, alive, *, expected_status, host_GiB, GPU_MiB, worker, expected_memory_profile=None):
    case=Path(case);placements=list(case.parent.glob('placement*/launch.json'))
    external=case/'residency_guard/external_interruption.json'
    if external.exists():
        read(external)
        raise ValueError('Independent watchdog recorded foreign residency; preserve this unqualified attempt')
    if len(placements)!=1:raise ValueError('One exact scientific placement required')
    stage=read(case/'status.json');launch=read(placements[0]);completion=read(placements[0].parent/'completion.json')
    guard=read(case/'residency_guard/owner.json');residency=read(case/'residency_guard/state.json')
    host=read(case/'host_guard/status.json');runtime=read(case/'tensor_runtime/status.json')
    validate_closed(stage,launch,completion,guard,residency,host,runtime,
        owner_alive=alive(stage['owner']),expected_status=expected_status,host_GiB=host_GiB,
        GPU_MiB=GPU_MiB,worker=worker,expected_memory_profile=expected_memory_profile)
    return stage
