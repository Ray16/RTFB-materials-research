"""Common closure, ownership and resource checks for independent QM reviews."""
import math
from pathlib import Path
import hashlib

SOURCE_BYTES=Path(__file__).read_bytes()


def freeze_source(case, source, snapshot=None):
    from redox.core.hashing import sha256_file
    source=Path(source);content=source.read_bytes() if snapshot is None else snapshot
    digest=hashlib.sha256(content).hexdigest()
    target=Path(case)/'review_sources'/digest/source.name
    if not target.exists():target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(content)
    if sha256_file(target)!=digest:raise ValueError('Frozen review source differs')
    return target


def validate_closed(stage, launch, completion, guard, residency, host, runtime, *,
                    owner_alive, expected_status, host_GiB, GPU_MiB, worker):
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
    quantities=[host.get(k) for k in ('peak_tree_RSS_bytes','self_VmHWM_bytes','free_RAM_bytes')]
    if (any(type(v) not in (int,float) or not math.isfinite(v) or v<0 for v in quantities)
            or max(quantities[:2])>host_GiB*2**30 or quantities[2]<32*2**30):
        raise ValueError('Finite host-memory resource checks required')
    owner=stage['owner']
    if (not owner['pid']==guard['pid']==host['pid']==runtime['pid']
            or not owner['process_start']==guard['process_start']==host['process_start']
            or not owner['host']==guard['host']==host['host']==launch['host']
            or not guard['CPU_affinity']==host['CPU_affinity']==runtime['CPU_affinity']==launch['CPU_affinity']
            or len(set(host['CPU_affinity']))!=2 or len(host['CPU_affinity'])!=2):
        raise ValueError('Exact ownership and disjoint two-CPU execution chain required')


def review_closure(case, read, alive, *, expected_status, host_GiB, GPU_MiB, worker):
    case=Path(case);placements=list(case.parent.glob('placement*/launch.json'))
    if len(placements)!=1:raise ValueError('One exact scientific placement required')
    stage=read(case/'status.json');launch=read(placements[0]);completion=read(placements[0].parent/'completion.json')
    guard=read(case/'residency_guard/owner.json');residency=read(case/'residency_guard/state.json')
    host=read(case/'host_guard/status.json');runtime=read(case/'tensor_runtime/status.json')
    validate_closed(stage,launch,completion,guard,residency,host,runtime,
        owner_alive=alive(stage['owner']),expected_status=expected_status,host_GiB=host_GiB,
        GPU_MiB=GPU_MiB,worker=worker)
    return stage
