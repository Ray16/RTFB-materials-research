#!/usr/bin/env python3
"""Dispatch prepared native states and independently review completed records.

The growing preparation stream is restricted to a frozen electronic cohort.
Every scientific worker retains the original inner guard and numerical gates.
Optional matched-class profiles add a prospectively tested stricter outer guard.
Existing externally owned namespaces are observed separately, never relaunched.
"""
import argparse
import fcntl
import importlib.util
import json
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import time

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.namespace_sources import read_checked_review
from redox.workflows.archived_native_namespace import verify_spec
from redox.workflows.derived_native_namespace import namespace_command
from redox.workflows.host_local_review import command_for_host
from redox.workflows.submission_refusal import release_prelaunch_refusal
from redox.workflows.native_resource_profiles import MatchedNativeProfiles,first_fitting,output_room


def load(path, name):
    definition = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(definition)
    definition.loader.exec_module(module)
    return module


def ready_inputs(preparation, allowed, bindings):
    """Reject changing identities and duplicates before any GPU submission."""
    rows, keys, namespaces = [], set(), set()
    for row in preparation['entries']:
        if row['status'] not in ('native_inputs_prepared', 'matching_native_namespace_reused'):
            continue
        if row['energy_namespace'] not in allowed or row['key'] in keys:
            raise ValueError('Preparation introduced an undeclared or duplicate electronic task')
        path = row['native_namespace']
        digest = sha256_file(path)
        if path in namespaces or digest != row['native_namespace_sha256'] or bindings.get(path, digest) != digest:
            raise ValueError('Prepared native namespace identity changed or duplicated')
        spec = verify_spec(Path(path))
        if spec['key'] != row['key'] or spec['maximum_host_RSS_GiB'] != 80:
            raise ValueError('Exact prepared native task and unchanged 80-GiB class required')
        keys.add(row['key'])
        namespaces.add(path)
        bindings[path] = digest
        rows.append((row, spec))
    return rows


def inspect(spec, alive):
    case = Path(spec['attempt'])
    stage = case/'experimental_native/status.json'
    if stage.exists():
        raw = json.loads(stage.read_text())
        if 'finished' in raw:
            return ('finished_review_pending' if raw['status'] ==
                    'experimental_native_finished_independent_review_pending' else 'failed_preserved'), raw
        try:
            live = alive(raw['owner'])
        except (OSError, ValueError, subprocess.SubprocessError):
            live = None
        return ('running' if live is True else 'stopped_attempt_preserved' if live is False
                else 'ownership_unverified_preserved'), raw
    if (case/'residency_guard/owner.json').exists():
        return 'guard_present_start_pending', json.loads((case/'residency_guard/owner.json').read_text())
    claim = case/'native_dispatch_claim'
    if claim.exists():
        return 'submission_claim_preserved', json.loads((claim/'owner.json').read_text())
    return 'pending', {}


def physical_review(path, spec, output, alive):
    """Wait for closed scientific ownership, then run the original reviewer."""
    case = Path(spec['attempt'])
    state = json.loads((case/'experimental_native/status.json').read_text())
    if alive(state['owner']):
        return None
    guard = case/'residency_guard/owner.json'
    if not guard.exists() or not json.loads(guard.read_text()).get('guard_stopped'):
        return None
    report_path = case/'independent_native_review.json'
    if not report_path.exists():
        # An interrupted observation can leave the source-host reviewer alive.
        # Its existing log is a handoff record, never permission to duplicate it.
        log_path = output/(spec['key']+'_review.log')
        if log_path.exists():
            return None
        profile = json.loads((Path(spec['archive'])/'environments/cutensor22/runtime.json').read_text())
        execution=case/'effective_execution.json'
        if execution.exists():
            from redox.workflows.effective_native_namespace import verify,mount
            wrapped=json.loads(execution.read_text());profile_path=Path(wrapped['profile'])
            resource_profile,checked_spec=verify(profile_path)
            if checked_spec!=spec or sha256_file(profile_path)!=wrapped['profile_sha256']:
                raise ValueError('Wrapped execution has a different native namespace/profile')
            snapshot=profile_path.parent/'qc_sources'/str(time.time_ns())/'review_effective_native.py'
            snapshot.parent.mkdir(parents=True)
            shutil.copyfile(ROOT/'scripts/validation/review_effective_native.py',snapshot)
            command=mount(resource_profile,spec)+['env','PYTHONPATH='+spec['prior_root']+'/src',
                profile['python'],str(snapshot),'--campaign',spec['original_campaign'],
                '--profile',str(profile_path),'--out',spec['original_campaign']+'/independent_native_review.json']
        else:
            command = namespace_command(spec)
            command[0] = str(ROOT.parent/'redox_polymer_continuations/runtime/bin/bwrap')
            command += ['env', 'PYTHONPATH='+spec['prior_root']+'/src', profile['python'],
                        spec['prior_root']+'/scripts/validation/review_smooth_polymer_native.py',
                        '--campaign', spec['original_campaign'], '--out',
                        spec['original_campaign']+'/independent_native_review.json']
        command = command_for_host(command, state['owner']['host'])
        with log_path.open('x') as log:
            completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, timeout=300)
        if completed.returncode:
            raise ValueError('Independent native review failed; original attempt preserved')
    checked,binding=read_checked_review(path,'native')
    execution=case/'effective_execution.json'
    if execution.exists():
        from redox.workflows.effective_native_namespace import review
        profile_path=Path(json.loads(execution.read_text())['profile'])
        resource_path=case/'effective_resource_review.json'
        if not resource_path.exists():review(profile_path,resource_path)
        resource=json.loads(resource_path.read_text())
        if (resource['status']!='effective_native_resource_profile_checked'
                or resource['normal_native_binding']!=binding
                or not resource['tested_run_profile_qualified']):
            raise ValueError('Native run lacks its matching checked effective resource result')
        binding=dict(binding,effective_resource_review=str(resource_path),
                     effective_resource_review_sha256=sha256_file(resource_path),
                     effective_bound_GiB=resource['effective_bound_GiB'])
    return checked,binding


def submit(path, spec, host, gpu, output, interpreter, effective_profile=None):
    """Preflight on the target node, then detach its existing guarded launcher."""
    directory = output/spec['key']/str(time.time_ns())
    directory.mkdir(parents=True, exist_ok=False)
    caps = ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'OPENBLAS_NUM_THREADS',
            'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'XTB_THREADS')
    prefix = ['env', *(name+'=1' for name in caps), 'PYTHONPATH='+str(ROOT/'src'),
              'PATH='+str(ROOT.parent/'redox_polymer_continuations/runtime/bin')+':/usr/local/bin:/usr/bin:/bin',
              interpreter, '-m', ('redox.workflows.derived_native_namespace' if effective_profile is None
                                 else 'redox.workflows.effective_native_namespace')]
    selected=['--spec',str(path)] if effective_profile is None else ['--profile',str(effective_profile)]
    ssh = ['ssh', '-S', 'none', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=8', host]
    probe = subprocess.run(ssh+[shlex.join(prefix+['probe',*selected])],
                           capture_output=True, text=True, timeout=240)
    (directory/'preflight_stdout.txt').write_text(probe.stdout)
    (directory/'preflight_stderr.txt').write_text(probe.stderr)
    if probe.returncode:
        return dict(status='preflight_failed_no_launch', host=host, returncode=probe.returncode)
    checked = json.loads(probe.stdout)
    expected=('derived_native_original_paths_and_CPU_dependencies_checked' if effective_profile is None
              else 'effective_native_node_sources_checked')
    key='namespace_sha256' if effective_profile is None else 'profile_sha256'
    if (checked['host'] != host or checked[key] != sha256_file(path if effective_profile is None else effective_profile)
            or checked['status'] != expected):
        raise ValueError('Node preflight has another identity or incomplete dependency check')
    argv = prefix+['launch',*selected,'--gpu',str(gpu),'--out',str(directory/'placement')]
    code = ('import json,os,pathlib,subprocess; '
            'pathlib.Path('+repr(str(directory))+').mkdir(parents=True,exist_ok=True); '
            'f=pathlib.Path('+repr(str(directory/'controller.log'))+').open("x"); '
            'p=subprocess.Popen('+repr(argv)+',stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True); '
            's=pathlib.Path(f"/proc/{p.pid}/stat").read_text().rsplit(")",1)[1].split()[19]; '
            'print(json.dumps(dict(host=os.uname().nodename,pid=p.pid,process_start=s)))')
    response = subprocess.run(ssh+[shlex.join(['/usr/bin/python3', '-c', code])],
                              capture_output=True, text=True, timeout=30)
    (directory/'spawn_stdout.txt').write_text(response.stdout)
    (directory/'spawn_stderr.txt').write_text(response.stderr)
    response.check_returncode()
    receipt = dict(status='controller_spawned_scientific_start_unverified', owner=json.loads(response.stdout),
                   argv=argv, host=host, GPU=gpu, namespace=str(path), namespace_sha256=sha256_file(path))
    atomic_write_json(directory/'controller.json', receipt)
    return receipt


def release_gate_refusal(spec, output, alive):
    """Retry an explicit atomic-gate refusal only after both owners are gone."""
    from redox.validation.namespace_sources import NamespaceSources
    case = Path(spec['attempt'])
    claim = case/'native_dispatch_claim'
    if not (claim/'owner.json').exists() or not (claim/'launch.json').exists():
        return False
    owner = json.loads((claim/'owner.json').read_text())
    receipt = json.loads((claim/'launch.json').read_text())
    if owner.get('queue') != str(output) or receipt.get('status') != 'controller_spawned_scientific_start_unverified':
        return False
    placement = Path(receipt['argv'][receipt['argv'].index('--out')+1])
    if not (placement/'completion.json').exists():
        return False
    terminal = json.loads((placement/'completion.json').read_text())
    refusal = f'gpu_reserve: gpu{terminal["GPU"]} BUSY/reserved -> not running'
    if (terminal.get('returncode') in (None, 0) or refusal not in (placement/'command.log').read_text().splitlines()
            or (case/'residency_guard').exists() or (case/'experimental_native').exists()
            or NamespaceSources(spec).path(spec['original_raw_root']).exists()):
        return False
    stopped = [receipt['owner'], dict(host=terminal['host'], pid=terminal['pid'],
                                     process_start=terminal['process_start'])]
    try:
        if any(alive(p) for p in stopped):
            return False
    except (OSError, ValueError, subprocess.SubprocessError):
        return False
    archive = output/'gate_contention_history'/spec['key']/str(time.time_ns())
    shutil.copytree(claim, archive)
    atomic_write_json(archive/'resolution.json', dict(status='gate_refused_no_scientific_attempt',
                      stopped_owners=stopped, placement=str(placement), checked=time.time(),
                      inference_attempt_consumed=False))
    for p in claim.iterdir():
        p.unlink()
    claim.rmdir()
    return True


def run(cohort, preparation_path, output, interpreter, hours, resume, recheck_hosts=(),effective_reference=None):
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    helpers = load(ROOT/'scripts/pipeline/run_namespace_energy_queue.py', 'checked_native_capacity')
    owners = load(ROOT/'scripts/pipeline/observe_archived_native_namespace.py', 'checked_native_owners')
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4)
        os.sched_setaffinity(0, cpu_allocation(1))
    config = json.loads(cohort.read_text())
    allowed = set(config['specs'])
    if not allowed or len(allowed) != len(config['specs']):
        raise ValueError('Distinct frozen electronic cohort required')
    # Discover through system policy; the cohort limits work, never host policy.
    import re
    policy = ROOT.parent/'.config/gpu_reserve/config.env'
    hosts = re.search(r'^ALLOWED_HOSTS="([^"]+)"', policy.read_text(), re.M)[1].split()
    output.mkdir(parents=True, exist_ok=True)
    frozen = {str(p): sha256_file(p) for p in [Path(__file__).resolve(), cohort, policy,
        ROOT/'scripts/pipeline/run_namespace_energy_queue.py',
        ROOT/'scripts/pipeline/observe_archived_native_namespace.py',
        ROOT/'src/redox/workflows/derived_native_namespace.py',
        ROOT/'src/redox/workflows/submission_refusal.py',
        ROOT/'src/redox/workflows/native_resource_profiles.py',
        ROOT/'src/redox/workflows/host_local_review.py', ROOT/'src/redox/validation/namespace_sources.py']}
    profiles=None if effective_reference is None else MatchedNativeProfiles(effective_reference)
    if profiles is not None:
        frozen.update(profiles.source_sha256)
        for name in ('src/redox/workflows/native_resource_profiles.py','src/redox/workflows/effective_native_namespace.py',
                     'scripts/validation/review_effective_native.py'):
            frozen[str(ROOT/name)]=sha256_file(ROOT/name)
    with (output/'queue.lock').open('a+') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        previous = json.loads((output/'status.json').read_text()) if (output/'status.json').exists() else {}
        if previous:
            if not resume or owners.alive(previous['owner']):
                raise ValueError('Previous dispatcher must be proven stopped before explicit resume')
            history = output/'history'
            history.mkdir(exist_ok=True)
            shutil.copyfile(output/'status.json', history/(str(time.time_ns())+'.json'))
        owner = dict(host=os.uname().nodename, pid=os.getpid(),
                     process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19])
        record = dict(status='native_conveyor_running', owner=owner, started=time.time(),
                      deadline=time.time()+hours*3600, source_sha256=frozen, declared_energy_states=len(allowed),
                      thermochemistry_inferred_from_energy=False,effective_memory_profile_used=False,
                      effective_reference=None if effective_reference is None else str(effective_reference))
        bindings = previous.get('namespaces_sha256', {})
        done = {row['key']: row for row in previous.get('entries', [])
                if row['status'] in ('native_independently_checked', 'failed_preserved', 'independent_review_failed_preserved')}
        stopped_hosts, accounted, failures = set(previous.get('stopped_hosts', [])), set(), {}
        if recheck_hosts and (not resume or not set(recheck_hosts) <= stopped_hosts or not set(recheck_hosts) <= set(hosts)):
            raise ValueError('Only explicitly resumed, previously stopped policy hosts can be rechecked')
        stopped_hosts.difference_update(recheck_hosts)
        record['explicit_hosts_to_recheck'] = list(recheck_hosts)
        try:
            while time.time() < record['deadline']:
                if any(sha256_file(p) != digest for p, digest in {**frozen, **bindings}.items()):
                    raise ValueError('Native dispatcher source or admitted namespace changed')
                preparation = json.loads(preparation_path.read_text())
                if preparation['source_sha256'][str(cohort)] != frozen[str(cohort)]:
                    raise ValueError('Preparation stream belongs to another frozen cohort')
                ready = ready_inputs(preparation, allowed, bindings)
                pending, rows, unsettled, newly_terminal = [], [], set(), []
                pending_profiles={}
                for entry, spec in ready:
                    key, path = spec['key'], Path(entry['native_namespace'])
                    if entry['status'] == 'matching_native_namespace_reused':
                        rows.append(dict(key=key, spec=str(path), status='external_namespace_preserved_separate_collector'))
                        continue
                    if key in done:
                        rows.append(done[key])
                        continue
                    release_gate_refusal(spec, output, owners.alive)
                    release_prelaunch_refusal(spec, output, owners.alive, hosts, 'native', frozen)
                    state, raw = inspect(spec, owners.alive)
                    row = dict(key=key, spec=str(path), status=state)
                    if state == 'finished_review_pending':
                        try:
                            checked = physical_review(path, spec, output, owners.alive)
                            if checked is not None:
                                report, binding = checked
                                row.update(status='native_independently_checked', binding=binding,
                                           G_by_cutoff_eV=report['G_by_cutoff_eV'])
                                done[key] = row
                        except (OSError, subprocess.SubprocessError) as error:
                            row.update(status='review_observation_unverified_preserved', error=repr(error))
                        except ValueError as error:
                            row.update(status='independent_review_failed_preserved', error=repr(error))
                            done[key] = row
                    elif state == 'pending':
                        pending.append((path, spec))
                        if profiles is not None:
                            matched=profiles.for_pending(path)
                            if matched is not None:pending_profiles[str(path)]=matched
                    elif state == 'failed_preserved':
                        done[key] = row
                    if state in ('guard_present_start_pending', 'submission_claim_preserved'):
                        if 'host' in raw and not (Path(spec['attempt'])/'host_guard/status.json').exists():
                            unsettled.add(raw['host'].split('.')[0])
                    if key not in accounted and state in ('finished_review_pending', 'failed_preserved'):
                        newly_terminal.append((raw['finished'], key, state, raw['owner']['host'].split('.')[0]))
                    rows.append(row)
                for _, key, state, host in sorted(newly_terminal):
                    failures[host] = 0 if state == 'finished_review_pending' else failures.get(host, 0)+1
                    if failures[host] >= 3:
                        stopped_hosts.add(host)
                    accounted.add(key)
                record.update(updated=time.time(), entries=rows, namespaces_sha256=bindings,
                              prepared_states=len(ready), pending=len(pending), independently_checked=sum(
                                  r['status'] == 'native_independently_checked' for r in rows),
                              stopped_hosts=sorted(stopped_hosts), consecutive_scientific_host_failures=failures)
                record['pending_admissions']=[dict(namespace=str(p),maximum_host_RSS_GiB=(80 if str(p) not in pending_profiles
                    else profiles.bound)) for p,_ in pending]
                record['prospective_profiles']={} if profiles is None else profiles.profiles
                record['shared_storage_admission']=output_room(ROOT)
                atomic_write_json(output/'status.json', record)
                if len(ready) == len(allowed) and all(r['status'] in (
                    'native_independently_checked', 'failed_preserved', 'independent_review_failed_preserved',
                    'external_namespace_preserved_separate_collector') for r in rows):
                    record['status'] = 'native_conveyor_terminal_outcomes_separate'
                    return
                slots, errors = helpers.discover_slots(hosts)
                record['availability_observation_errors'] = errors
                atomic_write_json(output/'status.json', record)
                for host, gpu in slots:
                    record['shared_storage_admission']=output_room(ROOT)
                    if not record['shared_storage_admission']['admit_new_work']:break
                    if not pending:
                        break
                    if host in stopped_hosts or host in unsettled:
                        continue
                    try:
                        capacity = helpers.capacity(host, gpu, interpreter)
                    except (OSError, ValueError, subprocess.SubprocessError):
                        continue
                    index,effective_profile=first_fitting(pending,capacity,pending_profiles,helpers.fits)
                    if index is None:
                        continue
                    path,spec=pending[index]
                    claim = Path(spec['attempt'])/'native_dispatch_claim'
                    try:
                        claim.mkdir()
                    except FileExistsError:
                        pending.pop(index)
                        continue
                    atomic_write_json(claim/'owner.json', dict(owner=owner, host=host, GPU=gpu,
                        queue=str(output), namespace=str(path), namespace_sha256=bindings[str(path)],
                        resource_preflight=capacity, claimed=time.time()))
                    try:
                        receipt = submit(path,spec,host,gpu,output/'placements',interpreter,effective_profile)
                        if effective_profile is not None:record['effective_memory_profile_used']=True
                        atomic_write_json(claim/'launch.json', receipt)
                        if receipt['status'] == 'preflight_failed_no_launch':
                            stopped_hosts.add(host)
                    except (OSError, ValueError, subprocess.SubprocessError) as error:
                        atomic_write_json(claim/'uncertain_submission.json', dict(error=repr(error), preserved=True))
                        stopped_hosts.add(host)
                    pending.pop(index)
                    unsettled.add(host)
                time.sleep(20)
            record['status'] = 'native_conveyor_window_ended_workers_untouched'
        except BaseException as error:
            record.update(status='native_conveyor_stopped_workers_preserved', error=repr(error))
            raise
        finally:
            record['finished'] = time.time()
            atomic_write_json(output/'status.json', record)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cohort', type=Path, required=True)
    p.add_argument('--preparation', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--python', required=True)
    p.add_argument('--hours', type=float, default=24)
    p.add_argument('--resume', action='store_true')
    p.add_argument('--recheck-host', action='append', default=[])
    p.add_argument('--effective-reference',type=Path,help='Checked resource review for prospective matched-class outer guards')
    a = p.parse_args()
    if not 0 < a.hours <= 48:
        p.error('A bounded window up to 48 hours is required')
    run(a.cohort.resolve(),a.preparation.resolve(),a.out.resolve(),a.python,a.hours,a.resume,a.recheck_host,
        None if a.effective_reference is None else a.effective_reference.resolve())
