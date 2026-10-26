#!/usr/bin/env python
"""Run one exact external-control UMA resource pilot after Nitro parity review.

The first state establishes molecule-specific resource evidence before fan-out.
Unknown launches and failed scientific attempts are preserved, not retried.
"""
import argparse
import json
import os
from pathlib import Path
import shlex
import socket
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows import resources


SCRIPT = Path(__file__).resolve()
WORKER = ROOT / 'scripts/pipeline/uma_external_control.py'
INVENTORY = ROOT / 'scripts/pipeline/coordinate_native_inventory_resilient.py'
CPU_PREFLIGHT = ROOT / 'scripts/fleet/preflight_cpu_dependencies.py'
POLICY = ROOT / 'config/uma_execution.json'
GATE = ROOT.parent / 'bin/gpu_reserve'


def parity_ready(observer):
    status_path = observer / 'status.json'
    if not status_path.exists():
        return False
    status = json.loads(status_path.read_text())
    if status.get('status') != 'independent_fixed_device_pair_passed':
        return False
    path = Path(status['review'])
    if api.digest(path) != status['review_sha256']:
        raise ValueError('Independent Nitro parity review changed')
    review = json.loads(path.read_text())
    api.verify_sources(review['source_sha256'])
    if (review.get('status') != 'independent_fixed_device_pair_passed'
            or not review.get('all_energy_and_force_checks_pass')
            or not review.get('actual_nitro_and_reference_modes_checked')
            or not review.get('same_device_checked') or not review.get('exclusive_residency_checked')
            or review.get('diagnostic_only')):
        raise ValueError('Nitro parity evidence is incomplete')
    return True


def qualified_parity_device(observer):
    if not parity_ready(observer):
        return None
    status = json.loads((observer/'status.json').read_text())
    review = json.loads(Path(status['review']).read_text())
    paths = [Path(p) for p in review['source_sha256'] if p.endswith('/nitro/status.json')]
    if len(paths) != 1:
        raise ValueError('Expected exactly one source-bound Nitro device record')
    record = json.loads(paths[0].read_text())
    if record.get('mode') != 'nitro' or not record.get('device_name') or len(record.get('capability',[])) != 2:
        raise ValueError('Qualified Nitro device identity is incomplete')
    return dict(name=record['device_name'], capability=tuple(record['capability']))


def cpu_preflight(campaign, host, policy, out):
    path = out / (host + '_preflight.json')
    if path.exists() and time.time() - path.stat().st_mtime < 600:
        cached = json.loads(path.read_text())
        if (cached.get('ok') and cached.get('python') == policy['nitro_python']
                and cached.get('host','').split('.')[0] == host and cached.get('GPU_imports_blocked_in_CPU_preflight')
                and cached.get('GPU_modules_loaded') == []):
            return path
    args = [sys.executable, str(CPU_PREFLIGHT), host, '--python', policy['nitro_python'],
            '--module', 'torch', '--module', 'fairchem', '--module', 'ase', '--module', 'rdkit', '--module', 'psutil',
            '--library-path', policy['nitro_library_path'], '--path', str(campaign / 'plan.json'),
            '--min-available-ram-mib', str(48 * 1024)]
    result = subprocess.run(args, capture_output=True, text=True, timeout=40, check=True)
    data = json.loads(result.stdout)
    if (not data.get('ok') or data.get('python') != policy['nitro_python']
            or data.get('host','').split('.')[0] != host
            or not data.get('GPU_imports_blocked_in_CPU_preflight') or data.get('GPU_modules_loaded') != []):
        raise ValueError('CPU-only dependency preflight did not pass')
    atomic_write_json(path, data)
    return path


def launch_local(campaign, observer, out, gpu, preflight):
    adapter = api.load_module('external_UMA_exact_campaign_adapter', WORKER)
    adapter.validate_campaign(campaign)
    qualified_device = qualified_parity_device(observer)
    if qualified_device is None:
        raise ValueError('Independent Nitro parity review is not passed')
    policy = json.loads(POLICY.read_text())
    checked = json.loads(preflight.read_text())
    if (not checked.get('ok') or checked.get('python') != policy['nitro_python']
            or checked.get('host','').split('.')[0] != socket.gethostname().split('.')[0]
            or not checked.get('GPU_imports_blocked_in_CPU_preflight') or checked.get('GPU_modules_loaded') != []
            or time.time() - preflight.stat().st_mtime >= 600):
        raise ValueError('Fresh CPU-only node preflight required')
    if ((campaign/'claims/neu').exists() or (campaign/'completed/neu.json').exists()
            or (campaign/'failures/neu.json').exists()
            or (campaign/'output/nitro'/adapter.GROUP/'neu').exists()):
        raise ValueError('Existing pilot attempt or claim is preserved')
    with resources.allocation_guard():
        resources.host_capacity(2, minimum_ram_gib=48)
        cpus = resources.cpu_allocation(2)
        device = resources.gpu_suitability(gpu)
        capability = tuple(map(int, device['compute_capability'].split('.')))
        if (capability != qualified_device['capability']
                or device['name'] != qualified_device['name']):
            return dict(status='unsupported_nitro_device', device=device)
        host = socket.gethostname().split('.')[0]
        if (host, gpu) not in resources.free_slots([host]):
            return dict(status='atomic_placement_unavailable', host=host, GPU=gpu)
        attempt = out / 'attempt_001'
        attempt.mkdir(parents=True, exist_ok=False)
        caps = [f'{key}=2' for key in ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS',
                                      'OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')]
        argv = ['taskset','--cpu-list',','.join(map(str,cpus)),str(GATE),'run',str(gpu),
                '--label','redox-ndi-external-UMA-pilot','--','env',*caps,
                'PYTHONPATH='+str(ROOT/'src'),'LD_LIBRARY_PATH='+policy['nitro_library_path'],
                'FAIRCHEM_CACHE_DIR='+policy['fairchem_cache'],policy['nitro_python'],
                str(ROOT/'scripts/fleet/exclusive_python_worker.py'),'--guard-out',str(attempt/'residency_guard'),
                str(ROOT/'scripts/fleet/bounded_host_python_worker.py'),'--out',str(attempt/'host_guard'),
                '--maximum-rss-gib','16','--minimum-free-gib','32',
                str(ROOT/'scripts/fleet/record_uma_gpu_memory.py'),'--out',str(attempt/'gpu_memory'),
                str(WORKER),'worker',
                '--campaign',str(campaign),'--state','neu','--require-nitro']
        with (attempt/'worker.log').open('xb') as log:
            process = subprocess.Popen(argv, cwd=ROOT, stdin=subprocess.DEVNULL, stdout=log,
                                       stderr=subprocess.STDOUT, start_new_session=True)
        stat = Path('/proc')/str(process.pid)/'stat'
        start = stat.read_text().rsplit(')',1)[1].split()[19] if stat.exists() else None
        return dict(status='gate_started_scientific_confirmation_pending', host=host, GPU=gpu,
                    owner=dict(host=host,pid=process.pid,process_start=start), device=device,
                    CPU_affinity=cpus, argv=argv, attempt=str(attempt), started=time.time())


def monitor_launched(campaign, result, status, status_path, deadline, proof):
    owners = api.load_module('external_UMA_pilot_owner_check', ROOT/'scripts/pipeline/run_polymer_spacing.py')
    attempt = Path(result['attempt'])
    handle = result['owner']
    while time.time() < deadline:
        api.verify_sources(proof)
        try:
            live = owners.owner_alive(handle)
        except (OSError,RuntimeError,subprocess.SubprocessError) as error:
            status.update(status='owned_worker_observation_unknown_preserve',error=repr(error),checked=time.time())
            atomic_write_json(status_path,status);time.sleep(15);continue
        claim = campaign/'claims/neu/owner.json'
        guard_path = attempt/'residency_guard/state.json'
        guard_owner_path = attempt/'residency_guard/owner.json'
        if claim.exists() and guard_path.exists() and guard_owner_path.exists():
            scientific = json.loads(claim.read_text())
            guard = json.loads(guard_path.read_text())
            guard_owner = json.loads(guard_owner_path.read_text())
            if (scientific['pid'] == guard_owner['pid'] and scientific['process_start'] == guard_owner['process_start']
                    and guard.get('own_residency_observed') and not guard.get('foreign_resident_PIDs')):
                try:
                    started = owners.owner_alive(scientific)
                except (OSError,RuntimeError,subprocess.SubprocessError):
                    started = False
                if started:
                    status.update(status='UMA_pilot_scientifically_started',scientific_starts=1,scientific_owner=scientific)
        status['checked'] = time.time()
        atomic_write_json(status_path,status)
        if not live:
            for _ in range(5):
                if (campaign/'completed/neu.json').exists() or (campaign/'failures/neu.json').exists():
                    break
                time.sleep(2)
            if (campaign/'completed/neu.json').exists():
                reviewer = api.load_module('external_UMA_independent_pilot_review',
                                           ROOT/'scripts/validation/review_external_uma_pilot.py')
                review_path = attempt/'independent_review.json'
                reviewed = reviewer.review(campaign,attempt,review_path)
                status.update(status=reviewed['status'],receipt=str(campaign/'completed/neu.json'),
                              independent_review=str(review_path),review_sha256=api.digest(review_path))
            elif (campaign/'failures/neu.json').exists():
                status.update(status='UMA_pilot_failed_preserved',failure=str(campaign/'failures/neu.json'))
            else:
                status.update(status='gate_ended_without_verified_scientific_result_preserve',
                              worker_log=str(attempt/'worker.log'))
            return
        time.sleep(15)
    status['status'] = 'observer_window_ended_owned_worker_not_restarted'


def coordinate(campaign, observer, out, hosts, hours):
    inventory = api.load_module('external_UMA_pilot_inventory', INVENTORY)
    if not hosts or len(set(hosts)) != len(hosts) or not set(hosts) <= inventory.permitted_hosts(inventory.POLICY):
        raise ValueError('Use distinct currently permitted hosts')
    if not 0 < hours <= 48 or out.exists():
        raise ValueError('Use a new output and a finite window up to 48 hours')
    adapter = api.load_module('external_UMA_pilot_adapter', WORKER)
    adapter.validate_campaign(campaign)
    with resources.allocation_guard():
        resources.host_capacity(1)
        os.sched_setaffinity(0, resources.cpu_allocation(1))
    out.mkdir(parents=True)
    policy = json.loads(POLICY.read_text())
    paths = (SCRIPT, WORKER, INVENTORY, CPU_PREFLIGHT, POLICY, campaign/'plan.json',
             ROOT/'scripts/fleet/exclusive_python_worker.py', ROOT/'scripts/fleet/bounded_host_python_worker.py',
             ROOT/'scripts/fleet/record_uma_gpu_memory.py',ROOT/'scripts/validation/review_external_uma_pilot.py')
    proof = {str(path):api.digest(path) for path in paths}
    owner = dict(host=socket.gethostname().split('.')[0],pid=os.getpid(),
                 process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    status = dict(status='waiting_for_independent_Nitro_parity',owner=owner,started=time.time(),
                  deadline=time.time()+hours*3600,source_sha256=proof,scientific_starts=0,
                  scope='One neutral UMA memory/resource pilot; no solution energy or redox accuracy inferred')
    status_path = out/'status.json'
    atomic_write_json(status_path,status)
    cooldown = {}
    try:
        while time.time() < status['deadline']:
            api.verify_sources(proof)
            if not parity_ready(observer):
                path = observer/'status.json'
                observed = json.loads(path.read_text()) if path.exists() else {}
                terminal = {'independent_fixed_device_pair_failed','independent_review_failed_preserved',
                            'waiter_not_live_no_reviewable_result','pair_finished_without_reviewable_result',
                            'observer_window_ended_without_reviewable_result'}
                if observed.get('status') in terminal:
                    status.update(status='Nitro_parity_not_qualified_pilot_not_started',
                                  parity_status=observed['status'])
                    return
                status.update(checked=time.time(),status='waiting_for_independent_Nitro_parity')
                atomic_write_json(status_path,status);time.sleep(15);continue
            slots, checks = inventory.query_all(hosts)
            status.update(checked=time.time(),status='waiting_for_free_suitable_Nitro_GPU',inventory=checks,available_slots=slots)
            atomic_write_json(status_path,status)
            for host,gpu in slots:
                if time.time() < cooldown.get(host,0):
                    continue
                try:
                    preflight = cpu_preflight(campaign,host,policy,out)
                except (OSError,ValueError,subprocess.SubprocessError) as error:
                    cooldown[host]=time.time()+120
                    status.setdefault('preflight_failures',[]).append(dict(host=host,error=repr(error),checked=time.time()))
                    atomic_write_json(status_path,status);continue
                intent = dict(host=host,GPU=gpu,owner=owner,checked=time.time(),outcome='requested')
                atomic_write_json(out/'dispatch_intent.json',intent)
                argv = ['env','-u','LD_LIBRARY_PATH','OMP_NUM_THREADS=1','OMP_THREAD_LIMIT=1',
                        'OPENBLAS_NUM_THREADS=1','MKL_NUM_THREADS=1','NUMEXPR_NUM_THREADS=1',
                        'PYTHONPATH='+str(ROOT/'src'),sys.executable,str(SCRIPT),'launch-local',
                        '--campaign',str(campaign),'--observer',str(observer),'--out',str(out),
                        '--gpu',str(gpu),'--preflight',str(preflight)]
                try:
                    response = subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8',host,shlex.join(argv)],
                        capture_output=True,text=True,timeout=45,check=True)
                    result = json.loads(response.stdout)
                except (OSError,ValueError,subprocess.SubprocessError) as error:
                    status.update(status='launch_observation_unverified_preserve',error=repr(error));return
                if result.get('status') == 'gate_started_scientific_confirmation_pending':
                    status.update(status='gate_started_scientific_confirmation_pending',launch=result)
                    atomic_write_json(status_path,status)
                    monitor_launched(campaign,result,status,status_path,status['deadline'],proof)
                    return
                if result.get('status') not in ('atomic_placement_unavailable','unsupported_nitro_device'):
                    status.update(status='unexpected_launch_response_preserved',launch=result);return
                atomic_write_json(out/'dispatch_intent.json',dict(intent,outcome='not_started'))
            time.sleep(15)
        status['status']='wait_window_ended_without_pilot_start'
    except BaseException as error:
        status.update(status='pilot_coordinator_failed_existing_attempt_preserved',error=repr(error))
        raise
    finally:
        status['finished']=time.time();atomic_write_json(status_path,status)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('coordinate','launch-local'))
    p.add_argument('--campaign',type=Path,required=True)
    p.add_argument('--observer',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True)
    p.add_argument('--hosts')
    p.add_argument('--max-hours',type=float,default=48)
    p.add_argument('--gpu',type=int)
    p.add_argument('--preflight',type=Path)
    a = p.parse_args()
    if a.action == 'coordinate':
        if not a.hosts: p.error('coordinate requires --hosts')
        coordinate(a.campaign.resolve(),a.observer.resolve(),a.out.resolve(),a.hosts.split(','),a.max_hours)
    else:
        if a.gpu is None or a.preflight is None: p.error('launch-local requires --gpu and --preflight')
        print(json.dumps(launch_local(a.campaign.resolve(),a.observer.resolve(),a.out.resolve(),a.gpu,a.preflight.resolve())))


if __name__=='__main__':main()
