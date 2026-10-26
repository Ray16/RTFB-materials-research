#!/usr/bin/env python
"""Recover only demonstrably unstarted occupied electronic dispatch intents.

The owned CPU controllers are briefly paused while their old intents are checked.
Any uncertainty fails closed. Existing scientific attempts are never replaced.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import shlex
import signal
import subprocess
import time

from redox.core.io import atomic_write_json


def digest(data):
    return hashlib.sha256(data).hexdigest()


def process_start(pid):
    p = Path('/proc') / str(pid) / 'stat'
    return p.read_text().rsplit(')', 1)[1].split()[19] if p.exists() else None


def no_attempt(case, raw):
    forbidden = ('launch.json', 'dispatch_attempts', 'experimental_energy', 'residency_guard',
                 'host_guard', 'tensor_runtime')
    return not raw.exists() and not any((case / name).exists() for name in forbidden)


def remote_check(host, cases, old_owners):
    program = r'''import json,os,pathlib,sys
data=json.load(sys.stdin); matches=[]; owners={}
for p in pathlib.Path('/proc').iterdir():
 if not p.name.isdigit(): continue
 try:
  args=(p/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')
  hits=[key for key,path in data['cases'].items() if path in args]
  if hits: matches.append({'pid':int(p.name),'keys':hits})
  if int(p.name) in data['old_pids']:
   owners[p.name]=(p/'stat').read_text().rsplit(')',1)[1].split()[19]
 except (OSError,ProcessLookupError): pass
print(json.dumps({'matches':matches,'old_pid_starts':owners}))
'''
    payload = json.dumps({'cases': cases, 'old_pids': old_owners})
    # Program is transmitted as an argument; paths and process identities go via stdin.
    result = subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=5', host,
                             shlex.join(['python3', '-c', program])], input=payload, text=True,
                            capture_output=True, timeout=15, check=True)
    return json.loads(result.stdout)


def reconcile(launch_path, out):
    launch = json.loads(launch_path.read_text())
    jobs = [j for j in launch['jobs'] if j['label'].startswith('electronic_')]
    if len(jobs) != 3 or out.exists():
        raise ValueError('Expected three owned controllers and a new audit directory')
    for j in jobs:
        if process_start(j['pid']) != j['process_start'] or os.getpgid(j['pid']) != j['pid']:
            raise ValueError('Controller PID/start identity or process group changed')
        cmd = (Path('/proc') / str(j['pid']) / 'cmdline').read_bytes()
        if b'coordinate_electronic_inventory_resilient.py' not in cmd:
            raise ValueError('Controller command changed')
    paused = []
    try:
        for j in jobs:
            os.killpg(j['pid'], signal.SIGSTOP)
            paused.append(j['pid'])
        for j in jobs:
            deadline = time.monotonic() + 5
            while time.monotonic() < deadline:
                status = (Path('/proc') / str(j['pid']) / 'status').read_text()
                if '\nState:\tT' in status:
                    break
                time.sleep(0.05)
            else:
                raise ValueError('Controller did not enter a stopped state')
        cases, originals, old_starts = {}, {}, {}
        for j in jobs:
            batch = Path(j['argv'][j['argv'].index('--batch') + 1])
            dispatch = json.loads((batch / 'dispatch.json').read_text())
            if dispatch['owner']['pid'] != j['pid'] or any(e['status'] == 'running' for e in dispatch['entries']):
                raise ValueError('Current controller ownership or scientific state changed')
            declared = json.loads((batch / 'batch.json').read_text())
            held = {e['key'] for e in dispatch['entries'] if e['status'] == 'submission_ownership_unverified_preserved'}
            for entry in declared['entries']:
                if entry['key'] not in held:
                    continue
                case = Path(entry['campaign'])
                intent_path = case / 'dispatch_intent.json'
                original = intent_path.read_bytes()
                intent = json.loads(original)
                plan = json.loads((case / 'plan.json').read_text())
                if intent.get('outcome') is not None or intent['owner']['host'] != 'lambda2' or not no_attempt(case, Path(plan['raw_root'])):
                    raise ValueError('Intent is not demonstrably unstarted')
                cases[entry['key']] = str(case)
                originals[entry['key']] = original
                old_starts[str(intent['owner']['pid'])] = intent['owner']['process_start']
        if len(cases) != 12:
            raise ValueError('Held case set changed; inspect before retrying reconciliation')
        hosts = ('lambda1', 'lambda2', 'lambda11', 'lambda13')  # original failed controllers' declared hosts
        with ThreadPoolExecutor(max_workers=len(hosts)) as pool:
            observations = dict(zip(hosts, pool.map(lambda host: remote_check(host, cases, list(map(int, old_starts))), hosts)))
        if any(v['matches'] for v in observations.values()):
            raise ValueError('A process still references a held case')
        for pid, start in old_starts.items():
            if observations['lambda2']['old_pid_starts'].get(pid) == start:
                raise ValueError('An original controller remains live')
        gate = subprocess.run(['/nfs/lambda_stor_01/homes/rzhu/bin/gpu_reserve', 'status'],
                              capture_output=True, text=True, timeout=30, check=True)
        if any((' ' + pid + ' lambda2 ') in gate.stdout for pid in old_starts):
            raise ValueError('An original controller still owns a GPU reservation')
        out.mkdir(parents=True, exist_ok=False)
        record = dict(status='ownership_checked_reconciliation_in_progress', checked=time.time(),
                      launch_manifest=str(launch_path), launch_manifest_sha256=digest(launch_path.read_bytes()),
                      controller_PIDs=paused, prior_owner_starts=old_starts, hosts=observations,
                      gate_status_sha256=digest(gate.stdout.encode()), cases={}, scientific_jobs_submitted=0)
        atomic_write_json(out / 'audit.json', record)
        for key, case_name in cases.items():
            case = Path(case_name)
            intent_path = case / 'dispatch_intent.json'
            original = originals[key]
            plan = json.loads((case / 'plan.json').read_text())
            if intent_path.read_bytes() != original or not no_attempt(case, Path(plan['raw_root'])):
                raise ValueError('A held case changed during reconciliation')
            backup = out / (key + '.intent.json')
            with backup.open('xb') as stream:
                stream.write(original)
            intent = json.loads(original)
            intent.update(outcome='not_started', recovery_audit=str(out / 'audit.json'), recovered=time.time())
            atomic_write_json(intent_path, intent)
            record['cases'][key] = dict(case=str(case), original_intent_sha256=digest(original),
                                        archived_intent=str(backup), replacement_intent_sha256=digest(intent_path.read_bytes()))
            atomic_write_json(out / 'audit.json', record)
        record.update(status='unstarted_intents_archived_and_reopened', finished=time.time(), reopened=len(record['cases']))
        atomic_write_json(out / 'audit.json', record)
        return record
    finally:
        for pid in paused:
            try:
                os.killpg(pid, signal.SIGCONT)
            except ProcessLookupError:
                pass


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--launch-manifest', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    result = reconcile(a.launch_manifest.resolve(), a.out.resolve())
    print(json.dumps({'status': result['status'], 'reopened': result['reopened']}))


if __name__ == '__main__':
    main()
