#!/usr/bin/env python
"""Retire only owned idle electronic controllers after a placement-adapter fix.

All declared hosts are inspected before stopping controllers. Uncertain intents
are reopened only if their precise attempt and raw-output paths are absent.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import importlib.util
import json
import os
from pathlib import Path
import signal
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json


def load_helper():
    path = ROOT / 'scripts/pipeline/reconcile_unstarted_intents.py'
    spec = importlib.util.spec_from_file_location('retire_intent_checks', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def retire(launch_path, out):
    helper = load_helper()
    launch = json.loads(launch_path.read_text())
    jobs = [j for j in launch['jobs'] if j['label'].startswith('electronic_')]
    if len(jobs) != 3 or out.exists():
        raise ValueError('Expected three owned electronic controllers and a new audit directory')
    for j in jobs:
        if helper.process_start(j['pid']) != j['process_start'] or os.getpgid(j['pid']) != j['pid']:
            raise ValueError('Controller identity changed')
        cmd = (Path('/proc') / str(j['pid']) / 'cmdline').read_bytes()
        if b'coordinate_electronic_inventory_resilient.py' not in cmd:
            raise ValueError('Controller command changed')
    paused, terminated = [], False
    try:
        for j in jobs:
            os.killpg(j['pid'], signal.SIGSTOP)
            paused.append(j['pid'])
        for j in jobs:
            deadline = time.monotonic() + 5
            while time.monotonic() < deadline:
                if '\nState:\tT' in (Path('/proc') / str(j['pid']) / 'status').read_text():
                    break
                time.sleep(0.05)
            else:
                raise ValueError('Controller could not be paused safely')
        cases, originals, hosts = {}, {}, set()
        for j in jobs:
            batch = Path(j['argv'][j['argv'].index('--batch') + 1])
            hosts.update(j['argv'][j['argv'].index('--hosts') + 1].split(','))
            dispatch = json.loads((batch / 'dispatch.json').read_text())
            if dispatch['owner']['pid'] != j['pid'] or any(e['status'] == 'running' for e in dispatch['entries']):
                raise ValueError('Controller has running scientific work or changed ownership')
            for entry in dispatch['entries']:
                if entry['status'] != 'submission_ownership_unverified_preserved':
                    continue
                case = Path(entry['campaign'])
                path = case / 'dispatch_intent.json'
                data = path.read_bytes()
                intent = json.loads(data)
                plan = json.loads((case / 'plan.json').read_text())
                if (intent.get('outcome') is not None or intent['owner']['pid'] != j['pid']
                        or intent['owner']['process_start'] != j['process_start']
                        or not helper.no_attempt(case, Path(plan['raw_root']))):
                    raise ValueError('Uncertain case is not demonstrably unstarted')
                cases[entry['key']], originals[entry['key']] = str(case), data
        if len(cases) != 2:
            raise ValueError('Uncertain case count changed; inspect before retiring controllers')
        with ThreadPoolExecutor(max_workers=len(hosts)) as pool:
            observations = dict(zip(sorted(hosts), pool.map(
                lambda host: helper.remote_check(host, cases, []), sorted(hosts))))
        if any(v['matches'] for v in observations.values()):
            raise ValueError('A host still has a process referring to an uncertain case')
        out.mkdir(parents=True, exist_ok=False)
        record = dict(status='idle_controllers_checked_before_retirement', checked=time.time(),
                      launch_manifest=str(launch_path), launch_manifest_sha256=helper.digest(launch_path.read_bytes()),
                      old_controllers=[{'pid': j['pid'], 'start': j['process_start']} for j in jobs],
                      hosts=observations, cases={}, scientific_jobs_submitted=0)
        atomic_write_json(out / 'audit.json', record)
        for pid in paused:
            os.killpg(pid, signal.SIGTERM)
            os.killpg(pid, signal.SIGCONT)
        terminated = True
        deadline = time.monotonic() + 10
        while time.monotonic() < deadline and any(helper.process_start(j['pid']) == j['process_start'] for j in jobs):
            time.sleep(0.1)
        if any(helper.process_start(j['pid']) == j['process_start'] for j in jobs):
            raise ValueError('An old controller did not terminate; do not start a replacement')
        for key, name in cases.items():
            case = Path(name)
            intent_path = case / 'dispatch_intent.json'
            original = originals[key]
            plan = json.loads((case / 'plan.json').read_text())
            if intent_path.read_bytes() != original or not helper.no_attempt(case, Path(plan['raw_root'])):
                raise ValueError('Case changed during retirement; preserve for review')
            archive = out / (key + '.intent.json')
            with archive.open('xb') as stream:
                stream.write(original)
            intent = json.loads(original)
            intent.update(outcome='not_started', recovery_audit=str(out / 'audit.json'), recovered=time.time())
            atomic_write_json(intent_path, intent)
            record['cases'][key] = dict(case=str(case), original_intent_sha256=helper.digest(original),
                                        archived_intent=str(archive), replacement_intent_sha256=helper.digest(intent_path.read_bytes()))
            atomic_write_json(out / 'audit.json', record)
        record.update(status='idle_controllers_retired_two_unstarted_intents_reopened', finished=time.time())
        atomic_write_json(out / 'audit.json', record)
        return record
    finally:
        if not terminated:
            for pid in paused:
                try:
                    os.killpg(pid, signal.SIGCONT)
                except ProcessLookupError:
                    pass


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--launch-manifest', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    args = p.parse_args()
    result = retire(args.launch_manifest.resolve(), args.out.resolve())
    print(json.dumps({'status': result['status'], 'reopened': len(result['cases'])}))


if __name__ == '__main__':
    main()
