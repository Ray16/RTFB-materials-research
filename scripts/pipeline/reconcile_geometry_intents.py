#!/usr/bin/env python
"""Reopen geometry intents only after proving no scientific attempt started.

An unresolved intent cannot be placed by the existing queue. Its original bytes
are preserved before publishing not_started. Live controllers can then consume
that new state. No controller or scientific process is signaled.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import shlex
import subprocess
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api

FORBIDDEN = ('launch.json', 'dispatch_attempts', 'geometry', 'residency_guard',
             'host_guard', 'tensor_runtime')


def no_attempt(case, raw):
    return not raw.exists() and not any((case / name).exists() for name in FORBIDDEN)


def remote_check(host, cases, old_pids):
    program = r'''import json,os,pathlib,sys
request=json.load(sys.stdin); matches=[]; starts={}; errors=[]
for p in pathlib.Path('/proc').iterdir():
 if not p.name.isdigit():continue
 try:
  args=(p/'cmdline').read_bytes().replace(b'\0',b' ').decode(errors='replace')
  hits=[key for key,paths in request['cases'].items() if any(path in args for path in paths)]
  if hits:matches.append({'pid':int(p.name),'keys':hits})
  if int(p.name) in request['old_pids']:
   starts[p.name]=(p/'stat').read_text().rsplit(')',1)[1].split()[19]
 except (FileNotFoundError,ProcessLookupError):pass
 except OSError as error:errors.append({'pid':p.name,'error':repr(error)})
print(json.dumps({'host':os.uname().nodename.split('.')[0],'matches':matches,'old_pid_starts':starts,'errors':errors}))
'''
    payload = json.dumps(dict(cases=cases, old_pids=old_pids))
    result = subprocess.run(['ssh', '-o', 'BatchMode=yes', '-o', 'ConnectTimeout=5', host,
                             shlex.join(['python3', '-c', program])], input=payload,
                            capture_output=True, text=True, timeout=20, check=True)
    row = json.loads(result.stdout)
    if row.get('host') != host or row.get('errors') or row.get('matches'):
        raise ValueError('Cannot prove absence of matching processes on ' + host + ': ' + repr(row))
    return row


def collect(audit):
    cases = {}
    for row in audit['batches']:
        batch = Path(row['batch'])
        manifest = json.loads((batch / 'batch.json').read_text())
        api.verify_sources(manifest['source_sha256'])
        if api.digest(batch / 'batch.json') != row['batch_sha256'] or manifest['stage'] != 'geometry':
            raise ValueError('Audited geometry batch changed')
        dispatch = json.loads((batch / 'dispatch.json').read_text())
        held = {e['key'] for e in dispatch['entries'] if e['status'] == 'ownership_unverified_preserved'}
        for entry in manifest['entries']:
            if entry['key'] not in held:
                continue
            case = Path(entry['campaign']); raw_plan = (case / 'plan.json').read_bytes()
            plan = json.loads(raw_plan); raw = Path(plan['raw_root'])
            intent_bytes = (case / 'dispatch_intent.json').read_bytes(); intent = json.loads(intent_bytes)
            if (api.digest(case / 'plan.json') != entry['plan_sha256'] or intent.get('outcome') is not None
                    or not no_attempt(case, raw)):
                raise ValueError('Held case is not a source-bound unstarted intent')
            if entry['key'] in cases:
                raise ValueError('Duplicate case identity')
            cases[entry['key']] = dict(case=case, raw=raw, intent_bytes=intent_bytes,
                                      intent=intent, plan_sha256=api.digest(case / 'plan.json'))
    return cases


def reconcile(audit_path, out):
    if out.exists():
        raise ValueError('Use a new immutable recovery directory')
    inventory = api.load_module('geometry_intent_policy', ROOT / 'scripts/pipeline/coordinate_native_inventory_resilient.py')
    audit = json.loads(audit_path.read_text()); cases = collect(audit)
    if not cases:
        raise ValueError('No held intents in the audited batches')
    hosts = sorted(inventory.permitted_hosts(inventory.POLICY))
    owners = [v['intent']['owner'] for v in cases.values()]
    if not {o['host'] for o in owners} <= set(hosts) or not {v['intent']['host'] for v in cases.values()} <= set(hosts):
        raise ValueError('Original owner or placement host cannot be inspected under current policy')
    query = {k:[str(v['case']), str(v['raw'])] for k,v in cases.items()}
    with ThreadPoolExecutor(max_workers=len(hosts)) as pool:
        observations = dict(zip(hosts, pool.map(lambda h:remote_check(h, query, [o['pid'] for o in owners]), hosts)))
    for owner in owners:
        if observations[owner['host']]['old_pid_starts'].get(str(owner['pid'])) == owner['process_start']:
            raise ValueError('An original dispatcher is still live')
    gate = subprocess.run([str(inventory.GATE), 'status'], capture_output=True, text=True, timeout=30, check=True)
    for line in gate.stdout.splitlines():
        fields = line.split()
        if len(fields) >= 3 and any(fields[1:3] == [str(o['pid']), o['host']] for o in owners):
            raise ValueError('An original dispatcher still owns a GPU reservation')
    out.mkdir(parents=True, exist_ok=False)
    record = dict(status='all_host_absence_checked_reconciliation_in_progress', checked=time.time(),
                  audit=str(audit_path), audit_sha256=api.digest(audit_path), hosts=observations,
                  gate_status_sha256=hashlib.sha256(gate.stdout.encode()).hexdigest(),
                  source_sha256={str(p):api.digest(p) for p in (Path(__file__).resolve(), inventory.POLICY, inventory.GATE)},
                  cases={}, scientific_jobs_submitted=0)
    # Archive every original before publishing any reopened case.
    for key, value in cases.items():
        (out / (key + '.intent.json')).write_bytes(value['intent_bytes'])
    atomic_write_json(out / 'audit.json', record)
    for key, value in cases.items():
        case = value['case']; path = case / 'dispatch_intent.json'
        if (path.read_bytes() != value['intent_bytes'] or not no_attempt(case, value['raw'])
                or api.digest(case / 'plan.json') != value['plan_sha256']):
            raise ValueError('A held case changed after process verification')
        intent = dict(value['intent'], outcome='not_started', recovery_audit=str(out / 'audit.json'), recovered=time.time())
        atomic_write_json(path, intent)
        record['cases'][key] = dict(case=str(case), archived_intent=str(out / (key + '.intent.json')),
                                   original_intent_sha256=api.digest(out / (key + '.intent.json')),
                                   replacement_intent_sha256=api.digest(path), plan_sha256=value['plan_sha256'])
        atomic_write_json(out / 'audit.json', record)
    record.update(status='unstarted_geometry_intents_archived_and_reopened', reopened=len(cases), finished=time.time())
    atomic_write_json(out / 'audit.json', record)
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--audit', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    record = reconcile(args.audit.resolve(), args.out.resolve())
    print(json.dumps(dict(status=record['status'], reopened=record['reopened'])))


if __name__ == '__main__':
    main()
