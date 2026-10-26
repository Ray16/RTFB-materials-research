#!/usr/bin/env python
"""Retire only idle owned electronic controllers before a gate-recovery upgrade."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import importlib.util
import json
import os
from pathlib import Path
import signal
import subprocess
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json


def helper_module():
    path=ROOT/'scripts/pipeline/reconcile_unstarted_intents.py'
    spec=importlib.util.spec_from_file_location('gate_transition_owner_checks',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def transition(launch_path,out):
    helper=helper_module();launch=json.loads(launch_path.read_text())
    jobs=[j for j in launch['jobs'] if j['label'].startswith('electronic_')]
    if len(jobs)!=3 or out.exists():raise ValueError('Expected three controllers and a new audit directory')
    live=[]
    for job in jobs:
        start=helper.process_start(job['pid'])
        if start==job['process_start']:
            cmd=(Path('/proc')/str(job['pid'])/'cmdline').read_bytes()
            if b'coordinate_electronic_inventory_resilient_v2.py' not in cmd or os.getpgid(job['pid'])!=job['pid']:
                raise ValueError('Controller command or process group changed')
            live.append(job)
    if len(live)!=2:raise ValueError('Unexpected live-controller count; inspect before transition')
    paused=[];terminated=False
    try:
        for job in live:
            os.killpg(job['pid'],signal.SIGSTOP);paused.append(job['pid'])
        for job in live:
            deadline=time.monotonic()+5
            while time.monotonic()<deadline:
                if '\nState:\tT' in (Path('/proc')/str(job['pid'])/'status').read_text():break
                time.sleep(.05)
            else:raise ValueError('Owned controller could not be paused')
        cases,stopped,hosts={},[],set()
        for job in jobs:
            batch=Path(job['argv'][job['argv'].index('--batch')+1]);hosts.update(job['argv'][job['argv'].index('--hosts')+1].split(','))
            status=json.loads((batch/'dispatch.json').read_text())
            if status['owner']['pid']!=job['pid'] or any(e['status']=='running' for e in status['entries']):
                raise ValueError('Controller ownership changed or a scientific worker is running')
            manifest=json.loads((batch/'batch.json').read_text())
            cases.update({e['key']:e['campaign'] for e in manifest['entries'] if not e['existing_checked_reference']})
            for entry in status['entries']:
                if entry['status']=='stopped_attempt_preserved_review_required':
                    case=Path(entry['campaign']);handle=json.loads((case/'launch.json').read_text())
                    log=Path(handle['attempt'])/'worker.log';plan=json.loads((case/'plan.json').read_text())
                    expected=f"gpu_reserve: gpu{handle['GPU']} BUSY/reserved -> not running"
                    if (log.read_text().strip()!=expected or Path(plan['raw_root']).exists()
                            or any((case/name).exists() for name in ('experimental_energy','residency_guard','host_guard','tensor_runtime'))):
                        raise ValueError('Stopped attempt is not an unstarted gate refusal')
                    stopped.append(dict(key=entry['key'],owner=handle['owner'],log=str(log),log_sha256=helper.digest(log.read_bytes())))
        if len(stopped)!=2:raise ValueError('Unexpected gate-refusal set; inspect before transition')
        old_pids=[entry['owner']['pid'] for entry in stopped]
        with ThreadPoolExecutor(max_workers=len(hosts)) as pool:
            observations=dict(zip(sorted(hosts),pool.map(lambda host:helper.remote_check(host,cases,old_pids),sorted(hosts))))
        if any(v['matches'] for v in observations.values()):raise ValueError('A host still has a matching worker')
        for entry in stopped:
            owner=entry['owner'];actual=observations[owner['host']]['old_pid_starts'].get(str(owner['pid']))
            if actual==owner['process_start']:raise ValueError('A prior gate owner remains live')
        gate=subprocess.run([str(ROOT.parent/'bin/gpu_reserve'),'status'],capture_output=True,text=True,timeout=30,check=True)
        for entry in stopped:
            owner=entry['owner']
            if f" {owner['pid']} {owner['host']} " in gate.stdout:
                raise ValueError('A prior gate owner still holds a GPU reservation')
        out.mkdir(parents=True,exist_ok=False)
        audit=dict(status='idle_owned_controllers_checked_before_transition',checked=time.time(),
            launch_manifest=str(launch_path),launch_manifest_sha256=helper.digest(launch_path.read_bytes()),
            live_controllers=[{'pid':j['pid'],'start':j['process_start']} for j in live],
            terminal_controllers=[{'pid':j['pid'],'start':j['process_start']} for j in jobs if j not in live],
            stopped_gate_refusals=stopped,hosts=observations,gate_status_sha256=helper.digest(gate.stdout.encode()),scientific_workers_stopped=0)
        atomic_write_json(out/'audit.json',audit)
        for pid in paused:
            os.killpg(pid,signal.SIGTERM);os.killpg(pid,signal.SIGCONT)
        terminated=True;deadline=time.monotonic()+10
        while time.monotonic()<deadline and any(helper.process_start(j['pid'])==j['process_start'] for j in live):time.sleep(.1)
        if any(helper.process_start(j['pid'])==j['process_start'] for j in live):
            raise ValueError('Owned controller did not terminate; do not start replacement')
        audit.update(status='idle_controllers_retired_gate_refusals_preserved',finished=time.time())
        atomic_write_json(out/'audit.json',audit);return audit
    finally:
        if not terminated:
            for pid in paused:
                try:os.killpg(pid,signal.SIGCONT)
                except ProcessLookupError:pass


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--launch-manifest',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();r=transition(a.launch_manifest.resolve(),a.out.resolve());print(json.dumps({'status':r['status'],'gate_refusals':len(r['stopped_gate_refusals'])}))
