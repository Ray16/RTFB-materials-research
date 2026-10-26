#!/usr/bin/env python
"""Resume the native coordinator with per-host, failure-visible GPU discovery.

Only availability discovery is adapted. The existing native worker, physical
checks, source-bound plans and atomic GPU reservation remain unchanged.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import json
import os
from pathlib import Path
import re
import subprocess
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows import resources

GATE=ROOT.parent/'bin/gpu_reserve'
ORIGINAL=ROOT/'scripts/pipeline/coordinate_smooth_native_batch.py'
POLICY=ROOT.parent/'.config/gpu_reserve/config.env'


def permitted_hosts(path):
    text=path.read_text()
    allowed=re.search(r'^ALLOWED_HOSTS="([^"]*)"\s*$',text,re.MULTILINE)
    reserved=re.search(r'^RESERVED_HOSTS="([^"]*)"\s*$',text,re.MULTILINE)
    if not allowed or not reserved:raise ValueError('Cannot parse authoritative GPU policy')
    return set(allowed.group(1).split())-set(reserved.group(1).split())


def query_host(host,gate=GATE,run=subprocess.run):
    try:
        result=run([str(gate),'list','--hosts',host],capture_output=True,text=True,timeout=15,check=True)
        slots=[]
        for line in result.stdout.splitlines():
            fields=line.split()
            if len(fields)!=2 or fields[0]!=host or not fields[1].isdigit():
                raise ValueError('Malformed or misattributed reservation inventory')
            slots.append((host,int(fields[1])))
        if len(set(slots))!=len(slots):raise ValueError('Duplicate GPU in reservation inventory')
        return dict(host=host,status='verified',slots=slots)
    except (OSError,ValueError,subprocess.SubprocessError) as error:
        return dict(host=host,status='inventory_unverified',error=repr(error),slots=[])


def query_all(hosts,query=query_host):
    if len(hosts)!=len(set(hosts)) or any(not h or h.startswith('-') for h in hosts):
        raise ValueError('Distinct explicit hosts required')
    with ThreadPoolExecutor(max_workers=len(hosts)) as executor:
        checks=list(executor.map(query,hosts))
    return sorted(slot for check in checks for slot in check['slots']),checks


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch',type=Path,required=True);parser.add_argument('--hosts',required=True)
    parser.add_argument('--max-hours',type=float,default=24);a=parser.parse_args()
    hosts=a.hosts.split(',');batch=a.batch.resolve()
    if not 0<a.max_hours<=48:parser.error('Finite runtime of at most 48 hours required')
    if not set(hosts)<=permitted_hosts(POLICY):parser.error('Requested host is outside the current GPU policy')
    profile=json.loads((ROOT/'environments/cutensor22/runtime.json').read_text())
    if Path(os.sys.executable)!=Path(profile['python']) or os.environ.get('LD_LIBRARY_PATH'):
        parser.error('Use the declared Python and host loader environment')
    manifest=json.loads((batch/'batch.json').read_text());api.verify_sources(manifest['source_sha256'])
    out=batch/'resilient_inventory';out.mkdir(exist_ok=True)
    owner=dict(host=os.uname().nodename.split('.')[0],pid=os.getpid(),
               process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    start=dict(status='resilient_inventory_adapter_started',owner=owner,hosts=hosts,started=time.time(),
        batch_sha256=api.digest(batch/'batch.json'),adapter_sha256=api.digest(__file__),
        original_coordinator_sha256=api.digest(ORIGINAL),gate_sha256=api.digest(GATE),policy_sha256=api.digest(POLICY),
        scope='Discovery only: unresolved hosts are never treated as verified free or busy')
    owner_file=out/('owner_'+str(owner['pid'])+'_'+owner['process_start']+'.json')
    atomic_write_json(owner_file,start)
    original_free=resources.free_slots
    def resilient_free(requested=None):
        chosen=hosts if requested is None else list(requested)
        slots,checks=query_all(chosen)
        atomic_write_json(out/'status.json',dict(owner=owner,checked=time.time(),checks=checks,
            available_slots=slots,unverified_hosts=[r['host'] for r in checks if r['status']!='verified']))
        return slots
    resources.free_slots=resilient_free
    try:
        original=api.load_module('resilient_original_native_coordinator',ORIGINAL)
        original.coordinate(batch,hosts,a.max_hours,True,original.native_placement(profile),profile)
    finally:
        resources.free_slots=original_free
        start['finished']=time.time();start['status']='resilient_inventory_adapter_stopped'
        atomic_write_json(owner_file,start)


if __name__=='__main__':main()
