"""Stop only this campaign's live V100 radical controllers, preserving all output."""
import json
import os
from pathlib import Path
import signal
import socket
import time

import psutil
from redox.core.io import atomic_write_json

base=Path(__file__).resolve().parent;host=socket.gethostname().split('.')[0]
assert host in ('lambda1','lambda11')
requests=json.loads((base/'resource_routing/deferrals.json').read_text())['stop_requests']
launch=json.loads((base/f'launch_{host}.json').read_text());results=[]
for request in requests:
    if request['host']!=host:continue
    ownership=[]
    for path in (base/'workers').glob(f'{host}_gpu*.json'):
        owner=json.loads(path.read_text())
        if owner['pid']==request['queue_pid'] and owner['process_start']==request['queue_start']:ownership.append(owner)
    assert len(ownership)==1;owner=ownership[0]
    if owner.get('active_key')!=request['key']:
        results.append(dict(key=request['key'],status='owner_changed_activity_not_stopped'));continue
    worker=next(w for w in launch['workers'] if w['gpu']==owner['gpu'])
    controller=psutil.Process(worker['controller_pid'])
    start=Path(f'/proc/{controller.pid}/stat').read_text().rsplit(')',1)[1].split()[19]
    assert start==worker['process_start']
    tree=[controller]+controller.children(recursive=True)
    assert any(p.pid==owner['pid'] for p in tree)
    identities={p.pid:Path(f'/proc/{p.pid}/stat').read_text().rsplit(')',1)[1].split()[19] for p in tree}
    controller.send_signal(signal.SIGTERM)
    deadline=time.time()+35
    while time.time()<deadline:
        live=[]
        for pid,start in identities.items():
            try:
                f=Path(f'/proc/{pid}/stat').read_text().rsplit(')',1)[1].split()
                if f[19]==start and f[0]!='Z':live.append(pid)
            except FileNotFoundError:pass
        if not live:break
        time.sleep(.5)
    assert not live,live
    results.append(dict(key=request['key'],status='only_original_owned_tree_verified_terminal',process_identities=identities,
                        numerical_thresholds_changed=False,outputs_deleted=False,checked=time.time()))
atomic_write_json(base/f'resource_routing/stopped_{host}.json',dict(status='resource_migration_stop_verified',results=results))
print(results)
