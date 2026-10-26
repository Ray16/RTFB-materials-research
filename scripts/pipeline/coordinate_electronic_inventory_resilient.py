#!/usr/bin/env python
"""Use per-host reservation inventory with the bound-aware electronic controller.

The numerical producer and scientific plans are unchanged. Failed host queries
are preserved as unknown and never interpreted as free GPU capacity.
"""
import argparse
import json
import os
from pathlib import Path
import runpy
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows import resources


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--batch', type=Path, required=True)
    parser.add_argument('--hosts', required=True)
    parser.add_argument('--max-hours', type=float, default=24)
    parser.add_argument('--resume', action='store_true')
    args = parser.parse_args()
    inventory = api.load_module('per_host_resilient_electronic_inventory',
                                ROOT / 'scripts/pipeline/coordinate_native_inventory_resilient.py')
    hosts = args.hosts.split(',')
    if (not hosts or len(set(hosts)) != len(hosts) or not set(hosts) <= inventory.permitted_hosts(inventory.POLICY)
            or not 0 < args.max_hours <= 48):
        parser.error('Use distinct policy-permitted hosts and a finite runtime up to 48 hours')
    profile = json.loads((ROOT / 'environments/cutensor22/runtime.json').read_text())
    if Path(sys.executable) != Path(profile['python']) or os.environ.get('LD_LIBRARY_PATH'):
        parser.error('Use the declared Python and host loader environment')
    batch = args.batch.resolve()
    manifest = json.loads((batch / 'batch.json').read_text())
    api.verify_sources(manifest['source_sha256'])
    wrapper = ROOT / 'scripts/pipeline/coordinate_with_capacity_preflight.py'
    proof = {str(p): api.digest(p) for p in (Path(__file__).resolve(), Path(inventory.__file__), wrapper,
                ROOT / 'scripts/pipeline/dispatch_smooth_energy_capacity.py',
                ROOT / 'scripts/pipeline/coordinate_smooth_electronic_batch.py', inventory.GATE, inventory.POLICY)}
    out = batch / 'resilient_inventory'
    out.mkdir(exist_ok=True)
    owner = dict(host=os.uname().nodename.split('.')[0], pid=os.getpid(),
                 process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19])
    owner_path = out / ('owner_' + str(owner['pid']) + '_' + owner['process_start'] + '.json')
    record = dict(status='electronic_resilient_inventory_adapter_started', owner=owner, hosts=hosts,
                  started=time.time(), batch_sha256=api.digest(batch / 'batch.json'), source_sha256=proof,
                  scope='Reservation discovery only; unknown hosts are not considered available')
    atomic_write_json(owner_path, record)
    previous = resources.free_slots

    def free_slots(requested=None):
        api.verify_sources(proof)
        slots, checks = inventory.query_all(hosts if requested is None else list(requested))
        atomic_write_json(out / 'status.json', dict(owner=owner, checked=time.time(), checks=checks,
                          available_slots=slots, unverified_hosts=[c['host'] for c in checks if c['status'] != 'verified']))
        return slots

    resources.free_slots = free_slots
    try:
        sys.argv = [str(wrapper), 'coordinate', '--batch', str(batch), '--hosts', ','.join(hosts),
                    '--max-hours', str(args.max_hours)] + (['--resume'] if args.resume else [])
        runpy.run_path(str(wrapper), run_name='__main__')
    finally:
        resources.free_slots = previous
        record.update(status='electronic_resilient_inventory_adapter_stopped', finished=time.time())
        atomic_write_json(owner_path, record)


if __name__ == '__main__':
    main()
