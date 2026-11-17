#!/usr/bin/env python3
"""Give ready native thermochemistry first access to suitable free capacity.

Only admission order changes. Running electronic workers, their resource bounds,
the atomic GPU gate, and all scientific calculations remain unchanged.
"""
import importlib.util
import json
import os
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.workflows.native_resource_profiles import output_room


def native_status(path):
    if not path.exists():
        raise ValueError('The native priority dispatcher must publish its initial state first')
    status = json.loads(path.read_text())
    owner = status['owner']
    if owner['host'].split('.')[0] != os.uname().nodename.split('.')[0]:
        raise ValueError('Priority and electronic dispatchers must share a host for exact owner checks')
    process = Path('/proc')/str(int(owner['pid']))/'stat'
    try:
        fields = process.read_text().rsplit(')', 1)[1].split()
    except FileNotFoundError:
        return None
    if fields[19] != owner['process_start'] or fields[0] == 'Z' or 'finished' in status:
        return None
    return status


def native_waiting(path):
    status=native_status(path)
    return status is not None and status.get('pending',0)>0


def native_capacity_waiting(path,resources,fits):
    """Reserve capacity only for an actual pending task that fits this host."""
    state=native_status(path)
    if state is None or state.get('pending',0)<=0:return False
    admissions=state.get('pending_admissions')
    if admissions is None:admissions=[dict(maximum_host_RSS_GiB=80)]
    return any(fits(row,resources) for row in admissions)


def main():
    index = sys.argv.index('--native-priority')
    path = Path(sys.argv[index+1]).resolve()
    del sys.argv[index:index+2]
    original = ROOT/'scripts/pipeline/run_namespace_energy_queue.py'
    definition = importlib.util.spec_from_file_location('unchanged_energy_queue_with_stage_priority', original)
    module = importlib.util.module_from_spec(definition)
    definition.loader.exec_module(module)
    proof = {str(p): sha256_file(p) for p in [Path(__file__).resolve(), original,
             ROOT/'scripts/pipeline/run_checked_native_queue.py',
             ROOT/'src/redox/workflows/native_resource_profiles.py']}
    output = Path(sys.argv[sys.argv.index('--out')+1]).resolve()
    atomic_write_json(output.with_name(output.name+'_native_priority_adapter.json'),
                      dict(status='admission_order_only_running_workers_unchanged',
                           native_priority_status=str(path), source_sha256=proof))
    ordinary_fits = module.fits

    def fits(spec, resources):
        if any(sha256_file(p) != digest for p, digest in proof.items()):
            raise ValueError('Priority admission source changed')
        return (output_room(ROOT)['admit_new_work'] and ordinary_fits(spec,resources)
                and not native_capacity_waiting(path,resources,ordinary_fits))

    module.fits = fits
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--python', required=True)
    parser.add_argument('--hours', type=float, default=24)
    parser.add_argument('--resume', action='store_true')
    a = parser.parse_args()
    if not 0 < a.hours <= 48:
        parser.error('A bounded window up to 48 hours is required')
    module.coordinate(a.manifest.resolve(), a.out.resolve(), a.python, a.hours, a.resume)


if __name__ == '__main__':
    main()
