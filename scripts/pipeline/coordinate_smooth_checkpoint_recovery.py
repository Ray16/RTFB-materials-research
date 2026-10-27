#!/usr/bin/env python
"""Apply the bounded geometry queue to the explicit checkpoint-retry worker."""
import json
import os
from pathlib import Path
import sys

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api
from redox.workflows import resources


def main():
    os.chdir(ROOT);batch=Path(sys.argv[sys.argv.index('--batch')+1]).resolve()
    profile=json.loads((ROOT/'environments/cutensor22/runtime.json').read_text())
    manifest=json.loads((batch/'batch.json').read_text())
    if (Path(sys.executable)!=Path(profile['python']) or os.environ.get('LD_LIBRARY_PATH')
            or manifest['stage']!='geometry' or manifest['worker']!=str(ROOT/'scripts/pipeline/run_smooth_checkpoint_recovery.py')
            or manifest['primary_candidate_admission']):raise ValueError('Exact checkpoint-recovery worker/runtime required')
    api.verify_sources(manifest['source_sha256'])
    inventory=api.load_module('checkpoint_recovery_inventory',ROOT/'scripts/pipeline/coordinate_native_inventory_resilient.py')
    hosts=sys.argv[sys.argv.index('--hosts')+1].split(',') if '--hosts' in sys.argv else [os.uname().nodename.split('.')[0]]
    if len(hosts)!=len(set(hosts)) or not set(hosts)<=inventory.permitted_hosts(inventory.POLICY):raise ValueError('Current-policy host set required')
    base=api.load_module('checkpoint_recovery_geometry_inventory',ROOT/'scripts/pipeline/coordinate_smooth_geometry_resilient.py')
    bounded=api.load_module('checkpoint_recovery_bounded_placement',ROOT/'scripts/pipeline/coordinate_smooth_geometry_bounded.py')
    queue=api.load_module('checkpoint_recovery_geometry_queue',ROOT/'scripts/pipeline/coordinate_geometry_batch.py')
    bounded.adapt(queue,profile,base.adapt);queue.SCRIPT=Path(__file__).resolve()
    queue.coordinate.__globals__['SCRIPT']=queue.SCRIPT
    proof={str(p):api.digest(p) for p in (Path(__file__).resolve(),Path(base.__file__),Path(bounded.__file__),inventory.GATE,inventory.POLICY)}
    original=resources.free_slots
    def free_slots(requested=None):
        api.verify_sources(proof)
        slots,_=inventory.query_all(hosts if requested is None else list(requested));return slots
    resources.free_slots=free_slots
    try:queue.main()
    finally:resources.free_slots=original


if __name__=='__main__':main()
