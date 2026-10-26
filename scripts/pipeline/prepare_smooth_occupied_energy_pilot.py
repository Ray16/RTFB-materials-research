#!/usr/bin/env python
"""Fresh electronic/RAM pilot for an independently checked occupied graph."""
import argparse
import json
import os
from pathlib import Path
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--geometry-review',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    from redox.workflows.resources import allocation_guard,cpu_allocation,host_capacity
    with allocation_guard():host_capacity(1);os.sched_setaffinity(0,cpu_allocation(1))
    review_path=a.geometry_review.resolve();review=json.loads(review_path.read_text());geometry=Path(review['campaign'])
    old=api.read_plan(geometry);task=old['tasks'][0]
    if (task['candidate_id']!='exp_phthalimide_control' or not task['occupied_sites']
            or task['atoms'] not in (105,119)):
        raise ValueError('Pilot is limited to existing occupied neutral phthalimide controls')
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);case=out/task['key']
    original=api.load_module('prospective_occupied_exact_energy_geometry',ROOT/'scripts/pipeline/smooth_polymer_energy.py')
    original.prepare(geometry,review_path,case);plan=api.read_plan(case)
    plan.update(execution_driver=str(ROOT/'scripts/pipeline/smooth_polymer_energy_bounded.py'),
        maximum_host_RSS_GiB=256,prospective_resource_pilot=True,
        resource_scope='First occupied-graph hybrid/RAM measurement; bounded 256-GiB ceiling with 32-GiB headroom, not a transferred resource certificate')
    plan['producer_sha256'][str(Path(__file__).resolve())]=api.digest(__file__)
    atomic_write_json(case/'plan.json',plan)
    sources=[Path(__file__).resolve(),ROOT/'scripts/pipeline/smooth_polymer_energy.py',
        ROOT/'scripts/pipeline/smooth_polymer_energy_bounded.py',ROOT/'scripts/pipeline/dispatch_smooth_energy_safe.py',
        ROOT/'scripts/pipeline/dispatch_smooth_energy_pilot.py',ROOT/'scripts/pipeline/coordinate_geometry_batch.py',
        ROOT/'scripts/fleet/preflight_cpu_dependencies.py',ROOT/'environments/cutensor22/runtime.json']
    atomic_write_json(out/'batch.json',dict(status='occupied_graph_electronic_resource_pilot_prepared',stage='electronic',
        entries=[dict(key=task['key'],campaign=str(case),plan_sha256=api.digest(case/'plan.json'))],
        worker=str(ROOT/'scripts/pipeline/smooth_polymer_energy_bounded.py'),maximum_RSS_GiB=256,
        minimum_free_RAM_GiB=32,minimum_GPU_MiB=32768,primary_candidate_admission=False,
        source_sha256={str(p):api.digest(p) for p in sources},created=time.time(),
        resource_ceiling_prospective=True,quantum_values_transferred=False))
    print(json.dumps(dict(prepared=task['key'],atoms=task['atoms'],host_RSS_ceiling_GiB=256)))


if __name__=='__main__':main()
