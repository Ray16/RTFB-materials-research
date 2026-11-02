#!/usr/bin/env python3
"""Select all declared sampling starting geometries without importing old E/G."""
import argparse
import json
import os
from pathlib import Path

from redox.core.io import atomic_write_json
from redox.workflows.sampling_inventory import inventory


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('archive','census','energy-roster','out'):p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--reference-root',type=Path,action='append',default=[])
    a=p.parse_args()
    from redox.workflows.resources import allocation_guard,host_capacity,cpu_allocation
    with allocation_guard():host_capacity(1,minimum_ram_gib=4);os.sched_setaffinity(0,cpu_allocation(1))
    if a.out.exists():raise FileExistsError('Preserve earlier sampling inventories')
    result=inventory(a.archive.resolve(),a.census.resolve(),a.energy_roster.resolve(),[p.resolve() for p in a.reference_root]);atomic_write_json(a.out,result)
    print(json.dumps({k:result[k] for k in ('status','declared_inputs','selected_inputs','missing','source_kind_counts')}))
