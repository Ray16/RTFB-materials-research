#!/usr/bin/env python
"""Inventory exact screened DF storage before allocating a polymer SCF task.

This runs solvent/auxiliary construction under a reservation and intercepts host
DF allocation. It computes no SCF energy, optimized geometry, Hessian or G.
"""
import argparse
import json
import os
from pathlib import Path
import resource
import time

import psutil

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import verify_gate_owner, host_capacity, gpu_suitability


def required_host_bytes(shapes):
    if len(shapes) != 2:
        raise ValueError('Inventory both full and range-separated DF arrays')
    for item in shapes:
        shape = item['shape']
        if (len(shape) != 2 or any(type(n) is not int or n <= 0 for n in shape)
                or item['bytes'] != shape[0]*shape[1]*8):
            raise ValueError('Inconsistent FP64 DF shape inventory')
    return int(1.2*sum(row['bytes'] for row in shapes)+24*2**30)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--task-key', required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    verify_gate_owner()
    caps = ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
    if len(os.sched_getaffinity(0)) != 2 or any(os.environ.get(k) != '2' for k in caps):
        raise ValueError('Two disjoint CPUs and matching native caps are required')
    host_capacity(2, minimum_ram_gib=32)
    device = gpu_suitability(int(os.environ['CUDA_VISIBLE_DEVICES']))
    if device['memory_MiB'] < 32768:
        raise ValueError('The unqualified inventory requires a 32 GiB GPU')
    source = args.source.resolve()
    plan = api.read_plan(source)
    api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])
    tasks = [row for row in plan['tasks'] if row['key'] == args.task_key]
    if len(tasks) != 1 or (tasks[0]['charge'],tasks[0]['mult']) != (0,1):
        raise ValueError('A single exact neutral-singlet task is required')
    controls = ROOT/'results/validation/pcm_switch_curvature_guarded_2026-10-23'
    checker_path = ROOT/'scripts/pipeline/switch_curvature_polymer_campaign.py'
    checker = api.load_module('inventory_full_controls', checker_path)
    evidence = checker.checked_controls(controls, 28)
    args.out.mkdir(parents=True, exist_ok=False)
    (args.out/'producer_source.py').write_bytes(Path(__file__).read_bytes())
    report = dict(status='inventory_running_no_SCF', host=os.uname().nodename, pid=os.getpid(),
        process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
        task_key=args.task_key, source=str(source), CPU_affinity=sorted(os.sched_getaffinity(0)),
        device=device, started=time.time(), SCF_energy_computed=False, native_G_qualified=False,
        control_evidence=evidence, source_sha256={str(Path(__file__).resolve()):api.digest(__file__),
            str(source/'plan.json'):api.digest(source/'plan.json')})
    atomic_write_json(args.out/'status.json', report)
    try:
        driver = api.load_module('inventory_original_driver', ROOT/'scripts/pipeline/run_polymer_spacing.py')
        task = driver.task_input(source, tasks[0])
        implementation = ROOT/'scripts/validation/qualify_polymer_workspace.py'
        adapter_path = ROOT/'scripts/pipeline/smd_switch_curvature_workspace.py'
        inventory = api.load_module('screened_DF_shape_inventory', implementation)
        adapter = api.load_module('inventory_bounded_workspace', adapter_path)
        monitor = api.load_module('inventory_GPU_memory', ROOT/'scripts/pipeline/run_pi_stacking.py')
        from gpu4pyscf.df import df as installed_df
        report['source_sha256'][str(Path(installed_df.__file__))] = api.digest(installed_df.__file__)
        with monitor.MemoryMonitor() as memory, adapter.compact_smd_execution(args.out/'execution',28):
            shapes = inventory.tensor_shapes(task)
        report.update(DF_tensor_shapes=shapes,required_host_bytes=required_host_bytes(shapes),
                      observed_peak_GPU_MiB=memory.peak)
        atomic_write_json(args.out/'status.json', report)
        if not 0 < memory.peak <= 31744:
            raise ValueError('GPU telemetry is missing or exceeds the inventory bound: '+str(memory.peak))
        required = required_host_bytes(shapes)
        available = psutil.virtual_memory().available
        report.update(status='screened_DF_inventory_complete_full_SCF_not_run',DF_tensor_shapes=shapes,
            required_host_bytes=required,available_host_RAM_bytes=available,
            enough_current_RAM_before_other_worker_promises=required <= .75*available,
            observed_peak_GPU_MiB=memory.peak,
            scope='Storage forecast for this exact input; no energy or execution-profile qualification')
        report['source_sha256'].update({str(p):api.digest(p) for p in (implementation,adapter_path,checker_path)})
        api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])
    except BaseException as error:
        report.update(status='inventory_failed_preserved_no_SCF',error=repr(error))
        raise
    finally:
        report.update(finished=time.time(),peak_RSS_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        atomic_write_json(args.out/'status.json', report)


if __name__ == '__main__':
    main()
