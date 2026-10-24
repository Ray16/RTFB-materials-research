#!/usr/bin/env python
"""Bounded single-conformer SCF/gradient pilot from a final CREST queue.

The screened-array forecast remains visible. This separate low-RAM experiment
must finish within enforced RSS/free-RAM limits; it does not relax physical QC.
"""
import argparse
import json
import math
import os
from pathlib import Path
import shutil
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import verify_gate_owner, host_capacity, gpu_suitability

ADAPTER = ROOT/'scripts/pipeline/smd_switch_curvature_workspace.py'
CHECKER = ROOT/'scripts/pipeline/switch_curvature_polymer_campaign.py'


def prepare(source, inventory, out):
    old = api.read_plan(source)
    api.verify_sources(old['producer_sha256']); api.verify_sources(old['input_sha256'])
    proof = json.loads((inventory/'status.json').read_text())
    if proof['status'] != 'screened_DF_inventory_complete_full_SCF_not_run' or proof['source'] != str(source):
        raise ValueError('A completed inventory for the same source is required')
    api.verify_sources(proof['source_sha256'])
    tasks = [t for t in old['tasks'] if t['key'] == proof['task_key']]
    if len(tasks) != 1 or (tasks[0]['charge'],tasks[0]['mult']) != (0,1):
        raise ValueError('One exact neutral-singlet source task is required')
    controls = ROOT/'results/validation/pcm_switch_curvature_guarded_2026-10-23'
    witness = api.load_module('CREST_point_controls',CHECKER).checked_controls(controls,28)
    out.mkdir(parents=True,exist_ok=False)
    shutil.copytree(source/'controls',out/'controls')
    key = tasks[0]['key']; target = out/'live_inputs'/key
    shutil.copytree(source/'live_inputs'/key,target)
    ready = json.loads((target/'task.json').read_text())
    ready['input_xyz'] = str(target/'input.xyz')
    ready['input_sha256'][ready['input_xyz']] = api.digest(target/'input.xyz')
    atomic_write_json(target/'task.json',ready)
    for stage in ('energy','thermal'):
        for folder in ('tasks','done','claims','failures','workers','pilot','launches'):
            (out/stage/folder).mkdir(parents=True)
    inputs = dict(old['input_sha256'],**ready['input_sha256'])
    for p in (source/'plan.json',inventory/'status.json',target/'task.json'):
        inputs[str(p)] = api.digest(p)
    sources = dict(old['producer_sha256'],**witness)
    for p in (Path(__file__).resolve(),ADAPTER,CHECKER):sources[str(p)] = api.digest(p)
    new = dict(old,tasks=tasks,pilot_key=key,source_campaign=str(source),input_sha256=inputs,
        producer_sha256=sources,raw_root=str(ROOT/'calcs/validation'/out.parent.name/out.name),
        minimum_large_GPU_MiB=32768,maximum_memory_pilot_MiB=31744,
        execution_workspace_GiB=28,maximum_RSS_GiB=80,minimum_free_RAM_GiB=32,
        inventory_forecast_host_bytes=proof['required_host_bytes'],
        resource_experiment='Enforced 80 GiB RSS plus 32 GiB free-RAM floor; distinct from the conservative 1.2*DF+24 GiB forecast',
        purpose='Fixed-geometry gas/SMD SCF and gradient resource experiment; no minimum, native G or ensemble',
        externally_managed_stages=True,created=time.time())
    atomic_write_json(out/'plan.json',new)


def run(campaign):
    verify_gate_owner()
    caps=('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
    if len(os.sched_getaffinity(0)) != 2 or any(os.environ.get(k) != '2' for k in caps):
        raise ValueError('Two disjoint CPUs and native caps of two are required')
    host_capacity(2,minimum_ram_gib=112)
    if gpu_suitability(int(os.environ['CUDA_VISIBLE_DEVICES']))['memory_MiB'] < 32768:
        raise ValueError('The energy pilot requires a 32 GiB device')
    plan = api.read_plan(campaign)
    api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])
    # The external same-PID host guard is part of the execution contract.
    guard = json.loads((campaign/'host_guard/status.json').read_text())
    if guard['pid'] != os.getpid() or guard['maximum_RSS_GiB'] != 80 or guard['minimum_free_GiB'] != 32:
        raise ValueError('The bounded host guard must own this process')
    runtime = campaign/'resource_point'; runtime.mkdir(exist_ok=False)
    record = dict(status='CREST_resource_point_running',started=time.time(),host=os.uname().nodename,
        pid=os.getpid(),process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19],
        native_G_qualified=False,optimized_minimum_qualified=False)
    atomic_write_json(runtime/'status.json',record)
    try:
        driver = api.load_module('CREST_normal_energy_driver',ROOT/'scripts/pipeline/refine_polymer_conformers.py')
        task = driver.backend(campaign,'energy').task_input(campaign,plan['tasks'][0])
        adapter = api.load_module('CREST_bounded_point_execution',ADAPTER)
        with adapter.compact_smd_execution(Path(plan['raw_root'])/'execution',28):
            driver.calculate_energy(campaign,task,pilot=True)
        status = json.loads((campaign/'energy/pilot/status.json').read_text())
        data = json.loads((campaign/'energy/pilot/energy_and_gradient.json').read_text())
        if status['status'] != 'memory_pilot_passed' or not 0 < status['observed_peak_GPU_MiB'] <= 31744:
            raise ValueError('Original energy/gradient pilot and GPU bounds must pass')
        if any(not math.isfinite(data['result'][k]) for k in ('e_smd_eV','e_gas_eV')):
            raise ValueError('Both energies must be finite')
        api.verify_sources(status['output_sha256'])
        api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])
        record.update(status='CREST_point_finished_independent_resource_review_pending',
            point_sha256={str(p):api.digest(p) for p in (campaign/'energy/pilot/status.json',campaign/'energy/pilot/energy_and_gradient.json')})
    except BaseException as error:
        record.update(status='CREST_point_failed_preserved',error=repr(error));raise
    finally:
        record['finished']=time.time();atomic_write_json(runtime/'status.json',record)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','run'))
    for name in ('source','inventory','campaign'):p.add_argument('--'+name,type=Path)
    a=p.parse_args()
    if a.campaign is None:p.error('--campaign required')
    if a.action == 'prepare':
        if a.source is None or a.inventory is None:p.error('prepare requires source and inventory')
        prepare(a.source.resolve(),a.inventory.resolve(),a.campaign.resolve())
    else:run(a.campaign.resolve())


if __name__=='__main__':main()
