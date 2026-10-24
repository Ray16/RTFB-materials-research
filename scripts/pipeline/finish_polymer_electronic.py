#!/usr/bin/env python
"""Compute canonical CPU energies after a separately checked polymer geometry."""
import argparse
from contextlib import redirect_stdout, redirect_stderr
import json
import os
from pathlib import Path
import shutil
import signal
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def prepare(geometry, pilot, out):
    with allocation_guard():
        host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    original, control = api.read_plan(geometry), api.read_plan(pilot)
    for plan in (original, control):
        api.verify_sources(plan['producer_sha256'])
        api.verify_sources(plan['input_sha256'])
    completed = json.loads((geometry/'geometry/status.json').read_text())
    qualified = json.loads((pilot/'electronic/status.json').read_text())
    if completed['status'] != 'stationary_geometry_checked_active_parent_and_native_minimum_pending':
        raise ValueError('The geometry must complete graph, stationarity and derivative checks')
    if qualified['status'] != 'split_electronic_resource_and_reference_point_passed':
        raise ValueError('The full-size CPU/reference resource point must pass first')
    if control.get('CPU_PCM_memory_MiB') is None:
        raise ValueError('The bounded CPU PCM execution setting must be explicit')
    budget_path = pilot/'PCM_budget/status.json'
    budget = json.loads(budget_path.read_text())
    if (budget['status'] != 'CPU_PCM_budget_finished_requires_numerical_review'
            or budget['PCM_memory_MiB'] != control['CPU_PCM_memory_MiB']
            or min(budget['tiles'].values()) < 1 or budget['maximum_tile_points'] > 400
            or not budget['started'] <= qualified['started'] <= qualified['finished'] <= budget['finished']):
        raise ValueError('The completed large point must exercise the bounded PCM execution')
    api.verify_sources(budget['source_sha256'])
    for record in (completed, qualified):
        api.verify_sources(record['output_sha256'])
    owner = api.load_module('split_electronic_owner_check', ROOT/'scripts/pipeline/run_polymer_spacing.py')
    if owner.owner_alive(completed['owner']) or owner.owner_alive(qualified):
        raise ValueError('Both predecessor workers must finish before handoff')
    if original['tasks'] != [control['task']] or original['quantum_settings'] != control['quantum_settings']:
        raise ValueError('The CPU resource qualification belongs to another graph/protocol')
    task = original['tasks'][0]
    if (task['charge'], task['mult']) != (0, 1):
        raise ValueError('The full-size CPU resource profile currently covers neutral singlets')
    out.mkdir(parents=True, exist_ok=False)
    shutil.copytree(geometry/'controls', out/'controls')
    for stage in ('energy', 'thermal'):
        for name in ('tasks','done','claims','failures','workers','pilot','launches'):
            (out/stage/name).mkdir(parents=True)
    sources = dict(original['producer_sha256'], **control['producer_sha256'])
    for file in (Path(__file__).resolve(), ROOT/'src/redox/workflows/split_polymer.py'):
        sources[str(file)] = api.digest(file)
    inputs = dict(original['input_sha256'], **control['input_sha256'])
    for file in (geometry/'geometry/status.json', geometry/'plan.json', pilot/'plan.json', pilot/'electronic/status.json', budget_path):
        inputs[str(file)] = api.digest(file)
    inputs.update(completed['output_sha256'])
    inputs.update(qualified['output_sha256'])
    plan = dict(original, raw_root=str(ROOT/'calcs/polymer_models'/out.parent.name/out.name),
        producer_sha256=sources, input_sha256=inputs, source_geometry_report=completed,
        CPU_resource_policy=control['electronic'], CPU_qualification=str(pilot),
        CPU_PCM_memory_setting=control['CPU_PCM_memory_MiB'],
        CPU_PCM_memory_setting_unit='PySCF max_memory MB; inherited key suffix retained in pilot schema',
        geometry_only_execution=False, energy_only_execution=True, native_pipeline_enabled=False,
        externally_managed_stages=True, created=time.time(),
        purpose='Checked geometry plus canonical CPU active energy; native minimum/G is a separate stage')
    atomic_write_json(out/'plan.json', plan)
    print(json.dumps(dict(status='CPU_electronic_stage_prepared', tasks=1)))


def run(campaign, threads, scratch):
    import pyscf.lib
    from redox.core import protocol
    from redox.qm import sp
    from redox.workflows.checkpoints import verify_state_checkpoint
    from redox.workflows.split_polymer import composite_result

    plan = api.read_plan(campaign)
    api.verify_sources(plan['producer_sha256'])
    api.verify_sources(plan['input_sha256'])
    settings = plan['CPU_resource_policy']
    caps = ('OMP_NUM_THREADS','OMP_THREAD_LIMIT','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','XTB_THREADS')
    if threads < 1 or any(os.environ.get(k) != str(threads) for k in caps):
        raise ValueError('All native thread caps must match the CPU allocation')
    with allocation_guard():
        resources = host_capacity(threads, minimum_ram_gib=settings['minimum_host_RAM_GiB'])
        if settings['maximum_RSS_GiB']*2**30 > .75*resources['available_RAM_bytes']:
            raise ValueError('The host cannot cover the declared CPU memory bound')
        os.sched_setaffinity(0, cpu_allocation(threads))
    if not scratch.is_dir() or os.environ.get('PYSCF_TMPDIR') != str(scratch):
        raise ValueError('Task-owned local scratch must be configured before PySCF import')
    if shutil.disk_usage(scratch).free < (settings['maximum_scratch_GiB']+64)*2**30:
        raise ValueError('Insufficient task-local scratch capacity')
    task = plan['tasks'][0]
    key = task['key']
    directory = campaign/'energy/tasks'/key
    state = Path(plan['raw_root'])/'energy'/key/'rest'
    if directory.exists() or state.exists():
        raise ValueError('Preserve the existing attempt; use a reviewed fresh campaign')
    claim = campaign/'energy/claims'/key
    claim.mkdir()
    owner = dict(host=os.uname().nodename, pid=os.getpid(), process_start=Path('/proc/self/stat').read_text().rsplit(')',1)[1].split()[19])
    atomic_write_json(claim/'owner.json', owner)
    source = Path(plan['source_geometry_report']['state_dir'])
    report = dict(status='computing', task=task, state_dir=str(state), owner=owner,
        started=time.time(), CPU_affinity=sorted(os.sched_getaffinity(0)), native_G_computed=False,
        geometry_parent=str(source), CPU_resource_qualification=plan['CPU_qualification'])
    def stop(signum, frame):
        raise KeyboardInterrupt('Owned CPU stop; checkpoints preserved')
    previous_handler = signal.signal(signal.SIGTERM, stop)
    previous_versions = sp.software_versions
    try:
        directory.mkdir(exist_ok=False)
        state.mkdir(parents=True, exist_ok=False)
        shutil.copyfile(source/'opt.xyz', state/'opt.xyz')
        shutil.copyfile(source/'result.json', state/'geometry_parent.json')
        shutil.copyfile(source/'result.json', state/'result.json')
        atomic_write_json(directory/'status.json', report)
        adapter = api.load_module('canonical_polymer_CPU_execution', ROOT/'scripts/pipeline/cpu_quantum_execution.py')
        resource = api.load_module('canonical_polymer_CPU_resource_bound', ROOT/'scripts/validation/qualify_split_polymer.py')
        pcm = api.load_module('canonical_polymer_CPU_PCM_budget', ROOT/'scripts/pipeline/cpu_pcm_memory.py')
        pyscf.lib.num_threads(threads)
        with resource.CPU_resources(scratch, settings, directory), adapter.cpu_execution(state/'CPU_execution', settings['library_memory_MiB']), pcm.bounded_pcm_memory(state/'PCM_budget', plan['CPU_PCM_memory_setting']), (directory/'calculation.log').open('x') as log, redirect_stdout(log), redirect_stderr(log):
            # sp imported this function by value. Bind its metadata-only CPU
            # replacement too, so recording versions never initializes CUDA.
            sp.software_versions = protocol.software_versions
            _, active = sp.ensure(key, 'rest', state, task['charge'], task['mult'], backend='cpu')
            checked = protocol.load_record(state, task['charge'], task['mult'])
            if checked is None:
                raise ValueError('The canonical CPU active record failed normal admission')
            result = composite_result(json.loads((state/'geometry_parent.json').read_text()), checked,
                protocol.geom_sha1(state/'opt.xyz'), protocol.record_path(state, task['charge'], task['mult']))
            atomic_write_json(state/'result.json', result)
        resources = json.loads((directory/'resources.json').read_text())
        if (not resources['samples'] or (directory/'resource_stop.json').exists()
                or (directory/'resource_monitor_error.json').exists()):
            raise ValueError('CPU resource observation was incomplete or interrupted')
        api.verify_sources(plan['producer_sha256'])
        api.verify_sources(plan['input_sha256'])
        report.update(status='state_records_verified', E_smd_eV=checked['e_smd_eV'], physical_QC='passed',
            physical_QC_source='Inherited exact-byte geometry and force qualification',
            topology_QC='ok', gradient_max_Ha_Bohr=plan['source_geometry_report']['gradient_max_Ha_Bohr'],
            conservative_total_force_check=plan['source_geometry_report']['conservative_total_force_check'],
            output_sha256={str(p):api.digest(p) for p in state.rglob('*') if p.is_file() and p.suffix in ('.json','.xyz')}, finished=time.time())
        report['output_sha256'][str(directory/'resources.json')] = api.digest(directory/'resources.json')
        atomic_write_json(directory/'status.json', report)
        receipt = campaign/'energy/done'/(key+'.json')
        atomic_write_json(receipt, dict(key=key, status_path=str(directory/'status.json'), status_sha256=api.digest(directory/'status.json'), owner=owner))
        verify_state_checkpoint(receipt)
    except BaseException as error:
        report.update(status='failed_preserved', error=repr(error), finished=time.time())
        atomic_write_json(directory/'status.json', report)
        raise
    finally:
        sp.software_versions = previous_versions
        signal.signal(signal.SIGTERM, previous_handler)
        shutil.rmtree(claim)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare','run'))
    for name in ('geometry','pilot','campaign','scratch'):
        parser.add_argument('--'+name, type=Path)
    parser.add_argument('--threads', type=int, default=16)
    args = parser.parse_args()
    if args.campaign is None:
        parser.error('--campaign required')
    if args.action == 'prepare':
        if args.geometry is None or args.pilot is None:
            parser.error('prepare requires --geometry and --pilot')
        prepare(args.geometry.resolve(), args.pilot.resolve(), args.campaign.resolve())
    else:
        if args.scratch is None:
            parser.error('run requires --scratch')
        run(args.campaign.resolve(), args.threads, args.scratch.resolve())


if __name__ == '__main__':
    main()
