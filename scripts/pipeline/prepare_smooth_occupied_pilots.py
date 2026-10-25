#!/usr/bin/env python
"""Prepare one bounded geometry/resource pilot per occupied control graph.

These are prospective resource measurements, not propagation of the empty-chain
pilot's memory certificate. No source geometry, energy or free energy is admitted
by preparing a plan. Expand each class only after independent pilot review.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def choose_task(plan):
    tasks = plan['tasks']
    if not tasks or len({t['key'] for t in tasks}) != len(tasks):
        raise ValueError('Distinct declared conformer inputs are required')
    task = tasks[0]
    if (task['candidate_id'] != 'exp_phthalimide_control'
            or not task['occupied_sites'] or (task['charge'], task['mult']) != (0, 1)
            or task['atoms'] not in (105, 119)):
        raise ValueError('Only the existing neutral occupied phthalimide controls are in scope')
    return task


def prepare(sources, controls, out):
    model_path = ROOT / 'scripts/pipeline/smooth_cds_control.py'
    worker = ROOT / 'scripts/pipeline/run_smooth_geometry.py'
    coordinator = ROOT / 'scripts/pipeline/coordinate_smooth_geometry.py'
    runtime = ROOT / 'environments/cutensor22/runtime.json'
    model = api.load_module('occupied_pilot_declared_model', model_path)
    checker = ROOT / 'scripts/validation/review_smooth_model_panel.py'
    evidence = api.load_module('occupied_pilot_quantum_controls', checker).checked_controls(controls)
    identity = hashlib.sha256(json.dumps(model.MODEL, sort_keys=True).encode()).hexdigest()
    entries = []
    prepared = []
    # Validate the whole request before creating output directories.
    for source in sources:
        plan = api.read_plan(source)
        api.verify_sources(plan['producer_sha256']); api.verify_sources(plan['input_sha256'])
        task = choose_task(plan)
        prepared.append((source, plan, task))
    if len({t['model_key'] for _, _, t in prepared}) != len(prepared):
        raise ValueError('Each graph/configuration gets exactly one resource pilot')
    out.mkdir(parents=True, exist_ok=False)
    for source, original, task in prepared:
        case = out / task['key']; case.mkdir()
        ready = case / 'live_inputs' / task['key']
        shutil.copytree(source / 'live_inputs' / task['key'], ready)
        xyz = ready / 'input.xyz'
        record = json.loads((ready / 'task.json').read_text())
        record['input_xyz'] = str(xyz); record['input_sha256'][str(xyz)] = api.digest(xyz)
        atomic_write_json(ready / 'task.json', record)
        from ase.io import read
        symbols = read(xyz).get_chemical_symbols()
        if len(symbols) != task['atoms']:
            raise ValueError('The declared atom count differs from the actual pilot input')
        audit = api.load_module('occupied_pilot_graph_check', ROOT / 'scripts/validation/audit_polymer_geometries.py')
        graph = audit.audit(task, xyz)
        if graph['issues']:
            raise ValueError('Pilot input graph/stereochemistry failed')
        atomic_write_json(ready / 'initial_graph_review.json', graph)
        inputs = {**original['input_sha256'], **evidence}
        for p in (source / 'plan.json', xyz, ready / 'task.json', ready / 'initial_graph_review.json', runtime):
            inputs[str(p)] = api.digest(p)
        producers = {**original['producer_sha256'], **evidence}
        for p in (Path(__file__).resolve(), model_path, worker, coordinator, checker):
            producers[str(p)] = api.digest(p)
        for stage in ('energy', 'thermal'):
            for name in ('tasks', 'done', 'claims', 'failures', 'workers', 'pilot', 'launches'):
                (case / stage / name).mkdir(parents=True)
        hetero = next(i for i, symbol in enumerate(symbols) if symbol not in ('C', 'H'))
        plan = dict(original, tasks=[task], input_sha256=inputs, producer_sha256=producers,
            raw_root=str(ROOT / 'calcs/polymer_models' / out.name / task['key']),
            execution_driver=str(worker), runtime_profile=str(runtime), maximum_host_RSS_GiB=32,
            quantum_settings=dict(original['quantum_settings'], geometry_SCF_conv_tol=1e-12),
            geometry_resource_policy=dict(minimum_GPU_MiB=32768, maximum_peak_GPU_MiB=31744,
                                          minimum_host_RAM_GiB=64),
            experimental_model=model.MODEL, experimental_method_id=identity,
            full_quantum_controls=str(controls), extra_force_axes=[[0, 1], [hetero, 2]],
            primary_candidate_admission=False, native_pipeline_enabled=False,
            geometry_only_execution=True, energy_only_execution=True, externally_managed_stages=True,
            prospective_resource_pilot=True, measured_alternate_model_reference=None,
            resource_scope='One guarded 32-GB-class measurement; no resource certificate inherited from another graph',
            purpose='Occupied phthalimide method-control geometry/resource pilot; no synthesis eligibility',
            created=time.time())
        atomic_write_json(case / 'plan.json', plan)
        entries.append(dict(key=task['key'], campaign=str(case), plan_sha256=api.digest(case / 'plan.json')))
    sources = {str(p): api.digest(p) for p in (Path(__file__).resolve(), model_path, worker, coordinator,
        ROOT / 'scripts/pipeline/coordinate_geometry_batch.py', runtime)}
    atomic_write_json(out / 'batch.json', dict(status='prospective_occupied_geometry_pilots_prepared',
        stage='geometry', entries=entries, worker=str(worker), maximum_RSS_GiB=32,
        minimum_free_RAM_GiB=32, minimum_GPU_MiB=32768, experimental_model=model.MODEL,
        primary_candidate_admission=False, source_sha256=sources, quantum_results_transferred=False))
    print(json.dumps(dict(prepared_resource_pilots=len(entries), resource_qualification_inferred=False)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, action='append', required=True)
    parser.add_argument('--controls', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    import os
    with allocation_guard():
        host_capacity(1); os.sched_setaffinity(0, cpu_allocation(1))
    prepare([p.resolve() for p in args.source], args.controls.resolve(), args.out.resolve())


if __name__ == '__main__':
    main()
