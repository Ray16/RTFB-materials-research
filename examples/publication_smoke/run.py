#!/usr/bin/env python
"""Small CPU-only reproduction of exact design -> descriptors -> route gate -> figure."""
import argparse
from importlib.metadata import version
from itertools import islice
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.screening import route_stage, streaming_search
from redox.workflows import chunked
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def outputs(campaign):
    plan = chunked.read_plan(campaign)
    values = {}
    for part in plan['chunks']:
        report = chunked.verify_done(campaign / 'done' / (part['key'] + '.json'), plan=plan)
        if report['failed']:
            raise ValueError('A failed task cannot enter the example result')
        for line in Path(report['results']).read_text().splitlines():
            row = json.loads(line)
            values[row['task_id']] = row['value']
    if len(values) != plan['unique_tasks']:
        raise ValueError('Example output coverage differs from its exact input')
    return plan, values


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, required=True, help='New output directory inside this checkout')
    parser.add_argument('--plot-python', type=Path, default=Path(sys.executable), help='Interpreter with Matplotlib')
    args = parser.parse_args()
    for key in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'XTB_THREADS'):
        if os.environ.get(key) != '1':
            raise ValueError('Export all six native thread caps as 1 before running this example')
    destination = args.out.resolve()
    if not destination.is_relative_to(ROOT):
        raise ValueError('Use an output inside the checkout so figure provenance stays relative')
    with allocation_guard():
        capacity = host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    destination.mkdir(parents=True, exist_ok=False)
    sources = [Path(__file__), ROOT / 'src/redox/workflows/chunked.py',
               ROOT / 'src/redox/screening/streaming_search.py', ROOT / 'src/redox/screening/route_stage.py',
               ROOT / 'scripts/plotting/polymer_models/plot_conditional_grafting.py',
               ROOT / 'manuscript/evidence/2026-10-22/polymer_validation.json']
    identities = {str(p.relative_to(ROOT)):chunked.file_hash(p) for p in sources}
    # This is a declared small subset, not a new candidate search or synthesis shortlist.
    campaign = destination / 'descriptors'
    chunked.prepare(islice(streaming_search.family_records(), 32), campaign,
                    streaming_search.stage_spec('descriptors'), chunk_size=4)
    plan = chunked.read_plan(campaign)
    chunked.run_worker(campaign, streaming_search.handler_for(plan), checkpoint_every=2)
    plan, actual = outputs(campaign)
    streaming_search.verify(campaign)
    for part in plan['chunks']:
        for line in (campaign / part['path']).read_text().splitlines():
            row = json.loads(line)
            key = row.pop('task_id')
            if actual[key] != streaming_search.descriptor(row):
                raise ValueError('Checkpointed descriptors differ from serial replay')

    empty = destination / 'empty_route_evidence.json'
    atomic_write_json(empty, dict(candidates={}))
    routes = destination / 'routes'
    stage = route_stage.stage(campaign, empty, empty)
    chunked.prepare(route_stage.records(campaign, empty, empty), routes, stage, chunk_size=4)
    chunked.run_worker(routes, route_stage.screen, checkpoint_every=2)
    _, decisions = outputs(routes)
    route_report = route_stage.verify(routes)
    if any(value['route_status'] != 'hold' or value['proceed_to_full_redox'] for value in decisions.values()):
        raise ValueError('Missing exact route evidence improperly opened the quantum gate')

    evidence = ROOT / 'manuscript/evidence/2026-10-22/polymer_validation.json'
    manifest = json.loads((evidence.parent / 'manifest.json').read_text())
    if chunked.file_hash(evidence) != manifest['snapshot_sha256'][str(evidence.relative_to(ROOT))]:
        raise ValueError('Frozen scientific example evidence changed')
    figures = destination / 'figures'
    command = [str(args.plot_python), str(ROOT / 'scripts/plotting/polymer_models/plot_conditional_grafting.py'),
               '--review', str(evidence), '--out', str(figures)]
    with (destination / 'figure.log').open('w') as log:
        subprocess.run(command, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, check=True)
    reference = ROOT / 'results/figures/polymer_models/conditional_grafting.csv'
    if (figures / reference.name).read_bytes() != reference.read_bytes():
        raise ValueError('Regenerated conditional figure data differs from the published table')
    if any(chunked.file_hash(ROOT / p) != h for p, h in identities.items()):
        raise ValueError('Example source changed during execution; preserve outputs without a verified claim')
    report = dict(status='small_CPU_workflow_and_frozen_figure_reproduced', recipes=plan['input_records'],
                  exact_designs=plan['unique_tasks'], serial_descriptor_discrepancies=0,
                  route_decisions=route_report['counts'], quantum_jobs_submitted=0,
                  figure_data_matches=True, quantum_energies_recomputed=False,
                  python=platform.python_version(), executable=sys.executable,
                  dependencies={name:version(name) for name in ['rdkit','numpy','scipy','psutil']},
                  capacity=capacity, CPU_affinity=sorted(os.sched_getaffinity(0)),
                  root=str(ROOT), source_sha256=identities,
                  limits='Small CPU execution and published-figure reproduction; not a fresh quantum calculation or complete environment qualification')
    atomic_write_json(destination / 'summary.json', report)
    print(json.dumps(report))


if __name__ == '__main__':
    main()
