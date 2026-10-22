#!/usr/bin/env python
"""Measure actual descriptor-stage concurrency, including launch and idle workers."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import random
import statistics
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.screening import streaming_search as search
from redox.workflows import chunked, cpu_chunks
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity


def workers(campaign):
    values = []
    for directory in sorted((campaign / 'workers').iterdir()):
        if not (directory / 'status.json').exists():
            continue
        values.append(dict(owner=json.loads((directory / 'owner.json').read_text()),
                           status=json.loads((directory / 'status.json').read_text())))
    return values


def audit(campaign, expected, count):
    plan = chunked.read_plan(campaign)
    total = 0
    for part in plan['chunks']:
        receipt = chunked.verify_done(campaign / 'done' / (part['key'] + '.json'), plan=plan)
        if receipt['failed']:
            raise ValueError('Descriptor failures are not benchmark completions')
        for row in map(json.loads, Path(receipt['results']).read_text().splitlines()):
            actual = row['value']
            if any(actual.get(key) != value for key, value in expected[row['task_id']].items()):
                raise ValueError('Parallel output differs from the independent serial definition')
            if actual['solution_properties_status'] != 'not computed' or actual['synthesis_eligibility'] != 'not established':
                raise ValueError('Cheap descriptors acquired unsupported physical or synthesis QC')
            total += 1
    if total != len(expected) or total != plan['unique_tasks']:
        raise ValueError('Benchmark coverage differs')
    records = workers(campaign)
    if len(records) != count or any(r['status']['status'] != 'worker_finished' for r in records):
        raise ValueError('Requested workers have not all finished')
    for a in records:
        if len(a['owner']['cpus']) != 1 or any(v != '1' for v in a['owner']['thread_caps'].values()):
            raise ValueError('Actual worker CPU affinity or native thread caps differ')
        for b in records:
            if a is b or a['owner']['host'] != b['owner']['host']:
                continue
            overlap = max(a['status']['started'], b['status']['started']) < min(a['status']['finished'], b['status']['finished'])
            if overlap and set(a['owner']['cpus']) & set(b['owner']['cpus']):
                raise ValueError('Concurrent workers used overlapping CPU sets')
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config', type=Path, default=ROOT / 'config/cpu_concurrency_benchmark.json')
    parser.add_argument('--raw', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    with allocation_guard():
        capacity = host_capacity(1)
        os.sched_setaffinity(0, cpu_allocation(1))
    if any(os.environ.get(key) != '1' for key in ('OMP_NUM_THREADS', 'OMP_THREAD_LIMIT', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'NUMEXPR_NUM_THREADS', 'XTB_THREADS')):
        raise ValueError('All six native caps must be one')
    config = json.loads(args.config.read_text())
    if type(config['repeats']) is not int or config['repeats'] < 2 or config['chunk_size'] < 1:
        raise ValueError('Repeated, positive-size experiment required')
    if len({c['name'] for c in config['conditions']}) != len(config['conditions']):
        raise ValueError('Repeated condition name')
    source = [Path(__file__), args.config.resolve(), ROOT / 'src/redox/workflows/cpu_chunks.py',
              ROOT / 'src/redox/workflows/chunked.py', ROOT / 'src/redox/screening/streaming_search.py',
              ROOT / 'scripts/validation/audit_chunk_benchmark.py']
    identities = {str(p):chunked.file_hash(p) for p in source}
    raw, out = args.raw.resolve(), args.out.resolve()
    raw.mkdir(parents=True, exist_ok=False); out.mkdir(parents=True, exist_ok=False)
    atomic_write_json(out / 'plan.json', dict(config=config, source_sha256=identities,
                       role='Real descriptor execution benchmark; not quantum throughput', capacity=capacity))
    serial_path = ROOT / 'scripts/validation/audit_chunk_benchmark.py'
    spec = importlib.util.spec_from_file_location('independent_serial_descriptor', serial_path)
    serial = importlib.util.module_from_spec(spec); spec.loader.exec_module(serial)
    queue = [(condition, repeat) for condition in config['conditions'] for repeat in range(config['repeats'])]
    random.Random(config['seed']).shuffle(queue)
    results, expected = [], None
    atomic_write_json(raw / 'randomized_order.json', [dict(condition=c['name'], repeat=r) for c,r in queue])
    for condition, repeat in queue:
        if any(chunked.file_hash(p) != h for p,h in identities.items()):
            raise ValueError('Benchmark producers changed')
        name = condition['name'] + '_r' + str(repeat)
        campaign = raw / name
        plan = chunked.prepare(search.family_records(), campaign, search.stage_spec('descriptors'), chunk_size=config['chunk_size'])
        if expected is None:
            expected = {}
            for part in plan['chunks']:
                for line in (campaign / part['path']).read_text().splitlines():
                    row = json.loads(line); key = row.pop('task_id')
                    if chunked.content_hash(row) != key:
                        raise ValueError('Exact task hash mismatch')
                    expected[key] = serial.serial_descriptor(row)
            atomic_write_json(raw / 'independent_serial.json', expected)
        started = time.time()
        cpu_chunks.launch(campaign, condition['hosts'], raw / (name + '_launch'))
        count = sum(condition['hosts'].values())
        deadline = time.monotonic() + 180
        while True:
            records = workers(campaign)
            if len(records) == count and all(r['status'].get('status') == 'worker_finished' for r in records):
                break
            if any(r['status'].get('status') == 'stopped_preserving_outputs' for r in records):
                raise RuntimeError('A real worker failed; preserve campaign for review')
            if time.monotonic() >= deadline:
                raise RuntimeError('Observation deadline; retain actual worker state without restarting it')
            time.sleep(.1)
        observed_complete = time.time()
        records = audit(campaign, expected, count)
        search.verify(campaign)
        begins = [r['status']['started'] for r in records]
        ends = [r['status']['finished'] for r in records]
        result = dict(condition=condition['name'], repeat=repeat, hosts=condition['hosts'],
            requested_workers=count, observed_workers=len(records),
            workers_executing_chunks=sum(bool(r['status']['completed_chunks']) for r in records),
            tasks=plan['unique_tasks'], input_recipes=plan['input_records'], chunks=len(plan['chunks']),
            preparation_seconds=plan['prepared_seconds'],
            launch_to_observed_completion_seconds=observed_complete-started,
            worker_timestamp_span_seconds=max(ends)-min(begins),
            largest_worker_peak_RSS_KiB=max(r['status']['peak_RSS_KiB'] for r in records),
            serial_discrepancies=0, workers=records,
            raw_campaign=str(campaign), manifest_sha256=plan['manifest_sha256'])
        results.append(result)
        atomic_write_json(out / 'progress.json', dict(status='running', completed=len(results), expected=len(queue), latest=result))
        print(json.dumps(dict(completed=len(results), condition=result['condition'], repeat=repeat,
                             seconds=result['launch_to_observed_completion_seconds'])), flush=True)
    if any(chunked.file_hash(p) != h for p,h in identities.items()):
        raise ValueError('Benchmark source changed during execution')
    summary = []
    for condition in config['conditions']:
        selected = [r for r in results if r['condition'] == condition['name']]
        item = dict(condition=condition['name'], hosts=condition['hosts'], workers=sum(condition['hosts'].values()), repeats=len(selected))
        for key in ('launch_to_observed_completion_seconds','worker_timestamp_span_seconds','workers_executing_chunks','largest_worker_peak_RSS_KiB'):
            values=[r[key] for r in selected]
            item[key]=dict(median=statistics.median(values), minimum=min(values), maximum=max(values))
        summary.append(item)
    files = [p for p in raw.rglob('*') if p.is_file()]
    atomic_write_json(out / 'summary.json', dict(status='real_descriptor_concurrency_checked',
        source_sha256=identities, conditions=summary, runs=results,
        run_count=len(results), tasks_per_run=len(expected), serial_discrepancies=0,
        raw_root=str(raw), raw_file_count=len(files), raw_bytes=sum(p.stat().st_size for p in files),
        timing_scope='Launch-to-observed completion includes node preflight, SSH, cold scorer setup, work and <=0.1s controller polling; preparation excluded',
        worker_span_scope='Worker wall timestamps include staggered starts; multi-host clock uncertainty not independently calibrated',
        limits=['Small 837-design library; not million-molecule chemical or quantum throughput',
                'Requested workers and workers doing useful work reported separately',
                'Shared-node measurements, three repeats; min/max ranges are not confidence intervals',
                'Per-worker peaks are not a simultaneous aggregate memory peak'],
        quantum_jobs_submitted=0, quantum_accuracy_established=False))


if __name__ == '__main__':
    main()
