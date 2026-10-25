#!/usr/bin/env python
"""Prepare fresh continuations for proven foreign-residency interruptions.

Only a finite, converged-SCF geometry evaluation or the original declared seed
is reused. Neither a trial frame nor its energy becomes a qualified result.
"""
import argparse
from io import StringIO
import json
import math
from pathlib import Path
import shutil

from redox.core.common import ROOT
from redox.core.io import atomic_write_json
from redox.workflows import polymer_campaign as api


def evaluated_frames(trajectory):
    lines = trajectory.read_text().splitlines(); offset = 0; result = []
    while offset < len(lines):
        size = int(lines[offset])
        end = offset + size + 2
        if end > len(lines):
            break
        result.append('\n'.join(lines[offset:end]) + '\n')
        offset = end
    return result


def prepare(batch, sources):
    from ase.io import read
    batch.mkdir(parents=True, exist_ok=False); entries = []
    owner_check = api.load_module('geometry_interruption_source_owners', ROOT / 'scripts/pipeline/run_polymer_spacing.py')
    audit = api.load_module('geometry_interruption_start_audit', ROOT / 'scripts/validation/audit_polymer_geometries.py')
    worker = ROOT / 'scripts/pipeline/prepare_conformer_geometry_batch.py'
    for source in sources:
        old = api.read_plan(source)
        api.verify_sources(old['producer_sha256']); api.verify_sources(old['input_sha256'])
        status_path, interference_path = source / 'geometry/status.json', source / 'residency_guard/interruption.json'
        report, interference = (json.loads(p.read_text()) for p in (status_path, interference_path))
        if (report['status'] != 'geometry_stage_failed_preserved' or owner_check.owner_alive(report['owner'])
                or interference['pid'] != report['owner']['pid'] or not interference['foreign_resident_PIDs']
                or interference['foreign_processes_signaled']):
            raise ValueError('A stopped, source-bound foreign-residency interruption is required')
        if len(old['tasks']) != 1:
            raise ValueError('One exact geometry task is required')
        task = old['tasks'][0]; key = task['key']; out = batch / key; out.mkdir()
        shutil.copytree(source / 'controls', out / 'controls')
        shutil.copytree(source / 'live_inputs', out / 'live_inputs')
        ready_dir = out / 'live_inputs' / key; xyz = ready_dir / 'input.xyz'
        checkpoint = Path(old['raw_root']) / 'energy' / key / 'rest/geometry_checkpoints'
        trajectory, metadata = checkpoint / 'evaluated_trajectory.xyz', checkpoint / 'evaluations.jsonl'
        selected = dict(kind='original_declared_seed_no_complete_geometry_evaluation')
        witnesses = [source / 'plan.json', status_path, interference_path]
        if trajectory.exists() and metadata.exists():
            frames = evaluated_frames(trajectory)
            records = [json.loads(line) for line in metadata.read_text().splitlines() if line.strip()]
            candidates = [(row['energy_Ha'], i, row) for i, row in enumerate(records[:len(frames)])
                if row['SCF_converged'] and math.isfinite(row['energy_Ha']) and math.isfinite(row['gradient_max_Ha_Bohr'])]
            for _, index, row in sorted(candidates):
                atoms = read(StringIO(frames[index]), format='xyz')
                if len(atoms) != task['atoms']:
                    continue
                xyz.write_text(frames[index])
                result = audit.audit(task, xyz)
                if not result['issues']:
                    selected = dict(kind='lowest_electronic_energy_valid_converged_SCF_evaluation',
                        source=str(checkpoint), evaluation=row, optimizer_acceptance_established=False,
                        native_minimum_or_G_qualified=False)
                    witnesses += [trajectory, metadata]
                    break
            else:
                xyz.write_bytes((source / 'live_inputs' / key / 'input.xyz').read_bytes())
        geometry = audit.audit(task, xyz)
        if geometry['issues']:
            raise ValueError('Continuation starting geometry fails exact graph/stereo sanity')
        atomic_write_json(ready_dir / 'recovery_start_review.json', geometry)
        ready_path = ready_dir / 'task.json'; ready = json.loads(ready_path.read_text())
        ready.update(input_xyz=str(xyz), source_restart_checkpoint=selected)
        ready['input_sha256'][str(xyz)] = api.digest(xyz); atomic_write_json(ready_path, ready)
        inputs = dict(old['input_sha256'])
        inputs.update({str(p): api.digest(p) for p in (*witnesses, xyz, ready_path, ready_dir / 'recovery_start_review.json')})
        producers = dict(old['producer_sha256']); producers[str(Path(__file__).resolve())] = api.digest(__file__)
        for stage in ('energy', 'thermal'):
            for name in ('tasks', 'done', 'claims', 'failures', 'workers', 'pilot', 'launches'):
                (out / stage / name).mkdir(parents=True)
        plan = dict(old, input_sha256=inputs, producer_sha256=producers,
            raw_root=str(ROOT / 'calcs/polymer_models' / batch.name / key),
            recovery_from=str(source), recovery_start=selected,
            purpose='Same-model continuation after confirmed foreign GPU attachment; starting frame remains unqualified')
        atomic_write_json(out / 'plan.json', plan)
        entries.append(dict(key=key, campaign=str(out), source=str(source), plan_sha256=api.digest(out / 'plan.json')))
    atomic_write_json(batch / 'batch.json', dict(status='interrupted_geometry_continuations_prepared',
        entries=entries, worker=str(worker), maximum_RSS_GiB=16, minimum_free_RAM_GiB=32,
        minimum_GPU_MiB=32768, quantum_results_transferred=False,
        source_sha256={str(p): api.digest(p) for p in (Path(__file__).resolve(), worker)}))
    print(json.dumps(dict(prepared=len(entries), qualified_results_transferred=0)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, action='append', required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    prepare(args.out.resolve(), [p.resolve() for p in args.source])


if __name__ == '__main__':
    main()
