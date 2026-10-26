#!/usr/bin/env python
"""Isolate an exact external control while reusing the normal UMA entry point.

The project library and existing UMA caches are not modified. This adapter runs
only one declared state per reserved process and retains failed attempts.
"""
import argparse
from contextlib import contextmanager
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

from redox.core.common import ROOT
from redox.core.io import atomic_write_json


GROUP = 'bm_ndi_styrene_external'
STATES = ('neu', 'red1', 'red2')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def validate_rows(rows, smiles, count):
    expected = [('neu', 0, 1), ('red1', -1, 2), ('red2', -2, 1)]
    if ([(r['state'], int(r['charge']), int(r['mult_hint'])) for r in rows] != expected
            or any(r['id'] != GROUP or r['smiles'] != smiles or int(r['n_conf']) != count for r in rows)):
        raise ValueError('External-control manifest differs from its exact graph/state declaration')


def prepare(seed_dir, review_path, out):
    seed_dir, review_path, out = Path(seed_dir).resolve(), Path(review_path).resolve(), Path(out).resolve()
    plan_path, manifest_path = seed_dir / 'plan.json', seed_dir / 'manifest.csv'
    seed = json.loads(plan_path.read_text())
    reviewed = json.loads(review_path.read_text())
    if (reviewed['status'] != 'external_neutral_seed_geometry_review_passed_not_quantum'
            or reviewed['plan_sha256'] != digest(plan_path) or reviewed['selected'] != seed['selected']
            or reviewed['reviewer_sha256'] != digest(ROOT / 'scripts/validation/review_external_control_geometries.py')
            or reviewed['producer_sha256'] != digest(ROOT / 'scripts/validation/prepare_external_control_geometries.py')
            or seed['quantum_jobs_submitted'] or seed['charged_state_minima_computed']):
        raise ValueError('External seed independent review is missing or inconsistent')
    if digest(seed['challenge']) != seed['challenge_sha256'] or digest(manifest_path) != seed['proposed_manifest_sha256']:
        raise ValueError('Control source or state manifest changed')
    with manifest_path.open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    validate_rows(rows, seed['source_smiles'], seed['selected'])
    for e in seed['entries']:
        path = (seed_dir / e['path']).resolve()
        if not path.is_relative_to(seed_dir) or digest(path) != e['sha256']:
            raise ValueError('Control seed path or hash changed')
    if out.exists():
        raise FileExistsError('Preserve existing external-control campaign')
    library = out / 'library'
    (library / GROUP / 'conformers').mkdir(parents=True)
    shutil.copyfile(manifest_path, library / 'manifest.csv')
    for e in seed['entries']:
        shutil.copyfile(seed_dir / e['path'], library / GROUP / 'conformers' / Path(e['path']).name)
    source = [Path(__file__).resolve(), ROOT / 'src/redox/qm/uma.py', ROOT / 'src/redox/qm/uma_current.py',
              ROOT / 'src/redox/qm/uma_nitro.py', ROOT / 'src/redox/qm/uma_model_identity.py',
              ROOT / 'src/redox/core/common.py', ROOT / 'src/redox/core/io.py',
              ROOT / 'src/redox/build/geometry.py', ROOT / 'src/redox/workflows/resources.py',
              ROOT / 'src/redox/validation/protonation.py', ROOT / 'config/uma_execution.json',
              ROOT / 'scripts/validation/review_external_control_geometries.py',
              ROOT / 'scripts/validation/prepare_external_control_geometries.py', plan_path, review_path]
    staged = [library / 'manifest.csv', *sorted((library / GROUP / 'conformers').glob('conf_*.xyz'))]
    if len(staged) != seed['selected'] + 1:
        raise ValueError('Staged seed inventory differs')
    plan = dict(status='external_control_UMA_campaign_prepared_not_started', group=GROUP, states=list(STATES),
                created=time.time(), seed_plan=str(plan_path), independent_review=str(review_path),
                library=str(library), output_root=str(out / 'output'),
                source_sha256={str(p):digest(p) for p in source},
                staged_sha256={str(p):digest(p) for p in staged},
                role='Exact free-molecule method control; no synthesis admission',
                quantum_jobs_submitted=0, DFT_or_experimental_accuracy_inferred=False)
    atomic_write_json(out / 'plan.json', plan)
    return plan


def validate_campaign(campaign):
    plan = json.loads((campaign / 'plan.json').read_text())
    if plan['status'] != 'external_control_UMA_campaign_prepared_not_started' or plan['group'] != GROUP or plan['states'] != list(STATES):
        raise ValueError('Unexpected external-control plan')
    for path, expected in {**plan['source_sha256'], **plan['staged_sha256']}.items():
        if digest(path) != expected:
            raise ValueError('External-control source or staged input changed: ' + path)
    library = Path(plan['library']).resolve()
    if not library.is_relative_to(campaign.resolve()) or Path(plan['output_root']).resolve() != campaign.resolve() / 'output':
        raise ValueError('Campaign paths must remain isolated')
    with (library / 'manifest.csv').open(newline='') as stream:
        rows = list(csv.DictReader(stream))
    seed = json.loads(Path(plan['seed_plan']).read_text())
    validate_rows(rows, seed['source_smiles'], seed['selected'])
    actual = sorted((library / GROUP / 'conformers').glob('conf_*.xyz'))
    if {str(p) for p in actual} != set(plan['staged_sha256']) - {str(library / 'manifest.csv')}:
        raise ValueError('Staged seed inventory changed')
    return plan, rows


@contextmanager
def isolated_library(uma, library, rows):
    previous = uma.LIBRARY, uma.read_manifest
    uma.LIBRARY, uma.read_manifest = library, lambda: rows
    try:
        yield
    finally:
        uma.LIBRARY, uma.read_manifest = previous


def worker(campaign, state, require_nitro=False):
    from redox.workflows.resources import verify_gate_owner
    from redox.qm import uma_current
    if state not in STATES:
        raise ValueError('Unrecognized external-control state')
    plan, rows = validate_campaign(campaign)
    verify_gate_owner()
    gpu = os.environ['CUDA_VISIBLE_DEVICES']
    raw = subprocess.check_output(['nvidia-smi', '--id='+gpu, '--query-gpu=compute_cap',
                                   '--format=csv,noheader,nounits'], text=True, timeout=15).strip()
    profile = uma_current.select_backend(tuple(map(int, raw.split('.'))), require_nitro)
    policy = json.loads((ROOT / 'config/uma_execution.json').read_text())
    python = policy['nitro_python'] if profile == 'nitro' else policy['legacy_python']
    wanted = Path(python).parent.parent
    needs_lib = profile == 'nitro' and not os.environ.get('LD_LIBRARY_PATH', '').startswith(policy['nitro_library_path'])
    if Path(sys.prefix).resolve() != wanted.resolve() or needs_lib:
        raise RuntimeError('Select the matching interpreter and library path before reserving the GPU; do not replace a live guard process')
    os.environ['FAIRCHEM_CACHE_DIR'] = policy['fairchem_cache']
    failure = campaign / 'failures' / (state + '.json')
    if failure.exists():
        raise RuntimeError('Previous scientific failure is preserved for review')
    claim = campaign / 'claims' / state
    claim.parent.mkdir(parents=True, exist_ok=True)
    claim.mkdir(exist_ok=False)
    owner = dict(host=os.uname().nodename.split('.')[0], pid=os.getpid(),
                 process_start=Path('/proc/self/stat').read_text().rsplit(')', 1)[1].split()[19],
                 state=state, profile=profile, started=time.time())
    atomic_write_json(claim / 'owner.json', owner)
    policy = dict(policy, nitro_output=str(campaign / 'output/nitro'), legacy_output=str(campaign / 'output/legacy'))
    from redox.qm import uma
    try:
        with isolated_library(uma, Path(plan['library']), rows):
            uma_current.run_uma(['--only', GROUP + ':' + state, '--device', 'cuda'], profile, policy)
        result = campaign / 'output' / profile / GROUP / state / 'result.json'
        if not result.is_file():
            raise RuntimeError('UMA returned without the expected state result')
        receipt = dict(status='UMA_control_worker_returned_review_required', owner=owner,
                       result=str(result), result_sha256=digest(result), finished=time.time(),
                       solution_G_or_experimental_accuracy_inferred=False)
        atomic_write_json(campaign / 'completed' / (state + '.json'), receipt)
    except BaseException as error:
        atomic_write_json(failure, dict(status='scientific_attempt_failed_preserved', owner=owner,
                                       error=repr(error), finished=time.time()))
        raise
    finally:
        # A hard interruption leaves the claim for explicit owner verification.
        if (claim / 'owner.json').exists():
            (claim / 'owner.json').unlink()
        claim.rmdir()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('action', choices=('prepare', 'worker'))
    p.add_argument('--seed-dir', type=Path)
    p.add_argument('--review', type=Path)
    p.add_argument('--campaign', type=Path, required=True)
    p.add_argument('--state', choices=STATES)
    p.add_argument('--require-nitro', action='store_true')
    args = p.parse_args()
    campaign = args.campaign.resolve()
    if args.action == 'prepare':
        if args.seed_dir is None or args.review is None:
            p.error('prepare requires --seed-dir and --review')
        result = prepare(args.seed_dir, args.review, campaign)
        print(json.dumps({'status':result['status'], 'states':len(result['states'])}))
    else:
        if args.state is None:
            p.error('worker requires --state')
        worker(campaign, args.state, args.require_nitro)


if __name__ == '__main__':
    main()
