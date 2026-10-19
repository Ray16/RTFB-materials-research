"""Register only the four unambiguous parent structures; freeze this panel's inputs."""
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time

from redox.core.common import ROOT, DFT, load_config, read_manifest
from redox.core.io import atomic_write_json
from redox.build import validation as build
from redox.build.manifest import locked_manifest
from redox.qm.native_thermal import NATIVE_THERMAL
from redox.core.protocol import ACTIVE_SP
from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity

base = Path(__file__).resolve().parent
assert not (base / 'plan.json').exists(), 'Prepared panel already exists'
cores = [b for b in load_config('benchmark_mecn').BENCHMARK
         if b['events'][0]['doi'] == '10.1021/acs.chemmater.8b04770'
         and any(e['kind'] == 'E' and e['tier'] == 'A' for e in b['events'])]
assert len(cores) == 4
for core in cores:
    assert not (DFT / core['id']).exists() and not (ROOT / 'library' / core['id']).exists()
assert len({row['smiles'] for row in cores}) == 4
with allocation_guard():
    capacity = host_capacity(1)
    cpus = cpu_allocation(1)
    os.sched_setaffinity(0, cpus)
with locked_manifest(build.MANIFEST):
    shutil.copy2(build.MANIFEST, base / 'manifest_before.csv')
    build._build_selected(cores)
rows = [row for row in read_manifest() if row['id'] in {c['id'] for c in cores}]
assert len(rows) == 12
tasks = []
for row in rows:
    inputs = sorted((ROOT / 'library' / row['id'] / 'conformers').glob('conf_*.xyz'))
    assert inputs
    tasks.append(dict(key=row['id'] + '__' + row['state'], id=row['id'], state=row['state'],
                      manifest_row=row,
                      input_sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs}))
sources = set((ROOT / 'src/redox').rglob('*.py'))
sources.update((ROOT / 'config').glob('*.py'))
sources.update((ROOT / 'data/raw/validation/two_wave_mecn').glob('*.csv'))
sources.update([base / 'worker.py', Path(__file__).resolve()])
plan = dict(schema=1, status='prepared_not_launched', created=time.time(),
            git_head=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
            parent_count=4, state_count=12, strict_target_potentials=5,
            irreversible_second_waves_are_predictions_only=True,
            source_conflicted_parent_not_scheduled='B(MEEO)EPT',
            protocol=ACTIVE_SP, thermal_protocol=NATIVE_THERMAL,
            optimization_seeds='Up to three FF-ranked ETKDG seeds; charge-specific xTB-ALPB warm start; full-response DFT/SMD optimization',
            seed=0xC0FFEE, native_minimum_repair_scheduled=False,
            attempt_limit_per_state=1, canonical_outputs_must_be_absent_at_start=True,
            resources=dict(preparation_cpu=cpus, host_capacity=capacity, worker_threads=2,
                           minimum_GPU_memory_MiB=31 * 1024, native_memory_cap_GiB=6),
            independent_admission_required_before_benchmark_update=True,
            tasks=tasks, producer_sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                         for p in sorted(sources)})
atomic_write_json(base / 'plan.json', plan)
print('prepared', len(tasks), 'new parent/state jobs')
