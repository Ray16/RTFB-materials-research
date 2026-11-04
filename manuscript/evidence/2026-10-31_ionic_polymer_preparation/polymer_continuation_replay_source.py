#!/usr/bin/env python3
"""Reconstruct completed continuation Hessians and thermochemistry on one CPU.

The independent physical equations are shared with the retained-state replay.
Namespace aliases relocate paths only; all underlying source hashes and the
normal electronic, stationary-geometry, Hessian and ownership gates still apply.
"""
import argparse
import importlib.util
import json
import os
from pathlib import Path

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.validation.namespace_sources import NamespaceSources, read_checked_review


def load(path, name):
    specification = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def run(manifest_path, output):
    if output.exists():
        raise FileExistsError('Preserve earlier replay evidence')
    manifest_hash = sha256_file(manifest_path)
    manifest = json.loads(manifest_path.read_text())
    specs = [Path(path) for path in manifest['native_namespaces']]
    if not specs or len(specs) != len(set(specs)):
        raise ValueError('Nonempty distinct namespace inputs required')
    paths = [Path(__file__).resolve(), ROOT/'scripts/validation/replay_archived_polymer_native.py',
             ROOT/'scripts/validation/audit/native_partition.py',
             ROOT/'src/redox/validation/namespace_sources.py',
             ROOT/'src/redox/validation/archived_sources.py']
    producers = {str(path): sha256_file(path) for path in paths}
    physical = load(paths[1], 'continuation_physical_replay')
    equations = load(paths[2], 'continuation_nuclear_equations')
    rows, readers, bindings = [], [], []
    for path in specs:
        report, binding = read_checked_review(path, 'native')
        spec = json.loads(path.read_text())
        reader = NamespaceSources(spec)
        # The independent equations must agree with the frozen source used by
        # this exact state's original audit, even after checkout relocation.
        original = str(Path(spec['prior_root'])/'scripts/validation/audit/native_partition.py')
        if sha256_file(reader.path(original)) != producers[str(paths[2])]:
            raise ValueError('Independent nuclear equations differ from archived source')
        item = dict(key=spec['key'], review=binding['review'], review_sha256=binding['review_sha256'])
        row = physical.replay(item, reader, equations.nuclear_equations)
        if row['key'] in {r['key'] for r in rows}:
            raise ValueError('Repeated chemical/conformer identity across namespaces')
        rows.append(row)
        readers.append(reader)
        bindings.append(binding)
        atomic_write_json(output.with_name(output.stem+'_progress.json'),
                          dict(status='partial_replay_final_recheck_pending', records=rows,
                               expected=len(specs), completed=len(rows), manifest_sha256=manifest_hash))
    for reader in readers:
        reader.recheck()
    for binding in bindings:
        if any(sha256_file(binding[k]) != binding[k+'_sha256'] for k in ('namespace', 'review')):
            raise ValueError('Continuation namespace or review changed during reconstruction')
    if (sha256_file(manifest_path) != manifest_hash
            or any(sha256_file(path) != digest for path, digest in producers.items())):
        raise ValueError('Replay input or producer changed')
    result = dict(status='namespace_native_physical_records_replayed', states=len(rows), records=rows,
                  source_bindings=sum(len(reader.receipts) for reader in readers), bindings=bindings,
                  manifest_sha256=manifest_hash, producer_sha256=producers,
                  CPU_affinity=sorted(os.sched_getaffinity(0)), quantum_calculations_run=False,
                  primary_candidate_admission=False, global_ensemble_convergence_proven=False)
    atomic_write_json(output, result)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4)
        os.sched_setaffinity(0, cpu_allocation(1))
    result = run(a.manifest.resolve(), a.out.resolve())
    print(json.dumps({k: result[k] for k in ('status', 'states', 'source_bindings')}))


if __name__ == '__main__':
    main()
