#!/usr/bin/env python3
"""Sensitivity to missing minima within an independently reviewed finite set.

This is an exhaustive subset diagnostic, not an estimate of unseen conformers,
their probabilities, or ensemble convergence. It does not calculate a cycle.
"""
import argparse
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path

from redox.properties.conformer_ensemble import EV_KJ_MOL, partition

ROOT = Path(__file__).resolve().parents[2]


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def subset_diagnostics(records, cutoff, temperature=298.15):
    if len(records) < 2 or len({row['key'] for row in records}) != len(records):
        raise ValueError('At least two distinct conformer identities required')
    values = [float(row['G_by_cutoff_eV'][cutoff]) for row in records]
    if not all(math.isfinite(x) for x in values):
        raise ValueError('Non-finite free energy')
    full = partition(values, temperature)
    groups = {}
    for size in range(1, len(records)):
        entries = []
        for indices in itertools.combinations(range(len(records)), size):
            selected = [values[i] for i in indices]
            subset = partition(selected, temperature)
            entries.append(dict(keys=[records[i]['key'] for i in indices],
                                subset_minus_full_G_kJ_mol=(subset['G_eV']-full['G_eV'])*EV_KJ_MOL,
                                subset_minus_full_minimum_G_kJ_mol=(min(selected)-min(values))*EV_KJ_MOL))
        groups[str(size)] = dict(combinations=len(entries),
                                 min_partition_difference_kJ_mol=min(e['subset_minus_full_G_kJ_mol'] for e in entries),
                                 max_partition_difference_kJ_mol=max(e['subset_minus_full_G_kJ_mol'] for e in entries),
                                 max_minimum_difference_kJ_mol=max(e['subset_minus_full_minimum_G_kJ_mol'] for e in entries),
                                 subsets=entries)
    return dict(full_G_eV=full['G_eV'], full_minimum_G_eV=full['minimum_G_eV'],
                populations=full['conditional_populations'], subsets_by_size=groups)


def analyze(root, manifest_path, finite_path):
    root = root.resolve()
    checker_path = root/'scripts/validation/summarize_checked_empty_conformers.py'
    spec = importlib.util.spec_from_file_location('checked_empty_review', checker_path)
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
    manifest = json.loads(manifest_path.read_text())
    reviewed = checker.summarize(manifest, root)
    finite = json.loads(finite_path.read_text())
    if finite['status'] != 'finite_checked_empty_conformer_set_not_converged' or finite['manifest_sha256'] != sha(manifest_path):
        raise ValueError('Finite geometry analysis does not match the reviewed inputs')
    for path, digest in finite['source_sha256'].items():
        if sha(path) != digest:
            raise ValueError('Finite-set producer changed')
    records = reviewed['records']
    if set(r['key'] for r in records) != set(finite['keys']):
        raise ValueError('Finite-set identities differ')
    for row in records:
        bound = finite['inputs'][row['key']]
        for prefix in ('geometry', 'electronic_review', 'native_review'):
            if sha(bound[prefix]) != bound[prefix+'_sha256']:
                raise ValueError('Finite-set input changed')
    result = {}
    for rotor in ('25', '50', '100'):
        saved = finite['partitions_by_rmsd_cutoff_A']['0.1'][rotor]
        if sorted(sum(saved['duplicate_groups'], [])) != sorted(row['key'] for row in records) or any(len(g) != 1 for g in saved['duplicate_groups']):
            raise ValueError('Subset calculation requires independently distinct conformers')
        result[rotor] = subset_diagnostics(records, rotor)
        if abs(result[rotor]['full_G_eV']-saved['G_eV']) > 1e-8:
            raise ValueError('Full-set partition does not reproduce upstream result')
    return dict(status='finite_checked_empty_subset_sensitivity_not_converged',
                source_manifest=str(manifest_path), source_manifest_sha256=sha(manifest_path),
                finite_set=str(finite_path), finite_set_sha256=sha(finite_path),
                producer_sha256=sha(__file__), checked_records=[r['key'] for r in records],
                graph_sha256=reviewed['graph_sha256'], method_id=reviewed['method_id'],
                configuration=reviewed['configuration'], solvent=reviewed['solvent'],
                temperature_K=298.15, by_rotor_cutoff_cm=result,
                scope='Exhaustive subsets of four already reviewed structures in the alternate acetonitrile model; not an estimate of unseen conformers or polymer interactions.',
                global_convergence_proven=False, polymer_cycle_computed=False,
                physical_QC_repeated=False, quantum_jobs_submitted=0)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=ROOT)
    p.add_argument('--manifest', type=Path, required=True)
    p.add_argument('--finite', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    result = analyze(a.root, a.manifest.resolve(), a.finite.resolve())
    a.out.parent.mkdir(parents=True, exist_ok=True)
    with a.out.open('x') as stream:
        json.dump(result, stream, indent=2)
        stream.write('\n')


if __name__ == '__main__':
    main()
