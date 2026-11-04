#!/usr/bin/env python3
"""Duplicate-aware finite native partitions with explicit nested input sets.

Each graph/configuration is analyzed separately. Full-set duplicate groups are
restricted to the first input batch, so adding a search hit cannot redefine an
earlier well's multiplicity. Search-hit counts are never statistical weights.
"""
import argparse
import importlib.metadata
import hashlib
import json
import math
import os
from pathlib import Path

from redox.core.common import ROOT
from redox.core.hashing import sha256_file
from redox.core.io import atomic_write_json
from redox.properties.conformer_ensemble import batch_change, cluster_minima, summarize
from redox.properties.conformer_geometry import alignment_molecule, symmetry_rmsd_matrix
from redox.validation.archived_sources import ArchivedSources
from redox.validation.namespace_sources import NamespaceSources


def load_record(row, reader):
    report_path = reader.verify(row['review'], row['review_sha256'])
    report = json.loads(report_path.read_text())
    if (report['status'] != 'alternate_native_state_independently_reconstructed'
            or report['primary_candidate_admission'] or report['task']['key'] != row['key']
            or (report['task']['charge'], report['task']['mult']) != (0, 1)
            or report['experimental_model']['solvent'] != 'acetonitrile'
            or report['experimental_method_id'] != row['method_id']
            or row['graph_sha256'] != report['task']['graph_sha256']
            or row['occupied_sites'] != report['task']['occupied_sites']
            or row['backbone_stereo_fingerprint'] != report['task']['backbone_stereo_fingerprint']
            or hashlib.sha256(json.dumps(report['experimental_model'], sort_keys=True).encode()).hexdigest() != row['method_id']):
        raise ValueError('Exact checked neutral native state at one solvent method required')
    sources = report['source_sha256']
    reader.verify_all(sources)

    def bound(path):
        return reader.verify(str(path), sources[str(path)])

    campaign = Path(report['campaign'])
    plan = json.loads(bound(campaign/'plan.json').read_text())
    stage = json.loads(bound(campaign/'experimental_native/status.json').read_text())
    geometry = Path(stage['raw_root'])/'geometry.xyz'
    physical = bound(geometry)
    if (sha256_file(physical) != plan['scored_geometry_sha256'] or plan['tasks'] != [report['task']]
            or stage['status'] != 'experimental_native_finished_independent_review_pending'
            or row['status'] != 'relocated_native_physical_replay_passed'):
        raise ValueError('Scored/native geometry, task or finished physical replay differs')
    parent = json.loads(bound(Path(plan['electronic_review'])).read_text())
    if (not math.isfinite(row['E_eV']) or abs(row['E_eV']-parent['E_smd_eV']) > 1e-9
            or set(row['G_by_cutoff_eV']) != {'25', '50', '100'}
            or any(not math.isfinite(v) or abs(v-report['G_by_cutoff_eV'][k]) > 1e-9
                   for k, v in row['G_by_cutoff_eV'].items())):
        raise ValueError('Conformer energy differs from its exact checked parent')
    molecule = alignment_molecule(report['task']['smiles'], physical)
    return molecule, dict(review=str(report_path), review_sha256=row['review_sha256'],
                          geometry=str(physical), geometry_sha256=sha256_file(physical),
                          graph_sha256=report['task']['graph_sha256'],
                          backbone_stereo_fingerprint=report['task']['backbone_stereo_fingerprint'])


def finite_batches(rows, matrix, first_count):
    """Use one consistent well assignment for each RMSD sensitivity scenario."""
    scenarios = {}
    early = [row for row in rows if int(row['key'].split('conf')[-1]) <= first_count]
    early_keys = {row['key'] for row in early}
    for threshold in (.05, .1, .2):
        groups = cluster_minima(rows, matrix, rmsd_cutoff_A=threshold, energy_tolerance_kJ_mol=.2)
        first_groups = [[key for key in group if key in early_keys] for group in groups]
        first_groups = [group for group in first_groups if group]
        values = {}
        for rotor in ('25', '50', '100'):
            records = [dict(key=r['key'], E_eV=r['E_eV'], G_eV=r['G_by_cutoff_eV'][rotor]) for r in rows]
            first = summarize([r for r in records if r['key'] in early_keys], first_groups,
                              [r['key'] for r in early])
            expanded = summarize(records, groups, [r['key'] for r in rows])
            values[rotor] = dict(first=first, expanded=expanded, comparison=batch_change(first, expanded))
        scenarios[str(threshold)] = values
    return scenarios


def analyze(replays, prior, archive, configurations, sites, first_count, total_count, output):
    if output.exists():
        raise FileExistsError('Preserve previous scientific analyses')
    if not configurations or len(configurations) != len(set(configurations)) or sites != sorted(set(sites)):
        raise ValueError('Distinct configurations and ordered unique occupancy sites required')
    hashes = {str(p): sha256_file(p) for p in replays}
    for path in (Path(__file__).resolve(), ROOT/'src/redox/properties/conformer_geometry.py',
                 ROOT/'src/redox/properties/conformer_ensemble.py', ROOT/'src/redox/build/geometry.py',
                 ROOT/'src/redox/validation/protonation.py'):
        hashes[str(path)] = sha256_file(path)
    selected, readers, inputs = [], [], {}
    seen = set()
    for path in replays:
        data = json.loads(path.read_text())
        if (data['status'] not in ('archived_polymer_native_physical_records_replayed',
                                   'namespace_native_physical_records_replayed')
                or data['states'] != len(data['records']) or data['global_ensemble_convergence_proven']):
            raise ValueError('Completed finite physical replay required')
        namespaces = {}
        for binding in data.get('bindings', []):
            spec_path = Path(binding['namespace'])
            if sha256_file(spec_path) != binding['namespace_sha256']:
                raise ValueError('Native namespace changed after reconstruction')
            hashes[str(spec_path)] = binding['namespace_sha256']
            spec = json.loads(spec_path.read_text())
            namespaces[spec['key']] = spec
        for row in data['records']:
            if row['key'] in seen:
                raise ValueError('Repeated native state across replay inputs')
            seen.add(row['key'])
            if row['configuration'] not in configurations or row['occupied_sites'] != sites:
                continue
            reader = (NamespaceSources(namespaces[row['key']]) if row['key'] in namespaces
                      else ArchivedSources(prior, archive, relative_to_prior=True))
            molecule, binding = load_record(row, reader)
            inputs[row['key']] = binding
            selected.append((row, molecule))
            readers.append(reader)
            atomic_write_json(output.with_name(output.stem+'_progress.json'),
                              dict(checked_records=len(selected), complete=False))
    if len({row['method_id'] for row, _ in selected}) != 1:
        raise ValueError('One declared solvent method per analysis required')
    groups = []
    for configuration in configurations:
        pairs = sorted([p for p in selected if p[0]['configuration'] == configuration],
                       key=lambda p: int(p[0]['key'].split('conf')[-1]))
        rows, molecules = [p[0] for p in pairs], [p[1] for p in pairs]
        label = '_'.join(map(str, sites)) if sites else 'none'
        expected = [f'{configuration}__occ_{label}__conf{i:03d}' for i in range(1, total_count+1)]
        missing = sorted(set(expected)-{r['key'] for r in rows})
        if {r['key'] for r in rows}-set(expected):
            raise ValueError('Unexpected conformer input index')
        row = dict(configuration=configuration, occupied_sites=sites, checked_states=len(rows), missing_keys=missing)
        if missing:
            row.update(status='pending_complete_native_batch', partitions=None)
        else:
            if len({r['graph_sha256'] for r in rows}) != 1:
                raise ValueError('Different exact graphs cannot share a partition')
            matrix, matches = symmetry_rmsd_matrix(molecules)
            row.update(status='finite_native_batches_compared', rmsd_A=matrix.tolist(),
                       symmetry_matches=matches, partitions=finite_batches(rows, matrix, first_count))
        groups.append(row)
    for reader in readers:
        reader.recheck()
    if any(sha256_file(path) != digest for path, digest in hashes.items()):
        raise ValueError('Analysis input or producer changed')
    result = dict(status='source_bound_finite_native_batch_analysis', groups=groups, inputs=inputs,
                  source_sha256=hashes, first_count=first_count, total_count=total_count,
                  versions={p: importlib.metadata.version(p) for p in ('rdkit', 'numpy', 'ase')},
                  duplicate_policy='full-set complete-link groups restricted to each nested input batch',
                  alignment_policy='heavy atoms, chirality-preserving graph maps, proper rotations; no reflections',
                  temperature_K=298.15, solvent='acetonitrile',
                  configurations_are_pooled=False, global_ensemble_convergence_proven=False,
                  quantum_calculations_run=False, primary_candidate_admission=False)
    atomic_write_json(output, result)
    return result


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--replay', type=Path, action='append', required=True)
    p.add_argument('--prior-root', type=Path, default=ROOT)
    p.add_argument('--archive', type=Path, required=True)
    p.add_argument('--configuration', action='append', default=None)
    p.add_argument('--sites', nargs='*', type=int, default=[])
    p.add_argument('--first-count', type=int, default=10)
    p.add_argument('--total-count', type=int, default=20)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    if not 0 < a.first_count < a.total_count:
        p.error('Strictly nested positive input counts required')
    from redox.workflows.resources import allocation_guard, cpu_allocation, host_capacity
    with allocation_guard():
        host_capacity(1, minimum_ram_gib=4)
        os.sched_setaffinity(0, cpu_allocation(1))
    result = analyze([path.resolve() for path in a.replay], a.prior_root.resolve(), a.archive.resolve(),
                     a.configuration or ['all_same', 'alternating'], a.sites, a.first_count, a.total_count,
                     a.out.resolve())
    print(json.dumps(dict(status=result['status'], groups=[dict(configuration=g['configuration'],
                        status=g['status'], checked_states=g['checked_states']) for g in result['groups']])))
