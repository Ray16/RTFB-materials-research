"""Reconstruct the retained-state SI quantities from frozen replay records."""
from collections import Counter
import hashlib
import json
import math
import xml.etree.ElementTree as ET


def audit(folder):
    load = lambda name: json.loads((folder/name).read_text())
    data = load('polymer_native_replay.json')
    summary = load('polymer_native_replay_summary.json')
    figure = load('polymer_native_sampling_figure.json')
    rows = data['records']
    expected = {('all_same', ()): 20, ('alternating', ()): 16,
                ('all_same', (2,)): 14, ('alternating', (2,)): 10}
    if (data['status'] != 'archived_polymer_native_physical_records_replayed'
            or data['states'] != len(rows) or len({r['key'] for r in rows}) != 60
            or Counter((r['configuration'], tuple(r['occupied_sites'])) for r in rows) != expected
            or data['sources_verified'] != 12386
            or any(data[k] for k in ('quantum_calculations_run', 'current_production_runtime_qualified',
                                    'global_ensemble_convergence_proven'))):
        raise ValueError('Retained-state coverage or scope differs')
    for row in rows:
        if (row['status'] != 'relocated_native_physical_replay_passed'
                or not math.isfinite(row['minimum_frequency_cm']) or row['minimum_frequency_cm'] < -5
                or not 0 <= row['frequency_error_cm'] <= 1e-5
                or not 0 <= row['maximum_column_error'] <= .01
                or set(row['G_by_cutoff_eV']) != {'25', '50', '100'}
                or any(not math.isfinite(v) for v in row['G_by_cutoff_eV'].values())
                or any(not 0 <= v <= 1e-8 for v in row['nuclear_errors_eV'].values())):
            raise ValueError('Retained physical replay fails')
    conv = 1.602176634e-19*6.02214076e23/1000
    panels = {}
    if len(figure['panels']) != 4 or len({r['method_id'] for r in rows}) != 1:
        raise ValueError('Four panels and a single declared solvent method required')
    for signature in expected:
        selected = [r for r in rows if (r['configuration'], tuple(r['occupied_sites'])) == signature]
        if len({r['graph_sha256'] for r in selected}) != 1:
            raise ValueError('Different graphs mixed')
        indices = [int(r['key'].split('conf')[-1]) for r in selected]
        if len(set(indices)) != len(indices) or not set(indices) <= set(range(1, 21)):
            raise ValueError('Distinct declared indices required')
        values = [r['G_by_cutoff_eV']['50'] for r in selected]
        low = min(values)
        early = [v for i, v in zip(indices, values) if i <= 10]
        panels[signature] = dict(checked=len(selected), missing_indices=sorted(set(range(1, 21))-set(indices)),
                                observed_range_kJ_mol=(max(values)-low)*conv,
                                minimum_change_kJ_mol=(low-min(early))*conv)
    for record in figure['panels']:
        actual = panels[record['configuration'], tuple(record['occupied_sites'])]
        if (actual['checked'] != record['checked'] or actual['missing_indices'] != record['missing_indices']
                or abs(actual['observed_range_kJ_mol']-record['observed_range_kJ_mol']) > 1e-8
                or abs(actual['minimum_change_kJ_mol']-
                       record['lowest_checked_minus_lowest_checked_in_first_ten_kJ_mol']) > 1e-8):
            raise ValueError('Sampling-figure numbers do not reconstruct')
    kt = 1.380649e-23/1.602176634e-19*298.15
    empty = sorted([r for r in rows if r['configuration'] == 'all_same' and not r['occupied_sites']],
                   key=lambda r: int(r['key'].split('conf')[-1]))
    changes = {}
    for rotor in ('25', '50', '100'):
        partitions = []
        for n in (10, 20):
            values = [r['G_by_cutoff_eV'][rotor] for r in empty[:n]]
            low = min(values)
            partitions.append(low-kt*math.log(math.fsum(math.exp((low-v)/kt) for v in values)))
        change = (partitions[1]-partitions[0])*conv
        if abs(change-summary['all_same_empty_batch_ensemble_G_change_kJ_mol'][rotor]) > 1e-8:
            raise ValueError('Finite partition differs')
        changes[rotor] = change
    figure_hash = hashlib.sha256((folder/'polymer_native_sampling.pdf').read_bytes()).hexdigest()
    if figure_hash != next(v for p, v in figure['outputs_sha256'].items() if p.endswith('.pdf')):
        raise ValueError('Figure changed')
    if data['full_report_sha256'] != summary['raw_report_sha256']:
        raise ValueError('Raw report binding differs')
    suites = list(ET.parse(folder/'polymer_native_replay_tests.xml').getroot().iter('testsuite'))
    if [sum(int(s.get(k, '0')) for s in suites) for k in ('tests', 'failures', 'errors', 'skipped')] != [17, 0, 0, 0]:
        raise ValueError('Replay tests incomplete')
    return dict(retained_states=60, source_bindings=12386, sampling_panels=4,
                empty_ensemble_change_kJ_mol=changes,
                occupied_observed_minimum_change_kJ_mol=panels['all_same', (2,)]['minimum_change_kJ_mol'],
                complete_occupied_batch_convergence_proven=False)
