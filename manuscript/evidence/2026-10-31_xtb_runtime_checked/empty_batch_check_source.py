"""Independently reconstruct the two empty-backbone sampling comparisons."""
from collections import Counter
import hashlib
import json
import math

CONV = 1.602176634e-19*6.02214076e23/1000
KT = 1.380649e-23/1.602176634e-19*298.15


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(actual, expected):
    if not math.isfinite(actual) or abs(actual-expected) > 1e-8:
        raise ValueError('Frozen empty-backbone quantity does not reconstruct')


def partition(rows, groups, rotor):
    lookup = {r['key']: r for r in rows}
    representatives = [min(g, key=lambda k: (lookup[k]['G_by_cutoff_eV'][rotor], k)) for g in groups]
    values = [lookup[k]['G_by_cutoff_eV'][rotor] for k in representatives]
    low = min(values)
    weights = [math.exp((low-value)/KT) for value in values]
    total = math.fsum(weights)
    return dict(G_eV=low-KT*math.log(total), minimum_G_eV=low,
                minimum_E_eV=min(r['E_eV'] for r in rows), representative_keys=representatives,
                conditional_populations=[value/total for value in weights], unique_minima=len(groups))


def audit(folder):
    load = lambda name: json.loads((folder/name).read_text())
    names = ('polymer_native_replay.json', 'polymer_native_continuations.json',
             'native_addition_five.json', 'native_addition_three.json')
    datasets = [load(name) for name in names]
    rows = [r for data in datasets for r in data['records']]
    expected = {('all_same', ()): 20, ('alternating', ()): 20,
                ('all_same', (2,)): 18, ('alternating', (2,)): 16, ('all_same', (4,)): 1}
    if (len(rows) != 75 or len({r['key'] for r in rows}) != 75
            or Counter((r['configuration'], tuple(r['occupied_sites'])) for r in rows) != expected
            or len({r['method_id'] for r in rows}) != 1):
        raise ValueError('Current native evidence identities or coverage differ')
    for row in rows:
        if (row['status'] != 'relocated_native_physical_replay_passed'
                or not math.isfinite(row['E_eV'])
                or set(row['G_by_cutoff_eV']) != {'25', '50', '100'}
                or any(not math.isfinite(v) for v in row['G_by_cutoff_eV'].values())
                or not math.isfinite(row['minimum_frequency_cm']) or row['minimum_frequency_cm'] < -5
                or not 0 <= row['maximum_column_error'] <= .01
                or not 0 <= row['frequency_error_cm'] <= 1e-5):
            raise ValueError('Frozen native physical reconstruction fails')
    analysis = load('empty_batch_analysis.json')
    if (analysis['status'] != 'source_bound_finite_native_batch_analysis'
            or analysis['global_ensemble_convergence_proven'] or analysis['configurations_are_pooled']
            or analysis['first_count'] != 10 or analysis['total_count'] != 20
            or analysis['temperature_K'] != 298.15 or analysis['solvent'] != 'acetonitrile'
            or len(analysis['groups']) != 2
            or {g['configuration'] for g in analysis['groups']} != {'all_same', 'alternating'}):
        raise ValueError('Finite separate-backbone scope differs')
    source_hashes = set(analysis['source_sha256'].values())
    if datasets[0]['full_report_sha256'] not in source_hashes or any(sha(folder/name) not in source_hashes for name in names[1:]):
        raise ValueError('Batch analysis references different native replays')
    results, macros, unique = {}, [], {}
    for group in analysis['groups']:
        selected = sorted([r for r in rows if r['configuration'] == group['configuration'] and not r['occupied_sites']],
                          key=lambda r: r['key'])
        keys = [r['key'] for r in selected]
        if any(analysis['inputs'][r['key']]['review_sha256'] != r['review_sha256']
               or analysis['inputs'][r['key']]['graph_sha256'] != r['graph_sha256'] for r in selected):
            raise ValueError('Geometry analysis has a different checked native parent')
        first = selected[:10]
        first_keys = set(keys[:10])
        if group['occupied_sites'] or group['checked_states'] != 20 or group['missing_keys']:
            raise ValueError('Complete empty-backbone batches required')
        matrix = group['rmsd_A']
        if len(matrix) != 20 or any(len(r) != 20 for r in matrix):
            raise ValueError('Full pairwise geometry matrix required')
        for i in range(20):
            for j in range(20):
                if not math.isfinite(matrix[i][j]) or matrix[i][j] < 0 or abs(matrix[i][j]-matrix[j][i]) > 1e-10:
                    raise ValueError('Invalid pairwise geometry distances')
            close(matrix[i][i], 0)
        for threshold in ('0.05', '0.1', '0.2'):
            groups = []
            for index in sorted(range(20), key=lambda i: (selected[i]['E_eV'], keys[i])):
                eligible = next((g for g in groups if all(matrix[index][j] <= float(threshold)
                        and abs(selected[index]['E_eV']-selected[j]['E_eV'])*CONV <= .2 for j in g)), None)
                if eligible is None:
                    groups.append([index])
                else:
                    eligible.append(index)
            groups = [[keys[i] for i in g] for g in groups]
            early_groups = [[k for k in g if k in first_keys] for g in groups]
            early_groups = [g for g in early_groups if g]
            for rotor in ('25', '50', '100'):
                value = group['partitions'][threshold][rotor]
                a, b = partition(first, early_groups, rotor), partition(selected, groups, rotor)
                for label, expected_partition, expected_groups in (('first', a, early_groups), ('expanded', b, groups)):
                    actual = value[label]
                    if (actual['duplicate_groups'] != expected_groups or actual['representative_keys'] != expected_partition['representative_keys']
                            or actual['unique_minima'] != expected_partition['unique_minima'] or actual['search_hits_are_weights']
                            or len(actual['conditional_populations']) != len(expected_groups)):
                        raise ValueError('Duplicate wells or their statistical weights differ')
                    for field in ('G_eV', 'minimum_G_eV', 'minimum_E_eV'):
                        close(actual[field], expected_partition[field])
                    for actual_weight, expected_weight in zip(actual['conditional_populations'], expected_partition['conditional_populations']):
                        close(actual_weight, expected_weight)
                comparison = value['comparison']
                for field, term in (('minimum_E_change_kJ_mol', 'minimum_E_eV'), ('minimum_G_change_kJ_mol', 'minimum_G_eV'),
                                    ('ensemble_G_change_kJ_mol', 'G_eV')):
                    close(comparison[field], (b[term]-a[term])*CONV)
                if (comparison['tolerance_kJ_mol'] != 2 or not comparison['within_batch_tolerance']
                        or comparison['global_ensemble_convergence_proven']
                        or any(abs(comparison[k]) > 2 for k in ('minimum_E_change_kJ_mol', 'minimum_G_change_kJ_mol', 'ensemble_G_change_kJ_mol'))):
                    raise ValueError('Declared batch diagnostic differs')
        nominal = group['partitions']['0.1']['50']
        prefix = 'EmptyAllSame' if group['configuration'] == 'all_same' else 'EmptyAlternating'
        macros.extend(['\\newcommand{\\'+prefix+'Change}{'+format(nominal['comparison']['ensemble_G_change_kJ_mol'], '.2f')+'}',
                       '\\newcommand{\\'+prefix+'Unique}{'+str(nominal['expanded']['unique_minima'])+'}'])
        results[group['configuration']] = nominal['comparison']['ensemble_G_change_kJ_mol']
        unique[group['configuration']] = nominal['expanded']['unique_minima']
    if (folder/'empty_batch_values.tex').read_text() != '\n'.join(macros)+'\n':
        raise ValueError('Manuscript numerical macros differ from reconstructed quantities')
    figure = load('empty_batch_figure.json')
    if (figure['analysis_sha256'] != sha(folder/'empty_batch_analysis.json')
            or sha(folder/'empty_batch_stability.pdf') != next(v for p, v in figure['outputs_sha256'].items() if p.endswith('.pdf'))
            or figure['settings_ranges_are_confidence_intervals'] or figure['global_ensemble_convergence_proven']):
        raise ValueError('Frozen figure binding or interpretation differs')
    if len(figure['records']) != 2 or unique != {'all_same': 20, 'alternating': 19}:
        raise ValueError('Published distinct-minimum counts differ')
    for record in figure['records']:
        group = next(g for g in analysis['groups'] if g['configuration'] == record['configuration'])
        if record['unique_minima'] != unique[record['configuration']]:
            raise ValueError('Figure minimum count differs')
        for metric in ('minimum_E_change_kJ_mol', 'minimum_G_change_kJ_mol', 'ensemble_G_change_kJ_mol'):
            close(record['nominal'][metric], group['partitions']['0.1']['50']['comparison'][metric])
            values = [v['comparison'][metric] for cutoff in group['partitions'].values() for v in cutoff.values()]
            for actual, expected_range in zip(record['settings_ranges_kJ_mol'][metric], (min(values), max(values))):
                close(actual, expected_range)
    return dict(checked_native_states=75, empty_states=40, configurations=2,
                unique_minima=unique,
                nominal_ensemble_change_kJ_mol=results, all_settings_within_batch_criterion=True,
                global_ensemble_convergence_proven=False, occupied_ensemble_cycles_complete=False)
