"""Paired potential shifts from admitted candidate snapshots.

The collector supplies field-specific thermal, energy and spin admission. A pair
must describe the same charge/electron ladder. Missing events stay missing;
surface estimates and parent measurements never substitute for admitted values.
"""
import json
import math


def _finite(value):
    return type(value) in (int, float) and math.isfinite(value)


def _ladder(row):
    states = row.get('states')
    if isinstance(states, str):
        states = json.loads(states)
    if not states:
        raise ValueError('Missing charge/electron ladder')
    ladder = sorted((int(s[1]), int(s[3])) for s in states)
    if len({q for q, _ in ladder}) != len(ladder):
        raise ValueError('Ambiguous charge ladder')
    return ladder


def paired_potential_shifts(rows):
    """Compare declared controls within one audited snapshot, without ranking."""
    by_id = {}
    for row in rows:
        key = row.get('registry_id') or row['id']
        if key in by_id:
            raise ValueError('Duplicate comparison identity: ' + key)
        by_id[key] = row
    comparisons = []
    for row in rows:
        control_id = row.get('comparison_baseline')
        if not control_id:
            continue
        control = by_id.get(control_id)
        result = dict(id=row['id'], core=row.get('search_core'),
                      contrast=row.get('search_contrast'), baseline=control_id,
                      delta_E1_V=None, delta_E2_V=None, delta_wave_spacing_V=None,
                      status='pending control', experimental_improvement_verified=False)
        if control:
            if _ladder(row) != _ladder(control):
                raise ValueError('Incompatible comparison ladders: ' + row['id'])
            for source, target in (('E1_vs_Fc_V', 'delta_E1_V'),
                                   ('E2_vs_Fc_V', 'delta_E2_V'),
                                   ('wave_spacing_V', 'delta_wave_spacing_V')):
                if _finite(row.get(source)) and _finite(control.get(source)):
                    result[target] = row[source] - control[source]
            count = sum(result[k] is not None for k in ('delta_E1_V', 'delta_E2_V'))
            result['status'] = ('both potential shifts admitted' if count == 2 else
                                'one potential shift admitted' if count == 1 else
                                'pending admitted event pair')
        comparisons.append(result)
    return comparisons
