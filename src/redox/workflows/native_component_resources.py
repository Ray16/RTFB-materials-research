"""Declared native-component resource classes; physical tolerances stay separate.

This module is copied into each new larger pilot's frozen sources so its worker
and independent reviewer use the same bounds without importing another checkout.
"""
import math


LARGE_CLASS = 'double_graft_H100_dimension_pilot_v1'


def retained_estimate(atoms, nao, nocc):
    if any(type(v) is not int or v < 1 for v in (atoms, nao, nocc)) or nocc > nao:
        raise ValueError('Positive native array dimensions required')
    columns = 6
    points = 1 + 4 * columns
    return (3*atoms*3*nao*nocc*8 + points*(2*nao*nao*8+262144)
            + atoms*25000*32 + 5*(3*atoms)**2*8 + 64*2**20)


def larger_profile(atoms, nao, nocc, charge, mult, occupied_sites):
    # These exact dimensions define a prospective resource class, not a claim
    # that the smaller V100 pilot already validated a 137-atom Hessian.
    if ((atoms, nao, nocc, charge, mult) != (137, 1319, 249, 4, 1)
            or len(occupied_sites) != 2 or len(set(occupied_sites)) != 2):
        raise ValueError('The declared double-graft resource class requires its exact dimensions/state')
    estimate = retained_estimate(atoms, nao, nocc)
    budget = 6*2**30
    if estimate*1.2 > budget:
        raise ValueError('Retained arrays exceed the declared output envelope')
    return dict(native_resource_class=LARGE_CLASS, estimated_retained_bytes=estimate,
                retention_margin_factor=1.2, maximum_sampled_columns=6,
                maximum_durable_output_bytes=budget,
                output_profile='double_graft_full_native_6GiB',
                maximum_peak_GPU_MiB=88000, minimum_GPU_memory_MiB=80*1024,
                minimum_GPU_compute_capability=9.0)


def validate_plan(plan):
    if plan.get('native_resource_class') != LARGE_CLASS:
        raise ValueError('Unknown larger native resource class')
    profile = larger_profile(plan['atoms'], plan['response_nao'], plan['response_nocc'],
                             plan['charge'], plan['mult'], plan['prepared_case']['occupied_sites'])
    if (not plan.get('charged_polymer_native_component') or plan.get('matched_empty_polymer_reference')
            or any(plan.get(k) != v for k, v in profile.items())
            or plan.get('maximum_host_RSS_GiB') != 80
            or plan.get('maximum_wall_seconds') != 28800
            or plan.get('minimum_initial_storage_GiB') != 38
            or plan.get('minimum_live_storage_GiB') != 32):
        raise ValueError('The larger native resource bounds changed')
    return profile


def peak_limit(plan):
    if not plan.get('native_resource_class'):
        return 31744
    return validate_plan(plan)['maximum_peak_GPU_MiB']


def check_peak(plan, observed):
    if type(observed) not in (int, float) or not math.isfinite(observed) or not 0 < observed <= peak_limit(plan):
        raise ValueError('Native full-Hessian GPU ceiling failed')
