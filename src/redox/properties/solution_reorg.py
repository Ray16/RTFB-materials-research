"""Arithmetic and validity contract for solution relaxation descriptors.

Electrostatic/composite continuum relaxation is an electronic descriptor.
It is not an experimentally validated rate, an isolated inner term, or a
finite-temperature activation free energy. Directional values remain available
even when an equal-curvature approximation is unsupported.
"""
from __future__ import annotations
import math


def directional_relaxation(e_initial_eq,e_final_at_initial_eq,e_final_neq,e_final_eq):
    values=[e_initial_eq,e_final_at_initial_eq,e_final_neq,e_final_eq]
    if any(not isinstance(x,(int,float)) or not math.isfinite(x) for x in values):
        raise ValueError('all four solution energies must be finite')
    return dict(total_solution_relaxation_eV=e_final_neq-e_final_eq,
        equilibrium_geometry_component_eV=e_final_at_initial_eq-e_final_eq,
        frozen_slow_increment_eV=e_final_neq-e_final_at_initial_eq,
        vertical_charge_change_eV=e_final_neq-e_initial_eq)


def summarize_relaxations(forward,reverse,*,endpoint_evidence,method_identity,
                          independent_validation=False,implementation_validation=False,equal_curvature_rtol=.05):
    """Preserve directionality and separate numerical from physical evidence.

    A five percent directional tolerance is a stated diagnostic criterion for
    treating the two local surfaces as nearly equal-curvature, not a universal
    physical law or a measured prediction uncertainty. Negative relaxation or
    missing endpoint minimum evidence blocks physical qualification.
    """
    if not method_identity or 'solvent_model' not in method_identity:
        raise ValueError('explicit method and solvent identity required')
    if not 0<equal_curvature_rtol<1:
        raise ValueError('invalid asymmetry threshold')
    totals=[r['total_solution_relaxation_eV'] for r in (forward,reverse)]
    if any(not math.isfinite(x) for x in totals):
        raise ValueError('nonfinite solution relaxation')
    average=sum(totals)/2
    asymmetry=abs(totals[0]-totals[1])/average if average>0 else None
    flags=[]
    if any(x<0 for x in totals):flags.append('negative_relaxation')
    implementation_verified=bool(implementation_validation or independent_validation)
    if not implementation_verified:flags.append('implementation_reference_unverified')
    if len(endpoint_evidence)!=2 or any(e.get('minimum_verified') is not True
                                       for e in endpoint_evidence):
        flags.append('endpoint_minimum_unverified')
    if asymmetry is None or asymmetry>equal_curvature_rtol:
        flags.append('equal_curvature_approximation_unverified')
    if any(e.get('conformer_policy') not in ('same_basin','ensemble_averaged')
           for e in endpoint_evidence):
        flags.append('conformer_ensemble_unresolved')
    return dict(lambda_reduction_eV=totals[0],lambda_oxidation_eV=totals[1],
        symmetric_single_body_approximation_eV=average,
        self_exchange_infinite_separation_eV=sum(totals),
        directional_relative_asymmetry=asymmetry,
        equal_curvature_diagnostic_tolerance=equal_curvature_rtol,
        physical_flags=flags,calculation_verified=bool(independent_validation),
        implementation_reference_verified=implementation_verified,
        endpoint_model_verified=not any(f in flags for f in [
            'negative_relaxation','implementation_reference_unverified',
            'endpoint_minimum_unverified','conformer_ensemble_unresolved']),
        empirical_accuracy_validated=False,production_promotion=False,
        warning='Do not add existing gas inner or outer reorganization terms')
