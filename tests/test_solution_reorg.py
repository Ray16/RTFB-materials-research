import pytest

from redox.properties.solution_reorg import directional_relaxation,summarize_relaxations


def test_cds_cancellation_and_directional_solution_bookkeeping():
    # Coordinate-specific CDS shifts affect a directional relaxation but cancel
    # in the two-direction symmetric mean; no extra outer energy may be added.
    f=directional_relaxation(-10.,-10.2,-9.5,-10.5)
    r=directional_relaxation(-10.5,-9.8,-9.,-10.)
    cds0,cds1=.03,.08
    fc=directional_relaxation(-10.+cds0,-10.2+cds0,-9.5+cds0,-10.5+cds1)
    rc=directional_relaxation(-10.5+cds1,-9.8+cds1,-9.+cds1,-10.+cds0)
    assert f['total_solution_relaxation_eV']==pytest.approx(1.)
    assert sum(x['total_solution_relaxation_eV'] for x in (f,r))==pytest.approx(
        sum(x['total_solution_relaxation_eV'] for x in (fc,rc)))


def test_asymmetry_and_endpoint_evidence_cannot_be_hidden_by_average():
    f={'total_solution_relaxation_eV':.4};r={'total_solution_relaxation_eV':1.2}
    evidence=[dict(minimum_verified=False,conformer_policy='different_minima')]*2
    summary=summarize_relaxations(f,r,endpoint_evidence=evidence,
                                 method_identity={'solvent_model':'SMD'})
    assert 'equal_curvature_approximation_unverified' in summary['physical_flags']
    assert not summary['endpoint_model_verified']
    assert not summary['production_promotion']


def test_valid_model_is_still_not_experimentally_validated():
    r={'total_solution_relaxation_eV':.8}
    e=[dict(minimum_verified=True,conformer_policy='same_basin')]*2
    summary=summarize_relaxations(r,r,endpoint_evidence=e,
                                 method_identity={'solvent_model':'SMD'})
    assert summary['endpoint_model_verified']
    assert not summary['empirical_accuracy_validated']
    assert summary['self_exchange_infinite_separation_eV']==1.6


def test_nonfinite_energies_and_negative_relaxations_fail_quality():
    with pytest.raises(ValueError):directional_relaxation(0.,float('nan'),0.,0.)
    e=[dict(minimum_verified=True,conformer_policy='same_basin')]*2
    summary=summarize_relaxations({'total_solution_relaxation_eV':-.1},
        {'total_solution_relaxation_eV':.2},endpoint_evidence=e,
        method_identity={'solvent_model':'SMD'})
    assert 'negative_relaxation' in summary['physical_flags']
    assert not summary['endpoint_model_verified']
