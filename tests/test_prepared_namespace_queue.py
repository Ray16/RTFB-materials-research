"""A waiting coordinator preserves active, uncertain and failed scientific work."""
import pytest
from redox.workflows.prepared_namespace_queue import disposition,run


def test_completed_inputs_wait_for_their_active_partner():
    assert disposition([dict(condition='completed_review_pending'),dict(condition='running')])=='observe_or_dispatch'
    assert disposition([dict(condition='completed_review_pending')]*2)=='finished_independent_review_pending'


@pytest.mark.parametrize('status',['scientific_failure_preserved','stopped_attempt_preserved',
    'guard_start_failed_preserved','uncertain_submission_preserved'])
def test_closed_or_uncertain_failure_requires_explicit_review(status):
    assert disposition([dict(condition=status)])=='attempt_preserved_for_explicit_review'
    assert disposition([dict(condition=status),dict(condition='completed_review_pending')])=='attempt_preserved_for_explicit_review'
    assert disposition([dict(condition=status),dict(condition='running')])=='observe_or_dispatch'
    assert disposition([dict(condition=status),dict(condition='pending')])=='observe_or_dispatch'


def test_transient_owner_uncertainty_is_observed_without_restarting():
    assert disposition([dict(condition='ownership_unverified')])=='observe_or_dispatch'
    with pytest.raises(ValueError):disposition([])


@pytest.mark.parametrize('hours',[-1,0,73])
def test_queue_window_is_bounded_before_any_probe_or_mutation(hours):
    with pytest.raises(ValueError):run([],None,None,None,hours,finished_status='done',validate=None,sources=[])
