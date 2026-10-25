from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


def test_extended_component_checks_must_finish_before_geometry_admission():
    module = api.load_module('geometry_repair_provisional_gate_test', ROOT / 'scripts/pipeline/repair_geometry_energy_consistency.py')
    original = dict(status='stationary_geometry_checked_active_parent_and_native_minimum_pending',
        conservative_total_force_check={'passed': True})
    pending = module.provisional(original)
    assert pending['status'] == 'base_geometry_checks_passed_extended_CDS_checks_pending'
    assert original['status'] == 'stationary_geometry_checked_active_parent_and_native_minimum_pending'


def test_existing_failure_cannot_be_promoted_by_repair_wrapper():
    module = api.load_module('geometry_repair_failure_test', ROOT / 'scripts/pipeline/repair_geometry_energy_consistency.py')
    assert module.provisional({'status': 'geometry_stage_failed_preserved'})['status'] == 'geometry_stage_failed_preserved'
