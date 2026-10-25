import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('first_site_claims', ROOT / 'manuscript/check_claims.py')
claims = importlib.util.module_from_spec(spec); spec.loader.exec_module(claims)


def load(name):
    path = (ROOT / 'results/reviews/2026-10-25/pmdi_first_site_grafting_checked/summary.json'
        if name == 'pmdi_first_site_grafting.json' else
        ROOT / 'results/polymer_models/campaigns/2026-10-24_GPU_electronic/baseline__design_pmdi_ethyl__occ_2_3/energy/tasks/baseline__design_pmdi_ethyl__occ_2_3__pose0/status.json')
    return json.loads(path.read_text())


def test_first_site_reactions_and_rank_sensitivity_reconstruct():
    report = claims.audit_first_site_grafting(load)
    assert report['checked_first_site_reactions'] == 4
    assert 5.3 < report['sites']['2']['polymer_minus_capped_kJ_mol'] < 5.4
    assert report['placement_deltaG_by_cutoff_kJ_mol']['25'] > 0
    assert report['placement_deltaG_by_cutoff_kJ_mol']['50'] < 0
    assert not report['paired_native_G_available']


@pytest.mark.parametrize('change', ['base_term', 'relative_penalty', 'missing_pair_filled', 'identity', 'balance', 'unqualified_pair'])
def test_changed_or_promoted_first_site_claims_fail(change):
    def changed(name):
        data = load(name)
        if name == 'pmdi_first_site_grafting.json':
            if change == 'base_term': data['reagent_states']['base']['G_by_cutoff_eV']['50'] += .1
            elif change == 'relative_penalty': data['singles'][0]['by_cutoff']['50']['polymer_minus_capped_kJ_mol'] += 1
            elif change == 'missing_pair_filled': data['paired_G_eV'] = -1
            elif change == 'identity': data['candidate_id'] = 'other_design'
            elif change == 'balance': data['singles'][0]['balances']['free_ions']['atom_balanced'] = False
        elif change == 'unqualified_pair': data['physical_QC'] = 'pending'
        return data
    with pytest.raises(ValueError): claims.audit_first_site_grafting(changed)
