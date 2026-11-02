import copy
import pytest
from redox.workflows.sampling_energy import profile, validate_inventory


def roster():
    rows = []
    for cfg in ('all_same', 'alternating'):
        for sites in ([], [2], [3], [4], [2, 3], [2, 4]):
            for start in range(20):
                key = cfg + '__' + str(sites) + '__' + str(start)
                task = dict(key=key, occupied_sites=sites, stereoisomer=cfg, charge=0, mult=1,
                    candidate_id='exp_phthalimide_control', role='interaction_method_control; not synthesis admission')
                rows.append(dict(key=key, task=task))
    return dict(states=rows, status='complete_corrected_sampling_start_inventory', declared_inputs=240,
        selected_inputs=240, missing=[], historical_values_transferred=False, experimental_occupancy_weights=False)


def test_all_twelve_declared_cells_and_size_profiles_are_retained():
    data = roster(); assert len(validate_inventory(data)) == 240
    for row in data['states']:
        assert profile(row['task'])['maximum_host_RSS_GiB'] == (240 if len(row['task']['occupied_sites']) == 2 else 160)


@pytest.mark.parametrize('bad', ['missing_cell', 'duplicate', 'charge', 'admission', 'weights', 'old_values'])
def test_partial_or_relabelled_roster_cannot_be_released(bad):
    data = copy.deepcopy(roster())
    if bad == 'missing_cell': data['states'][0]['task']['occupied_sites'] = [4]
    if bad == 'duplicate': data['states'][0] = data['states'][1]
    if bad == 'charge': data['states'][0]['task']['charge'] = 1
    if bad == 'admission': data['states'][0]['task']['role'] = 'synthesis candidate'
    if bad == 'weights': data['experimental_occupancy_weights'] = True
    if bad == 'old_values': data['historical_values_transferred'] = True
    with pytest.raises(ValueError): validate_inventory(data)
