from pathlib import Path

import pytest

from redox.workflows.archived_energy_namespace import validate_plan


def plan():
    return dict(tasks=[dict(charge=0, mult=1)], experimental_model=dict(solvent='acetonitrile'),
                execution_driver='/old/scripts/pipeline/smooth_polymer_energy_bounded.py',
                native_G_computed=False, primary_candidate_admission=False,
                energy_only_execution=True, maximum_peak_GPU_MiB=31744, maximum_host_RSS_GiB=128.)


def test_predeclared_bounded_energy_class_retains_scope():
    validate_plan(plan(), Path('/old'))


@pytest.mark.parametrize('field,value', [('native_G_computed', True), ('primary_candidate_admission', True),
                                      ('maximum_host_RSS_GiB', float('nan')),
                                      ('maximum_host_RSS_GiB', 256.), ('energy_only_execution', False),
                                      ('execution_driver', '/old/scripts/other.py')])
def test_rejects_scope_or_unqualified_resource_changes(field, value):
    row = plan()
    row[field] = value
    with pytest.raises(ValueError):
        validate_plan(row, Path('/old'))


def test_rejects_a_different_charge_state_or_solvent():
    row = plan()
    row['tasks'][0]['charge'] = -1
    with pytest.raises(ValueError):
        validate_plan(row, Path('/old'))
    row = plan()
    row['experimental_model']['solvent'] = 'water'
    with pytest.raises(ValueError):
        validate_plan(row, Path('/old'))


def test_absent_raw_batch_uses_private_mountpoints_and_read_only_siblings(tmp_path):
    from redox.workflows.archived_energy_namespace import energy_mount
    archive = tmp_path/'archive'
    sibling = archive/'calcs/polymer_models/preserved'
    sibling.mkdir(parents=True)
    value = dict(prior_root='/old', archive=str(archive),
                 original_campaign='/old/results/polymer_models/campaigns/new/key',
                 original_raw_parent='/old/calcs/polymer_models/new',
                 attempt=str(tmp_path/'new/campaign'), raw_parent=str(tmp_path/'new/raw'))
    argv = energy_mount(value)
    parent_index = argv.index('/old/calcs/polymer_models')
    assert argv[parent_index-1] == '--tmpfs'
    sibling_index = argv.index(str(sibling))
    assert argv[sibling_index-1] == '--ro-bind'
    assert argv.index('/old/calcs/polymer_models/new') > parent_index
