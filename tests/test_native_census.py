import importlib.util
import json
from pathlib import Path

import pytest

spec=importlib.util.spec_from_file_location('tested_native_census',Path(__file__).parents[1]/'scripts/validation/collect_native_census.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


@pytest.fixture
def replay(tmp_path):
    record=dict(key='all_same__occ_none__conf001',configuration='all_same',occupied_sites=[],status='relocated_native_physical_replay_passed',
        G_by_cutoff_eV={'25':-1.,'50':-1.,'100':-1.},maximum_column_error=.001,minimum_frequency_cm=5.,graph_sha256='graph',method_id='method')
    data=dict(status='namespace_native_physical_records_replayed',states=1,records=[record],producer_sha256={},quantum_calculations_run=False,global_ensemble_convergence_proven=False)
    path=tmp_path/'replay.json';path.write_text(json.dumps(data));return path


def test_missing_states_stay_missing_and_counts_do_not_imply_an_ensemble(replay):
    result=module.census([replay],['all_same'],[(),(2,)],20)
    assert result['checked_states']==1 and result['expected_states']==40
    assert result['groups'][1]['checked']==0 and len(result['groups'][0]['missing_keys'])==19
    assert not result['global_ensemble_convergence_proven']


@pytest.mark.parametrize('change',['duplicate','unknown_graph','bad_G','bad_FD','bad_minimum','promoted','outside_index'])
def test_changed_coverage_or_unqualified_records_are_rejected(replay,change):
    data=json.loads(replay.read_text());row=data['records'][0];paths=[replay]
    if change=='duplicate':paths=[replay,replay]
    if change=='unknown_graph':row['occupied_sites']=[3]
    if change=='bad_G':del row['G_by_cutoff_eV']['50']
    if change=='bad_FD':row['maximum_column_error']=.011
    if change=='bad_minimum':row['minimum_frequency_cm']=-6.
    if change=='promoted':data['global_ensemble_convergence_proven']=True
    if change=='outside_index':row['key']='all_same__occ_none__conf021'
    replay.write_text(json.dumps(data))
    with pytest.raises(ValueError):module.census(paths,['all_same'],[(),(2,)],20)
