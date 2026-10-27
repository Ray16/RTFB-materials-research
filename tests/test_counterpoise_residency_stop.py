import copy
import importlib.util
import json
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('cp_residency_review',ROOT/'scripts/validation/review_counterpoise_residency_stop.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
CASE=ROOT/'results/validation/local_counterpoise_2026-10-26_motif_v2/occ_2_3__pose0__cap109'


def records():
    paths=('coordinator.json','launch.json','host_guard/status.json','residency_guard/state.json',
           'residency_guard/interruption.json','tensor_runtime/status.json')
    return [json.loads((CASE/p).read_text()) for p in paths]


def test_recorded_pre_science_residency_stop_is_classified():
    result=module.classify(*records(),(CASE/'worker.log').read_text(),False,False)
    assert result['scientific_starts']==0 and not result['foreign_processes_signaled']


@pytest.mark.parametrize('mutation',['raw','stage','own_residency','no_foreign','changed_pid','missing_import'])
def test_started_or_ambiguous_attempts_cannot_be_reopened(mutation):
    rows=copy.deepcopy(records());log=(CASE/'worker.log').read_text();raw=stage=False
    if mutation=='raw':raw=True
    elif mutation=='stage':stage=True
    elif mutation=='own_residency':rows[3]['own_residency_observed']=True
    elif mutation=='no_foreign':rows[4]['foreign_resident_PIDs']=[]
    elif mutation=='changed_pid':rows[4]['pid']+=1
    else:log=''
    with pytest.raises(ValueError):module.classify(*rows,log,raw,stage)
