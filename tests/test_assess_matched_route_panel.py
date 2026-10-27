import importlib.util
import json
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('matched_route_panel',ROOT/'scripts/pipeline/assess_matched_route_panel.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


@pytest.mark.parametrize('directory',['pmdi_linker_matrix_2026-10-26','phenothiazine_linker_matrix_2026-10-26'])
def test_shared_panel_gate_stays_pending_and_rejects_foreign_identity(directory,tmp_path):
    source=ROOT/'results/discovery'/directory;policy=ROOT/'config/candidate_feasibility.json'
    result=module.assess(source/'matrix.json',source/'route_screen.json',policy)
    assert len(result['decisions'])==4 and all(d['status']=='pending' for d in result['decisions'])
    gate=json.loads((source/'route_screen.json').read_text());gate['candidates'][0]['evidence']['exact_design_key']='wrong'
    altered=tmp_path/'gate.json';altered.write_text(json.dumps(gate))
    with pytest.raises(ValueError,match='exact candidate'):
        module.assess(source/'matrix.json',altered,policy)
