import copy
import importlib.util
import json
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('aq_slide_build',ROOT/'scripts/presentations/build_colleague_workflow.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
SOURCE=ROOT/'results/discovery/anthraquinone_carboxylate_matrix_2026-10-26'


def test_aq_slide_evidence_is_source_bound():
    paths=module.validate_aq_panel()
    assert len(paths)==4


@pytest.mark.parametrize('mutation',['route','decision','property'])
def test_aq_slide_cannot_promote_foreign_or_unqualified_evidence(mutation):
    data=[json.loads((SOURCE/name).read_text()) for name in ('summary.json','matrix.json','route_screen.json','all_axis_decisions.json')]
    summary,matrix,gate,decisions=copy.deepcopy(data)
    if mutation=='route':gate['candidates'][0]['evidence']['exact_design_key']='different'
    elif mutation=='decision':decisions['decisions'][0]['status']='eligible'
    else:matrix['candidates'][0]['properties_computed']=True
    with pytest.raises(ValueError,match='exact design'):
        module.check_aq_panel(summary,matrix,gate,decisions)
