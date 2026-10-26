import importlib.util
import json
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('tempo_slide_builder_test',ROOT/'scripts/presentations/build_colleague_workflow.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_slide_gate_accepts_checked_evidence_and_rejects_uncomputed_grafting_claim(tmp_path):
    source=ROOT/'results/grafting/tempo_exact_support_2026-10-26'
    tempo,figure=module.validate_tempo(source/'summary.json')
    assert tempo['all_axis_status']=='pending' and figure.is_file()
    bad=json.loads((source/'summary.json').read_text());bad['actual_route_deltaG_computed']=True
    file=tmp_path/'summary.json';file.write_text(json.dumps(bad))
    with pytest.raises(ValueError,match='Unqualified TEMPO'):
        module.validate_tempo(file)
