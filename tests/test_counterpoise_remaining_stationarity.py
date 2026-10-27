import importlib.util
import json
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('remaining_stationarity',ROOT/'scripts/validation/check_counterpoise_remaining_stationarity.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_requires_independent_completed_term_review(tmp_path,monkeypatch):
    monkeypatch.setattr(module,'ROOT',tmp_path)
    raw=tmp_path/'calcs'/'term';raw.mkdir(parents=True);scratch=tmp_path/'calcs'/'scratch';scratch.mkdir()
    monkeypatch.setenv('PYSCF_TMPDIR',str(scratch))
    review=tmp_path/'review.json';review.write_text(json.dumps({'status':'waiting','physical_gradient_recomputed':False,'full_counterpoise_computed':False}))
    out=tmp_path/'calcs'/'stationarity'
    with pytest.raises(ValueError,match='independent term review'):
        module.check(raw,review,out)
    assert not out.exists()
