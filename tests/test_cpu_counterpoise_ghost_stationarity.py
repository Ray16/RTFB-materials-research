import importlib.util
import json
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('ghost_stationarity',ROOT/'scripts/validation/check_cpu_counterpoise_ghost_stationarity.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_requires_completed_independent_ghost_review(tmp_path,monkeypatch):
    monkeypatch.setattr(module,'ROOT',tmp_path)
    raw=tmp_path/'calcs'/'ghost';raw.mkdir(parents=True);scratch=tmp_path/'calcs'/'scratch';scratch.mkdir()
    monkeypatch.setenv('PYSCF_TMPDIR',str(scratch))
    review=tmp_path/'review.json';review.write_text(json.dumps({'status':'waiting','physical_gradient_recomputed':False,'full_counterpoise_computed':False}))
    out=tmp_path/'calcs'/'stationarity'
    with pytest.raises(ValueError,match='independent ghost checkpoint review'):
        module.check(raw,review,out)
    assert not out.exists()
