import importlib.util
import json
from pathlib import Path

import numpy as np
import pyscf
import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('cpu_ghost_review',ROOT/'scripts/validation/review_counterpoise_cpu_ghost.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_launch_arguments_must_be_unique():
    assert module.argument(['run','--own','/a'],'--own')==Path('/a')
    with pytest.raises(ValueError,match='Missing or duplicated'):
        module.argument(['run','--own','/a','--own','/b'],'--own')


def test_running_ghost_cannot_be_promoted(tmp_path,monkeypatch):
    monkeypatch.setattr(module,'ROOT',tmp_path)
    raw=tmp_path/'calcs'/'ghost';raw.mkdir(parents=True)
    (raw/'launch.json').write_text(json.dumps({'versions':{'pyscf':pyscf.__version__,'numpy':np.__version__}}))
    (raw/'status.json').write_text(json.dumps({'status':'CPU_ghost_pilot_SCF_running_unqualified'}))
    with pytest.raises(ValueError,match='Completed but unqualified'):
        module.review(raw,tmp_path/'review.json')
    assert not (tmp_path/'review.json').exists()
