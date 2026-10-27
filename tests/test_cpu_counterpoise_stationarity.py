import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('cpu_stationarity',ROOT/'scripts/validation/check_cpu_counterpoise_stationarity.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
BASE=ROOT/'calcs/validation/2026-10-26_counterpoise/CPU_fixed_auxiliary_factory'


def test_rejects_checkpoint_for_a_different_geometry():
    record=json.loads((BASE/'summary.json').read_text());coords=np.array(record['source_coordinates_A']);coords[0,0]+=0.01
    with pytest.raises(ValueError,match='Checkpoint molecule'):
        module.evaluate(record['source_symbols'],coords,record['A_indices'],BASE/'grid_coordinates_AU.npy',
                        BASE/'grid_weights.npy',BASE/'A_own/SCF.h5',0.)


def test_rejects_bad_grid_before_calculation(tmp_path):
    record=json.loads((BASE/'summary.json').read_text());coords=tmp_path/'coords.npy';weights=tmp_path/'weights.npy'
    np.save(coords,np.array([[0.,0.,float('nan')]]));np.save(weights,np.array([1.]))
    with pytest.raises(ValueError,match='Invalid common-grid'):
        module.evaluate(record['source_symbols'],record['source_coordinates_A'],record['A_indices'],
                        coords,weights,BASE/'A_own/SCF.h5',0.)


def test_pilot_requires_completed_independent_review(tmp_path,monkeypatch):
    monkeypatch.setattr(module,'ROOT',tmp_path)
    raw=tmp_path/'calcs'/'pilot';raw.mkdir(parents=True);scratch=tmp_path/'calcs'/'scratch';scratch.mkdir()
    monkeypatch.setenv('PYSCF_TMPDIR',str(scratch))
    review=tmp_path/'review.json';review.write_text(json.dumps(dict(status='waiting')))
    with pytest.raises(ValueError,match='independent pilot checkpoint review'):
        module.check_pilot(raw,review,tmp_path/'calcs'/'stationarity')
