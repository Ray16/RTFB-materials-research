import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('checkpoint_frames_test',ROOT/'scripts/pipeline/stage_interrupted_conformers.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_bounded_selection_retains_lowest_frames_and_source_indices():
    frames=['4.0','1.0','-2.0','3.0','-1.0']
    rows,count=module.lowest_frames(None,lambda _:iter(frames),float,maximum=2)
    assert count==5
    assert [(e,i) for e,i,_ in rows]==[(-2.0,2),(-1.0,4)]


def test_nonfinite_energy_cannot_seed_a_quantum_job():
    with pytest.raises(ValueError,match='nonfinite'):
        module.lowest_frames(None,lambda _:iter(['nan']),float)


def test_truncated_checkpoint_is_never_silently_treated_as_complete():
    def frames(_):
        yield '1.0'
        raise ValueError('Incomplete XYZ frame')
    with pytest.raises(ValueError,match='Incomplete'):
        module.lowest_frames(None,frames,float)
