"""The interaction/deformation split must retain the balanced association cycle."""
import importlib.util
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('fragment_energy',ROOT/'scripts/pipeline/polymer_fragment_interactions.py')
P=importlib.util.module_from_spec(spec);spec.loader.exec_module(P)


def test_frozen_interaction_excludes_the_cost_of_segment_deformation():
    result=P.decompose(pair=-20.,monomer=-9.5,a=-9.,b=-9.)
    assert result['frozen_interaction_deltaE_kJ_mol']==pytest.approx(-2*96.485)
    assert result['segment_deformation_deltaE_kJ_mol']==pytest.approx(96.485)
    assert result['relaxed_association_deltaE_kJ_mol']==pytest.approx(-96.485)
    assert abs(result['cycle_closure_error_kJ_mol'])<1e-8


def test_negative_composite_deformation_is_reported_instead_of_clipped():
    result=P.decompose(pair=-21.,monomer=-9.,a=-9.5,b=-9.5)
    assert result['segment_deformation_deltaE_kJ_mol']<0
    assert result['frozen_interaction_deltaE_kJ_mol']+result['segment_deformation_deltaE_kJ_mol']==pytest.approx(result['relaxed_association_deltaE_kJ_mol'])
