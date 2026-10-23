import importlib.util
from pathlib import Path
import pytest

PATH=Path(__file__).resolve().parents[1]/'scripts/plotting/candidates/plot_descriptor_pareto.py'
spec=importlib.util.spec_from_file_location('descriptor_front_plot',PATH);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def row(key,first,second,capacity):
    return dict(id=key,E1_vs_Fc_V=first,E2_vs_Fc_V=second,capacity_nominal_mAh_g=capacity)


def test_joint_front_is_not_a_projected_2D_boundary():
    # b wins in the first-potential/capacity projection, but a has the better
    # second potential. Both belong to the three-objective descriptor front.
    a=row('a',-1,-1,100);b=row('b',-1.5,-.5,110);c=row('c',-.4,-.4,90)
    assert module.nondominated([a,b,c])=={'a','b'}


def test_duplicate_identity_cannot_be_silently_overwritten():
    with pytest.raises(ValueError,match='Duplicate'):
        module.nondominated([row('a',-1,-1,100),row('a',-2,-2,200)])
