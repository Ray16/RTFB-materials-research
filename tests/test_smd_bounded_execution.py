import importlib.util
from pathlib import Path
import numpy as np
import pytest

SPEC=importlib.util.spec_from_file_location("bounded_SMD",Path(__file__).resolve().parents[1]/"scripts/pipeline/smd_bounded_execution.py")
MODULE=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


@pytest.mark.parametrize("transpose",[False,True])
@pytest.mark.parametrize("columns",[None,1,3])
def test_contracted_response_matches_dense_equation(transpose,columns):
    rng=np.random.default_rng(8);D=rng.normal(size=(7,7));A=rng.uniform(.1,2,7)
    vector=rng.normal(size=(7,) if columns is None else (7,columns))
    matrix=(D*A).T if transpose else D*A
    expected=-.9*(vector-matrix@vector/(2*np.pi))
    np.testing.assert_allclose(MODULE.contracted_R(np,D,A,vector,.9,transpose),expected,atol=1e-14,rtol=1e-14)


def test_workspace_budget_accounts_for_existing_allocations_and_reusable_pool():
    assert MODULE.budget_available(100,60,5,80)==45
    assert MODULE.budget_available(100,10,0,80)==0
