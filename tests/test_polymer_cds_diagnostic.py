import importlib.util
from pathlib import Path

import numpy as np
import pytest

spec=importlib.util.spec_from_file_location('tested_CDS_split',Path(__file__).parents[1]/'scripts/validation/diagnose_polymer_components.py')
diagnostic=importlib.util.module_from_spec(spec);spec.loader.exec_module(diagnostic)


def test_decomposition_distinguishes_surface_error_from_other_force_error():
    h=.001;column=np.array([2.,3.,4.]);cds=np.array([.1,.2,.3]);jump=np.array([1e-5,0,0])
    plus=h*column+jump;minus=-h*column
    values,*_=diagnostic.split_residual(plus,minus,column,h*cds+jump,-h*cds,cds,h)
    assert values['remaining_gradient_difference_residual_norm']<1e-16
    values,*_=diagnostic.split_residual(plus,minus,column,h*cds,-h*cds,cds,h)
    assert values['remaining_error_over_full_residual']==pytest.approx(1.)
