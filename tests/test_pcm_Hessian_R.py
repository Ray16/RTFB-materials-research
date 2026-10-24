"""Tensor/batch orientation agrees with dense R for both transpose directions."""
import ast
import importlib.util
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("PCM_Hessian_contracted_R",
    ROOT / "scripts/pipeline/smd_first_response_workspace_v2.py")
m = importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)


@pytest.mark.parametrize("transpose", [False, True])
def test_response_batch_orientation_matches_dense_asymmetric_R(transpose):
    rng = np.random.default_rng(20261024)
    D, A = rng.normal(size=(19,19)), rng.uniform(.1,1.,size=19)
    factor = .97
    R = -factor*(np.eye(19)-D*A/(2*math.pi))
    values = rng.normal(size=(5,3,19))
    def apply(vector, R_transpose=False):
        return -factor*(vector-(A[:,None]*(D.T@vector) if R_transpose else D@(A[:,None]*vector))/(2*math.pi))
    actual = m.apply_R_batch(SimpleNamespace(left_multiply_R=apply), values, transpose)
    expected = np.einsum('ij,Adj->Adi', R.T if transpose else R, values)
    np.testing.assert_allclose(actual,expected,atol=1e-12,rtol=1e-12)


def test_installed_templates_transform_only_supported_dense_R_patterns():
    import importlib.util
    package = Path(importlib.util.find_spec('gpu4pyscf').origin).parent
    source = package/'solvent/hessian/pcm.py'
    tree = ast.parse(source.read_text())
    functions = {n.name:ast.get_source_segment(source.read_text(),n)+'\n'
                 for n in tree.body if isinstance(n,ast.FunctionDef)}
    for name in ('get_dqsym_dx_fix_vgrids','get_dqsym_dx_fix_K_R','analytical_hess_solver'):
        transformed = m.transform(name,functions[name])
        assert '_apply_R_batch(' in transformed
        compile(transformed,'transformed_PCM_test','exec')
        with pytest.raises(ValueError):
            m.transform(name,functions[name].replace('cupy.eye','changed_eye'))


def test_unexpected_response_layout_is_not_silently_reshaped():
    with pytest.raises(ValueError,match='tensor'):
        m.apply_R_batch(None,np.zeros((3,19)),False)
