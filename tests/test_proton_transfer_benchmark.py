"""External relative basicity must correspond to an atom/charge-balanced transfer."""
import importlib.util
from pathlib import Path

import pytest


@pytest.fixture
def benchmark():
    path=Path(__file__).resolve().parents[1]/'scripts/validation/proton_transfer_benchmark.py'
    spec=importlib.util.spec_from_file_location('proton_transfer_test',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def test_declared_pyridine_transfer_is_balanced_without_HCl_or_chloride(benchmark):
    assert benchmark.balanced_transfer('c1ccncc1','c1cc[nH+]cc1',
        'C1CCC2=NCCCN2CC1','C1CCC2=[NH+]CCCN2CC1')


@pytest.mark.parametrize('product',['c1cc[nH+]cc1.[Cl-]','C[NH+](C)C'])
def test_an_extra_counterion_or_different_base_cannot_be_silently_substituted(benchmark,product):
    assert not benchmark.balanced_transfer('c1ccncc1',product,
        'C1CCC2=NCCCN2CC1','C1CCC2=[NH+]CCCN2CC1')
