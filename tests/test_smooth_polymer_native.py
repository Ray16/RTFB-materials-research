import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('alternate_native_test',ROOT/'scripts/pipeline/smooth_polymer_native.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_native_columns_keep_protocol_and_largest_force_coordinates():
    gradient=np.zeros((5,3));gradient[2,0]=.1;gradient[3,1]=-.2
    columns=module.fd_columns(gradient,['C','O','H','H','Cl'])
    assert set(columns)=={(0,1),(1,2),(2,0),(3,1)}


def test_native_no_oxygen_uses_last_atom_and_rejects_nonfinite_forces():
    assert (3,2) in module.fd_columns(np.zeros((4,3)),['C','H','H','Cl'])
    with pytest.raises(ValueError):module.fd_columns(np.full((4,3),np.nan),['C','H','H','Cl'])


@pytest.fixture
def partition():
    return dict(experimental_method_id='smooth',primary_candidate_admission=False,
        derivative_checks_passed=True,symmetry_source='xTB autodetection on exact geometry',
        nuclear={'nuclear_by_cutoff_eV':{'25':.1,'50':.2,'100':.3}},spin_partition_eV=0.,
        G_by_cutoff_eV={'25':-999.9,'50':-999.8,'100':-999.7})


def test_alternate_G_closes_from_exact_energy_parent(partition):
    module.validate_partition(partition,-1000.,'smooth')
    with pytest.raises(ValueError):module.validate_partition(partition,-1001.,'smooth')


@pytest.mark.parametrize('change',['model','primary','derivative','symmetry','G','nan'])
def test_unqualified_or_mixed_model_native_records_are_rejected(partition,change):
    if change=='model':partition['experimental_method_id']='legacy'
    elif change=='primary':partition['primary_candidate_admission']=True
    elif change=='derivative':partition['derivative_checks_passed']=False
    elif change=='symmetry':partition['symmetry_source']='imposed comparison number one'
    elif change=='G':partition['G_by_cutoff_eV']['50']+=.1
    else:partition['nuclear']['nuclear_by_cutoff_eV']['50']=float('nan')
    with pytest.raises(ValueError):module.validate_partition(partition,-1000.,'smooth')
