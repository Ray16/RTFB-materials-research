"""A failed force check cannot be relabeled as resource-only recovery."""
import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('returned_geometry',ROOT/'scripts/pipeline/finish_polymer_geometry.py')
MODULE=importlib.util.module_from_spec(SPEC);SPEC.loader.exec_module(MODULE)


def records():
    return (dict(status='geometry_stage_failed_preserved',error='cudaErrorMemoryAllocation: out of memory'),
            dict(optimized=True,geometry_optimization_converged=True,optimization_SCF_converged=True,
                 geometry_only_result=True,charge=0,mult=1,solvent='acetonitrile'))


def test_converged_resource_interruption_keeps_the_original_optimizer_flags():
    status,result=records();original=dict(result)
    MODULE.validate_predecessor(status,result)
    assert result==original


@pytest.mark.parametrize('field',['optimized','geometry_optimization_converged','optimization_SCF_converged','geometry_only_result'])
def test_missing_optimizer_or_scf_acceptance_is_not_filled(field):
    status,result=records();result[field]=False
    with pytest.raises(ValueError):MODULE.validate_predecessor(status,result)


def test_physical_failure_or_different_medium_cannot_be_resource_recovery():
    status,result=records();status['error']='Final force fails finite differences'
    with pytest.raises(ValueError):MODULE.validate_predecessor(status,result)
    status,result=records();result['solvent']='water'
    with pytest.raises(ValueError):MODULE.validate_predecessor(status,result)
