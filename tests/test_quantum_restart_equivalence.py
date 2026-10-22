"""Restart equivalence needs matching physical outputs, not just a process exit."""
from copy import deepcopy
import importlib.util
from pathlib import Path

from ase import Atoms
import numpy as np
import pytest


@pytest.fixture
def control():
    path=Path(__file__).resolve().parents[1]/'scripts/validation/quantum_restart_equivalence.py'
    spec=importlib.util.spec_from_file_location('restart_control_test',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def record():
    return [Atoms('COHH',positions=[[0,0,0],[1.2,0,0],[-.6,.9,0],[-.6,-.9,0]]),
            dict(e_smd_eV=-10.,e_gas_eV=-9.),dict(frequencies_cm=[100.,200.,300.,400.,500.,600.]),dict(G_eV=-8.)]


def test_rigid_rotation_and_translation_do_not_fake_a_restart_difference(control):
    a=record();b=deepcopy(a);b[0].positions=b[0].positions@np.array([[0,1,0],[-1,0,0],[0,0,1]])+[4,5,6]
    assert control.restart_differences(a,b)['aligned_RMSD_Angstrom']<1e-12


@pytest.mark.parametrize('field',['energy','free_energy','vibration','chemistry'])
def test_wrong_energy_vibration_or_chemistry_cannot_pass(control,field):
    a=record();b=deepcopy(a)
    if field=='energy':b[1]['e_smd_eV']+=1e-3
    if field=='free_energy':b[3]['G_eV']+=.001
    if field=='vibration':b[2]['frequencies_cm'][0]+=1
    if field=='chemistry':b[0].set_chemical_symbols('CNHH')
    with pytest.raises(ValueError):control.restart_differences(a,b)
