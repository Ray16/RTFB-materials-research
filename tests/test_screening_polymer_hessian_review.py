"""Independent harmonic-potential and corrupted-displacement checks."""
import copy
import importlib.util

import numpy as np
import pytest

from redox.core.common import ROOT


@pytest.fixture
def panel():
    spec=importlib.util.spec_from_file_location('independent_polymer_columns',ROOT/'scripts/validation/review_screening_polymer_hessian.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    matrix=np.array([[2.,.1,0.],[.1,3.,.2],[0.,.2,4.]])
    coordinates=np.array([[.01,.02,.03]])
    center=dict(coordinates_Bohr=coordinates.tolist(),gradient_Ha_Bohr=(matrix@coordinates.ravel()).reshape(1,3).tolist())
    points={};checks=[]
    for axis in range(3):
        for step in (.0005,.001):
            pair=[]
            for sign in (-1,1):
                xyz=coordinates.copy();xyz[0,axis]+=sign*step
                key=f'{axis}_{step}_{sign}'
                points[key]=dict(coordinates_Bohr=xyz.tolist(),gradient_Ha_Bohr=(matrix@xyz.ravel()).reshape(1,3).tolist(),SCF_converged=True,spin=dict(spin_contam=0.))
                pair.append(key)
            checks.append(dict(atom=0,axis=axis,step_Bohr=step,relative_error=0.,minus=pair[0],plus=pair[1]))
    return module.column_checks,matrix,center,checks,points


def test_complete_harmonic_hessian_is_reconstructed_from_both_steps(panel):
    checker,matrix,center,checks,points=panel
    errors=checker(matrix,center,[(0,0),(0,1),(0,2)],checks,lambda path:points[str(path)])
    assert len(errors)==6 and max(errors)<1e-12


@pytest.mark.parametrize('change',['coordinate','gradient','spin','matrix','missing','duplicate'])
def test_wrong_physical_evidence_cannot_pass_column_review(panel,change):
    checker,matrix,center,checks,points=copy.deepcopy(panel)
    point=points[checks[0]['plus']]
    if change=='coordinate':point['coordinates_Bohr'][0][0]+=.01
    if change=='gradient':point['gradient_Ha_Bohr'][0][1]+=.01
    if change=='spin':point['spin']['spin_contam']=.2
    if change=='matrix':matrix[0,0]+=.5
    if change=='missing':checks.pop()
    if change=='duplicate':checks.append(checks[0])
    with pytest.raises(ValueError):checker(matrix,center,[(0,0),(0,1),(0,2)],checks,lambda path:points[str(path)])
