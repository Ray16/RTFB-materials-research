import numpy as np
import pytest

from redox.validation.solution_pair_review import density_spin


def test_closed_and_open_shell_densities_reconstruct_known_spin():
    overlap=np.eye(2)
    closed=density_spin(np.diag([2.,0.]),overlap,(2,1))
    assert closed['s_squared']==pytest.approx(0.)
    doublet=density_spin(np.array([np.diag([1.,1.]),np.diag([1.,0.])]),overlap,(3,2))
    assert doublet['s_squared']==pytest.approx(.75)


@pytest.mark.parametrize('mutation',['charge','multiplicity','nonfinite','shape','contamination','nonhermitian','fractional'])
def test_wrong_populations_shapes_or_spin_are_refused(mutation):
    overlap=np.eye(2);dm=np.array([np.diag([1.,1.]),np.diag([1.,0.])]);identity=(3,2)
    if mutation=='charge':identity=(2,2)
    if mutation=='multiplicity':identity=(3,1)
    if mutation=='nonfinite':dm[0,0,0]=float('nan')
    if mutation=='shape':dm=dm[:,:1,:1]
    if mutation=='nonhermitian':dm[0,0,1]=.1
    if mutation=='fractional':dm[1]=np.eye(2)*.5
    if mutation=='contamination':
        dm=np.array([np.diag([1.,0.]),np.diag([0.,1.])]);identity=(2,1)
    with pytest.raises(ValueError):density_spin(dm,overlap,identity)
