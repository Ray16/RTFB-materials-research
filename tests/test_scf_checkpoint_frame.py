import pytest
pytest.importorskip('pyscf')
from pyscf import gto
from pyscf.lib import chkfile
from redox.validation.scf_frame import checkpoint_frame


def molecule(**changes):
    args=dict(atom='H 0 0 0; H 0 0 1.4',unit='Bohr',basis='sto-3g',charge=0,spin=0,verbose=0)
    args.update(changes);return gto.M(**args)


def test_actual_saved_molecule_is_checked(tmp_path):
    mol=molecule();path=tmp_path/'SCF.h5';chkfile.save_mol(mol,str(path))
    r=checkpoint_frame(path,charge=0,mult=1,basis='sto-3g',coordinates_Bohr=mol.atom_coords(),symbols=['H','H'])
    assert r['maximum_coordinate_error_Bohr']==0 and r['multiplicity']==1


@pytest.mark.parametrize('changes',[{'atom':'H 0 0 0; H 0 0 1.5'},{'charge':1,'spin':1},{'basis':'6-31g'}])
def test_metadata_cannot_hide_changed_actual_checkpoint(tmp_path,changes):
    expected=molecule();actual=molecule(**changes);path=tmp_path/'SCF.h5';chkfile.save_mol(actual,str(path))
    with pytest.raises(ValueError,match='checkpoint'):
        checkpoint_frame(path,charge=0,mult=1,basis='sto-3g',coordinates_Bohr=expected.atom_coords(),symbols=['H','H'])
