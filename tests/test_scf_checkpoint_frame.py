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


def test_actual_elements_cannot_be_replaced_by_metadata(tmp_path):
    expected=molecule();actual=molecule(atom='He 0 0 0; He 0 0 1.4');path=tmp_path/'SCF.h5';chkfile.save_mol(actual,str(path))
    with pytest.raises(ValueError,match='checkpoint'):
        checkpoint_frame(path,charge=0,mult=1,basis='sto-3g',coordinates_Bohr=expected.atom_coords(),symbols=['H','H'])


@pytest.mark.parametrize('value',[float('nan'),float('inf'),None,True,.10001,-.10001])
def test_missing_nonfinite_or_outside_gate_spin_is_rejected(value):
    from redox.validation.scf_frame import check_spin_contamination
    with pytest.raises(ValueError):check_spin_contamination(value)


def test_frozen_checker_is_content_bound(tmp_path):
    from redox.validation.scf_frame import freeze_reviewer
    frozen=freeze_reviewer(tmp_path);assert freeze_reviewer(tmp_path)==frozen
    frozen.write_text('changed')
    with pytest.raises(ValueError):freeze_reviewer(tmp_path)


def test_stricter_native_coordinate_gate_is_preserved(tmp_path):
    mol=molecule();path=tmp_path/'SCF.h5';chkfile.save_mol(mol,str(path));coordinates=mol.atom_coords();coordinates[0,0]+=3e-11
    args=dict(charge=0,mult=1,basis='sto-3g',coordinates_Bohr=coordinates,symbols=['H','H'])
    checkpoint_frame(path,**args)
    with pytest.raises(ValueError):checkpoint_frame(path,**args,maximum_error_Bohr=1e-12)


def spin_checkpoint(tmp_path,mol,coefficients,occupations):
    import h5py
    path=tmp_path/'spin.h5';chkfile.save_mol(mol,str(path))
    with h5py.File(path,'a') as f:
        f.create_dataset('scf/mo_coeff',data=coefficients);f.create_dataset('scf/mo_occ',data=occupations)
    return path


def orthonormal_orbitals(mol):
    import numpy as np
    values,vectors=np.linalg.eigh(mol.intor_symmetric('int1e_ovlp'))
    return (vectors/np.sqrt(values))@vectors.T


def test_restricted_saved_wavefunction_has_checked_singlet_spin(tmp_path):
    import numpy as np
    from redox.validation.scf_frame import checkpoint_spin
    mol=molecule();path=spin_checkpoint(tmp_path,mol,orthonormal_orbitals(mol),np.array([2.,0.]))
    r=checkpoint_spin(path);assert abs(r['s_squared'])<1e-12 and r['wavefunction_spin_reconstructed']


def test_unrestricted_one_electron_spin_is_reconstructed(tmp_path):
    import numpy as np
    from redox.validation.scf_frame import checkpoint_spin
    mol=molecule(atom='H 0 0 0',spin=1);c=orthonormal_orbitals(mol)
    path=spin_checkpoint(tmp_path,mol,np.array([c,c]),np.array([[1.],[0.]]))
    assert checkpoint_spin(path)['s_squared']==pytest.approx(.75)


def test_declared_singlet_cannot_hide_a_spin_contaminated_wavefunction(tmp_path):
    import numpy as np
    from redox.validation.scf_frame import checkpoint_spin
    mol=molecule();c=orthonormal_orbitals(mol)
    # Occupied alpha/beta orbitals are orthogonal: <S²>=1 despite spin=0 metadata.
    path=spin_checkpoint(tmp_path,mol,np.array([c,c[:,::-1]]),np.array([[1.,0.],[1.,0.]]))
    with pytest.raises(ValueError,match='spin contamination'):checkpoint_spin(path)


@pytest.mark.parametrize('bad',['normalization','electron_count','nonfinite'])
def test_wrong_saved_orbitals_or_occupations_cannot_pass(tmp_path,bad):
    import numpy as np
    from redox.validation.scf_frame import checkpoint_spin
    mol=molecule();c=orthonormal_orbitals(mol);occ=np.array([2.,0.])
    if bad=='normalization':c*=1.01
    if bad=='electron_count':occ[1]=2.
    if bad=='nonfinite':c[0,0]=float('nan')
    path=spin_checkpoint(tmp_path,mol,c,occ)
    with pytest.raises(ValueError):checkpoint_spin(path)
