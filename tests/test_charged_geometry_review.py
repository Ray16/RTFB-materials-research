"""Corrupt final-state binding and graph changes must fail independently."""
import copy
import importlib.util
from pathlib import Path

import numpy as np
import pytest
from pyscf.lib import param
from rdkit import Chem
from rdkit.Chem import AllChem


@pytest.fixture
def checker():
    path=Path(__file__).resolve().parents[1]/'scripts/validation/review_charged_polymer_optimization.py'
    spec=importlib.util.spec_from_file_location('test_ionic_geometry_review',path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result


@pytest.fixture
def physical():
    xyz=np.array([[0.,0.,0.],[1.,0.,0.]])
    plan=dict(atoms=2,charge=2,mult=1,optimizer_settings=dict(convergence_gmax=4.5e-4,convergence_grms=3e-4))
    returned=dict(converged=True,coordinates_A=xyz.tolist(),charge=2,mult=1)
    final=dict(coordinates_Bohr=(xyz/param.BOHR).tolist(),gradient_Ha_Bohr=[[1e-5,0.,0.],[0.,1e-5,0.]],
        optimizer_converged=True,SCF_converged=True,charge=2,mult=1,spin=dict(spin_contam=0.),
        native_G_qualified=False,E_Ha=-10.,gradient_max_Ha_Bohr=1e-5)
    return plan,returned,final,xyz.copy(),xyz.copy(),-10.


def test_fresh_final_state_binding_passes(checker,physical):
    result=checker.physical_checks(*physical)
    assert result['maximum_atom_gradient_Ha_Bohr']==1e-5


@pytest.mark.parametrize('change',['checkpoint_energy','nonfinite_checkpoint_energy','checkpoint_geometry','returned_geometry','charge','spin','force_norm','unconverged'])
def test_unbound_or_nonstationary_result_fails(checker,physical,change):
    values=list(copy.deepcopy(physical));plan,returned,final,xyz,chk,energy=values
    if change=='checkpoint_energy':values[-1]+=.1
    if change=='nonfinite_checkpoint_energy':values[-1]=float('nan')
    if change=='checkpoint_geometry':chk[0,0]+=.01
    if change=='returned_geometry':returned['coordinates_A'][0][0]+=.01
    if change=='charge':final['charge']=4
    if change=='spin':final['spin']['spin_contam']=.11
    if change=='force_norm':
        final['gradient_Ha_Bohr'][0]=[.0004]*3;final['gradient_max_Ha_Bohr']=.0004
    if change=='unconverged':returned['converged']=False
    with pytest.raises(ValueError):checker.physical_checks(*values)


def test_connectivity_and_stereochemistry_are_both_required(checker,tmp_path):
    smiles='C[C@H](F)Cl';mol=Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol,randomSeed=42)==0
    AllChem.MMFFOptimizeMolecule(mol)
    initial=tmp_path/'initial.xyz';final=tmp_path/'final.xyz'
    initial.write_text(Chem.MolToXYZBlock(mol));final.write_text(Chem.MolToXYZBlock(mol))
    assert checker.graph_checks(smiles,initial,final)['stereochemistry_retained']
    reflected=Chem.Mol(mol);conf=reflected.GetConformer()
    for i in range(reflected.GetNumAtoms()):
        point=conf.GetAtomPosition(i);conf.SetAtomPosition(i,(-point.x,point.y,point.z))
    final.write_text(Chem.MolToXYZBlock(reflected))
    with pytest.raises(ValueError,match='stereochemistry'):checker.graph_checks(smiles,initial,final)
    broken=Chem.Mol(mol);conf=broken.GetConformer();point=conf.GetAtomPosition(broken.GetNumAtoms()-1)
    conf.SetAtomPosition(broken.GetNumAtoms()-1,(point.x+5,point.y,point.z))
    final.write_text(Chem.MolToXYZBlock(broken))
    with pytest.raises(ValueError,match='graph'):checker.graph_checks(smiles,initial,final)
