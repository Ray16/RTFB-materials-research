import copy
import json
from pathlib import Path

import pytest
from rdkit import Chem

from redox.workflows.reaction_component_namespace import validate_component
from redox.workflows.molecular_native_namespace import command
from redox.validation.molecular_native import cross_device_uma_review


@pytest.fixture
def reagent():
    recipe=json.loads((Path(__file__).parents[1]/'config/viologen_grafting_component.json').read_text())
    route={'route_binding':{'precursor_smiles':'C[n+]1ccc(-c2ccncc2)cc1'}}
    return recipe,route


def test_reagent_has_its_own_graph_charge_and_role(reagent):
    recipe,route=reagent;mol=validate_component(recipe,route)
    assert Chem.GetFormalCharge(mol)==1
    assert Chem.AddHs(mol).GetNumAtoms()==24
    assert recipe['role']=='reaction_component'
    assert Chem.MolToInchiKey(mol)!='PBNXFHGCSKYNHC-UHFFFAOYSA-N'


@pytest.mark.parametrize('change',[{'role':'candidate'},{'smiles':'c1ccncc1'},{'charge':2},
    {'mult':2},{'source_id':'../candidate'},{'UMA_fmax_eV_A':.1},{'top_k':0}])
def test_invalid_reagent_or_uncontrolled_warm_start_is_rejected(reagent,change):
    recipe,route=copy.deepcopy(reagent);recipe.update(change)
    with pytest.raises(ValueError):validate_component(recipe,route)


def test_spectator_convention_is_atom_and_charge_balanced(reagent):
    recipe,_=reagent
    reactants=[recipe['smiles'],'Cc1ccc(CCl)cc1']
    products=['Cc1ccc(C[n+]2ccc(-c3cc[n+](C)cc3)cc2)cc1','[Cl-]']
    def inventory(smiles):
        from collections import Counter
        atoms=Counter();charge=0
        for text in smiles:
            mol=Chem.AddHs(Chem.MolFromSmiles(text));charge+=Chem.GetFormalCharge(mol)
            atoms.update(a.GetSymbol() for a in mol.GetAtoms())
        return atoms,charge
    assert inventory(reactants)==inventory(products)
    assert len(products)-len(reactants)==0  # No net 1 M conversion term.


def test_private_reagent_manifest_is_read_only_in_the_namespace(tmp_path):
    case=tmp_path/'case';case.mkdir();(case/'library').mkdir();sources=tmp_path/'sources';sources.mkdir()
    spec=dict(case=str(case),sources=str(sources),prior_root='/original',archive='/archive',library_overlay=str(case/'library'))
    args=command(spec)
    assert any(args[i:i+3]==['--ro-bind',str(case/'library'),'/original/library'] for i in range(len(args)-2))
    spec['library_overlay']=str(tmp_path)
    with pytest.raises(ValueError):command(spec)


@pytest.fixture
def cross_device():
    reference=dict(numbers=[1],energy_eV=-10.,forces_eV_A=[[.1,0.,0.]])
    check=dict(energy_tolerance_eV=1e-4,force_tolerance_eV_A=1e-4,
        repeats=[dict(energy_eV=-10.,forces_eV_A=[[.1,0.,0.]],energy_error_eV=0.,force_max_error_eV_A=0.) for _ in range(3)])
    return check,reference


def test_cross_device_control_reconstructs_fixed_geometry_arrays(cross_device):
    assert len(cross_device_uma_review(*cross_device))==3


@pytest.mark.parametrize('change',['energy','force','nan','shape','missing','relaxed_gate'])
def test_control_cannot_pass_with_forged_small_reported_errors(cross_device,change):
    check,reference=cross_device
    if change=='energy':check['repeats'][0]['energy_eV']+=.01
    if change=='force':check['repeats'][0]['forces_eV_A'][0][0]+=.01
    if change=='nan':check['repeats'][0]['forces_eV_A'][0][0]=float('nan')
    if change=='shape':check['repeats'][0]['forces_eV_A']=[]
    if change=='missing':check['repeats'].pop()
    if change=='relaxed_gate':check['force_tolerance_eV_A']=.1
    with pytest.raises(ValueError):cross_device_uma_review(check,reference)
