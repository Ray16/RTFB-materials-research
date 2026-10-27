import importlib.util
import json
from pathlib import Path

import pytest
from rdkit import Chem


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('ptz_source',ROOT/'scripts/pipeline/prepare_ptz_hydroxyethyl_candidate.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
IDENTITY=dict(CID=352306,InChIKey='FMENRHNHBUWRMN-UHFFFAOYSA-N',
              ConnectivitySMILES='C1=CC=C2C(=C1)N(C3=CC=CC=C3S2)CCO')


def test_exact_o_graft_and_single_oxidation_are_explicit():
    row=module.build(IDENTITY)
    assert row['reaction']['atom_balanced'] and row['reaction']['charge_balanced']
    assert row['states']==[['neutral',0,1,0],['oxidized',1,2,1]]
    assert not row['properties_computed']
    mol=Chem.MolFromSmiles(row['route_binding']['mapped_product_smiles'])
    maps={a.GetAtomMapNum():a.GetIdx() for a in mol.GetAtoms()};a,b=row['route_binding']['new_bond_atom_maps']
    assert {mol.GetAtomWithIdx(maps[a]).GetSymbol(),mol.GetAtomWithIdx(maps[b]).GetSymbol()}=={'C','O'}
    assert mol.GetBondBetweenAtoms(maps[a],maps[b]) is not None


@pytest.mark.parametrize('field,value',[('CID',137994),('InChIKey','wrong'),('ConnectivitySMILES','CCO')])
def test_rejects_mismatched_source_identity(field,value):
    with pytest.raises(ValueError):module.build({**IDENTITY,field:value})


def test_source_precursor_access_does_not_override_route_hold():
    route_spec=importlib.util.spec_from_file_location('ptz_route',ROOT/'scripts/pipeline/screen_grafting_routes.py')
    gate=importlib.util.module_from_spec(route_spec);route_spec.loader.exec_module(gate)
    row=module.build(IDENTITY);result=gate.screen(row,{row['id']:{'exact_precursor_access_reviewed':True}},{})
    assert result['route_status']=='hold' and not result['proceed_to_full_redox']
    assert 'attachment_selectivity_reviewed' in result['missing_route_evidence']


@pytest.mark.parametrize('directory',[
    'ptz_hydroxyethyl_source_2026-10-26','pmdi_hydroxyethyl_source_control_2026-10-26'])
def test_shared_all_axis_adapter_binds_product_and_stays_pending(directory,tmp_path):
    source=ROOT/'results/discovery'/directory
    assess_spec=importlib.util.spec_from_file_location('source_route_assess',ROOT/'scripts/pipeline/assess_source_route_candidate.py')
    assessor=importlib.util.module_from_spec(assess_spec);assess_spec.loader.exec_module(assessor)
    result=assessor.assess(source/'candidates.json',source/'route_screen.json',ROOT/'config/candidate_feasibility.json')
    assert result['decision']['status']=='pending' and len(result['decision']['pending_axes'])==10
    gate=json.loads((source/'route_screen.json').read_text());gate['candidates'][0]['evidence']['exact_design_key']='wrong'
    altered=tmp_path/'altered.json';altered.write_text(json.dumps(gate))
    with pytest.raises(ValueError,match='exact candidate'):
        assessor.assess(source/'candidates.json',altered,ROOT/'config/candidate_feasibility.json')
