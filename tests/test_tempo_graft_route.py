import importlib.util
import json
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('tempo_route_test',ROOT/'scripts/pipeline/review_tempo_graft_route.py')
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def source(tmp_path):
    patent=tmp_path/'patent.html'
    patent.write_text('<p>4-oxy-bound to a Merrifield polymer</p><p>Example A4 Preparation of Compound 12 '
        '4-hydroxy-TEMPO NaH Merrifield polymer 1% divinylbenzene 1.7 mmol of Cl/g '
        '0.9 mmol of N-oxyl/g DMF or THF 1 to 4 days Examples B:</p>')
    identity=tmp_path/'identity.json'
    identity.write_text(json.dumps({'PropertyTable':{'Properties':[{'CID':137994,'InChIKey':'UZFMOKQJFYMBGY-UHFFFAOYSA-N'}]}}))
    return patent,identity


def test_existing_tempo_graph_and_NaH_stoichiometry(tmp_path):
    patent,identity=source(tmp_path)
    result=module.build(patent,identity,ROOT/'library/manifest.csv')
    assert result['candidate']['id']=='tempo'
    assert result['NaH_overall']['atom_balanced'] and result['NaH_overall']['radical_electrons_balanced']
    assert result['NaH_overall']['radical_electrons']==1
    assert not result['full_material_eligibility']


def test_unbalanced_species_and_missing_source_are_rejected(tmp_path):
    with pytest.raises(ValueError,match='preserve'):
        module.balanced_smiles(['C'],['CC'])
    patent,identity=source(tmp_path);patent.write_text('<p>unrelated source</p>')
    with pytest.raises(ValueError,match='Patent lacks'):
        module.build(patent,identity,ROOT/'library/manifest.csv')
