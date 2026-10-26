import importlib.util
import json
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


native=load('tempo_native_test',ROOT/'scripts/validation/review_tempo_native_potential.py')
decision=load('tempo_decision_test',ROOT/'scripts/pipeline/assess_tempo_route_and_material.py')


def test_oxidation_reference_sign_and_common_energy_shift():
    assert native.potential(-10.,-15.,-20.,-24.) == pytest.approx(1.)
    assert native.potential(90.,85.,80.,76.) == pytest.approx(1.)
    with pytest.raises(ValueError,match='Finite'):
        native.potential(float('nan'),0.,0.,0.)


def test_route_proceed_does_not_admit_incomplete_material(tmp_path,monkeypatch):
    evidence={'exact':'source'}
    route={'candidate':{'product_smiles':'C'},'evidence':evidence}
    screen={'candidates':[{'id':'tempo','route_status':'proceed','proceed_to_full_redox':True,'evidence':evidence}]}
    potential={'status':'TEMPO_local_oxidation_native_partition_independently_reconstructed',
               'product_smiles':'C','product_inchikey':'test-key',
               'potential_V_vs_Fc_by_rotor_cutoff':{'50':0.58},'source_sha256':{}}
    paths=[]
    for name,value in [('route',route),('screen',screen),('potential',potential)]:
        path=tmp_path/(name+'.json');path.write_text(json.dumps(value));paths.append(path)
    result=decision.assess(*paths,ROOT/'config/candidate_feasibility.json')
    assert result['decision']['status']=='pending'
    assert not result['decision']['all_required_axes_pass']
    assert next(a for a in result['decision']['axes'] if a['axis']=='attachment_route')['status']=='pass'
    assert 'motif_association' in result['decision']['pending_axes']
