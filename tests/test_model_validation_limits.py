"""A failed second-wave model cannot certify window capacity or dominate a shortlist."""
import importlib.util
import json
from pathlib import Path

import pytest

from redox.core.common import ROOT,load_config
from redox.screening import pareto
from redox.validation.model_limits import selection_admission

PTZ='CCN1c2ccc(OC)cc2Sc2cc(OC)ccc21'


@pytest.mark.parametrize('charges,role', [([0,1],'catholyte'),([0,-1,-2],'anolyte'),
                                       ([0,1,2],'anolyte')])
def test_first_wave_and_unrelated_electrode_paths_are_not_held(charges,role):
    assert selection_admission(PTZ,charges,role)['model_primary_admissible']


def test_second_oxidation_is_scoped_by_graph_and_declared_ladder():
    held=selection_admission(PTZ,[0,1,2],'catholyte')
    assert not held['model_primary_admissible']
    assert held['model_max_observed_error_V']==.547
    # Even a one-step path from an already oxidized resting state uses the failed step.
    assert not selection_admission(PTZ,[1,2],'catholyte')['model_primary_admissible']
    assert selection_admission('O=C1C=CC(=O)C=C1',[0,1,2],'catholyte')['model_primary_admissible']


@pytest.mark.parametrize('mutation',['stale_evidence','different_protocol','unknown_status'])
def test_unreviewed_method_or_evidence_cannot_silently_clear_the_scope(tmp_path,mutation):
    review=json.loads((ROOT/'config/model_validation_limits.json').read_text())
    if mutation=='stale_evidence':review['rules'][0]['evidence_sha256']='0'*64
    elif mutation=='different_protocol':review['rules'][0]['tested_energy_protocol']['basis']='other'
    else:review['rules'][0]['status']='independently_cleared'
    path=tmp_path/'review.json';path.write_text(json.dumps(review))
    assert not selection_admission(PTZ,[0,1,2],'catholyte',review_path=path)['model_primary_admissible']


def test_missing_review_fails_explicitly(tmp_path):
    with pytest.raises(RuntimeError,match='unavailable'):
        selection_admission(PTZ,[0,1,2],'catholyte',review_path=tmp_path/'missing')


def test_table_snapshot_tracks_model_review_backing_evidence(tmp_path,monkeypatch):
    from redox.workflows import tables
    monkeypatch.setattr(tables,'ROOT',tmp_path);monkeypatch.setattr(tables,'CALCS',tmp_path/'calcs')
    evidence=tmp_path/'evidence.json';evidence.write_text('{"residual":.547}')
    config=tmp_path/'config';config.mkdir()
    (config/'model_validation_limits.json').write_text(json.dumps(dict(rules=[dict(evidence_path='evidence.json')])))
    before=tables.snapshot_inputs(tmp_path/'results')
    assert str(evidence) in before
    evidence.write_text('{"residual":.1}')
    assert tables.snapshot_inputs(tmp_path/'results')[str(evidence)]!=before[str(evidence)]


def test_held_axes_cannot_give_primary_pareto_dominance():
    cfg=load_config('scorecard_config')
    good=dict(id='good',pool='catholyte',E_V=.2,sigma_E_V=.01,
              specific_capacity_mAh_g=50,lambda_i_eV=.4,sigma_lambda_eV=.1,
              disp_applicable=False,model_primary_admissible=True)
    held=dict(good,id='held',E_V=1.,specific_capacity_mAh_g=200,lambda_i_eV=.1,
              model_primary_admissible=False)
    front,objectives,dominated,missing=pareto._pareto_front([good,held],'catholyte',cfg)
    assert front==['good'] and 'model_validation' in missing['held']
    assert 'held' not in dominated['good']
    assert 'voltage' not in objectives['held'] and 'capacity' not in objectives['held']


def test_candidate_exports_keep_model_estimates_but_not_qualified_capacity():
    path=ROOT/'scripts/pipeline/collect_verified_candidates.py'
    spec=importlib.util.spec_from_file_location('model_scope_collector',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    row=dict(states=json.dumps([('neu',0,1,0),('ox1',1,2,1),('ox2',2,1,2)]),
             exact_grafted_smiles=PTZ,E1_vs_Fc_V=.249,E2_vs_Fc_V=1.202,
             capacity_nominal_mAh_g=180,accessible_capacity_mAh_g=90)
    module.apply_model_selection_gates([row])
    assert row['accessible_capacity_mAh_g'] is None
    assert row['model_accessible_capacity_estimate_mAh_g']==90
    assert row['capacity_nominal_mAh_g']==180 and row['E2_vs_Fc_V']==1.202
