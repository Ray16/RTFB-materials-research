import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('dipea_experimental_anchor',ROOT/'scripts/validation/assess_dipea_experimental_anchor.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
EPOCH=ROOT/'manuscript/evidence/2026-10-26_native_replay'


def article(path,doi='10.1039/d3sc02611f',paragraph='NMR competition experiments estimate 18.7 for DIPEAHPF6 hexafluorophosphate in MeCN.'):
    path.write_text('<article><article-meta><article-id pub-id-type="doi">'+doi+'</article-id>'
        '<title-group><article-title>Controlling product selectivity during dioxygen reduction with Mn complexes</article-title>'
        '</title-group></article-meta><body><p>'+paragraph+'</p></body></article>')


def evaluate(path):
    return module.assess(path,EPOCH/'base_assisted.json',EPOCH/'acid_base.json',EPOCH/'basicity_reference_origin.json')


def test_signs_and_hybrid_sensitivity_reconstruct_without_promoting_a_reaction_prediction(tmp_path):
    path=tmp_path/'synthetic.xml';article(path);result=evaluate(path);c=result['comparison']
    assert c['pKa_implied_capture_kJ_mol']==pytest.approx(-47.94727982)
    assert c['DFT_minus_pKa_implied_capture_kJ_mol']==pytest.approx(-23.02815513)
    assert c['pKa_implied_DIPEA_DBU_transfer_kJ_mol']==pytest.approx(32.02193331)
    assert c['DFT_minus_pKa_implied_relative_kJ_mol']==pytest.approx(-.22035737)
    row=next(r for r in result['candidate_hybrid_sensitivities'] if r['id']=='exp_bq_oxy_methoxy')
    assert row['hybrid_capture_sensitivity_kJ_mol']==pytest.approx(9.68718072)
    assert not result['calibrated_full_route_prediction'] and not result['candidate_decisions_modified']
    assert next(r for r in result['candidate_hybrid_sensitivities'] if r['id']=='exp_aq_2oxy_6methoxy')['hybrid_capture_sensitivity_kJ_mol'] is None


@pytest.mark.parametrize('doi,text',[
    ('10.0000/other','NMR competition experiments estimate 18.7 for DIPEAHPF6 hexafluorophosphate in MeCN.'),
    ('10.1039/d3sc02611f','A predicted value of 18.7 is listed for DIPEA.'),
])
def test_wrong_source_or_unsupported_measurement_is_refused(tmp_path,doi,text):
    path=tmp_path/'synthetic.xml';article(path,doi,text)
    with pytest.raises(ValueError):evaluate(path)
