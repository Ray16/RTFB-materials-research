import importlib.util
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('dipea_uv_anchor',ROOT/'scripts/validation/assess_dipea_uv_anchor.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
EPOCH=ROOT/'manuscript/evidence/2026-10-26_dipea_nmr'


def source_files(tmp_path,monkeypatch,page):
    article=tmp_path/'synthetic.xml';article.write_text('<article><article-meta><article-id pub-id-type="doi">10.1039/d1ra04738h</article-id></article-meta><body><p>Protonation constants of DIPEA in MeCN were measured at 25 °C by bromocresol green.</p></body></article>')
    si=tmp_path/'synthetic.pdf';si.write_bytes(b'fixture')
    monkeypatch.setattr(module,'PdfReader',lambda path:SimpleNamespace(pages=[None]*42+[SimpleNamespace(extract_text=lambda:page)]))
    return article,si


def test_25c_experimental_relative_and_capture_comparisons_preserve_scope(tmp_path,monkeypatch):
    article,si=source_files(tmp_path,monkeypatch,'Figure S51. DIPEA in MeCN at 25.0 C. log KH = 18.64(1)')
    r=module.assess(article,si,EPOCH/'base_assisted.json',EPOCH/'acid_base.json',EPOCH/'basicity_reference_origin.json',EPOCH/'dipea_experimental_anchor.json')
    assert r['comparison']['pKa_implied_capture_kJ_mol']==pytest.approx(-47.60479925)
    assert r['comparison']['DFT_minus_pKa_implied_capture_kJ_mol']==pytest.approx(-23.37063570)
    assert r['comparison']['pKa_implied_relative_base_transfer_kJ_mol']==pytest.approx(32.36441388)
    assert r['comparison']['DFT_minus_pKa_implied_relative_kJ_mol']==pytest.approx(-.56283794)
    assert not r['calibrated_full_route_prediction'] and not r['candidate_decisions_modified']
    assert next(x for x in r['candidate_hybrid_sensitivities'] if x['id']=='exp_aq_2oxy_6methoxy')['hybrid_capture_sensitivity_kJ_mol'] is None


@pytest.mark.parametrize('page',[
    'Figure S41. DIPEA in DMSO at 25.0 C. log KH = 9.23(1)',
    'Figure S51. DIPEA in MeCN at 25.0 C. log KH = 18.6(1)',
    'Figure S51. DIPEA in DMSO at 25.0 C. log KH = 18.64(1)',
])
def test_wrong_figure_solvent_or_value_is_rejected(tmp_path,monkeypatch,page):
    article,si=source_files(tmp_path,monkeypatch,page)
    with pytest.raises(ValueError):module.extract_observation(article,si)
