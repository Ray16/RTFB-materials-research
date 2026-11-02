"""Actual frozen ionic sums and scope cannot silently drift in the manuscript."""
import importlib.util
import json
import shutil
import pytest
from redox.core.common import ROOT


def checker():
    path=ROOT/'manuscript/check_ionic_native_cycles.py'
    spec=importlib.util.spec_from_file_location('tested_ionic_paper_cycle',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module);return module


def fixture(tmp_path):
    folder=ROOT/'manuscript/evidence/2026-11-02_first_ionic_G_cycle'
    for name in ('ionic_native_cycles.json','ionic_native_cycle_preparation.json','ionic_native_cycle_values.tex','ionic_native_cycle_table.tex'):
        shutil.copyfile(folder/name,tmp_path/name)
    return tmp_path


def test_current_ionic_cycle_reconstructs_with_all_pending_partners_explicit(tmp_path):
    result=checker().audit(fixture(tmp_path))
    assert result['complete_finite_cycles']==1 and result['pending_cycles']==3
    assert abs(result['delta_delta_G_by_cutoff_kJ_mol']['50']+3.9365010618683995)<1e-9
    assert not result['conformer_ensemble_qualified'] and not result['electronic_coupling_computed']


@pytest.mark.parametrize('bad',['component','cycle','scope','charge','macro','zero_for_missing'])
def test_corrupted_scalar_identity_or_scope_cannot_enter_the_paper(tmp_path,bad):
    folder=fixture(tmp_path);path=folder/'ionic_native_cycles.json';data=json.loads(path.read_text())
    if bad=='component':data['records']['all_same__occ_2_3']['G_by_cutoff_eV']['50']+=.1
    if bad=='cycle':data['cycles'][0]['delta_delta_G_by_cutoff_kJ_mol']['50']=0.
    if bad=='scope':data['conformer_ensemble_qualified']=True
    if bad=='charge':data['records']['all_same__occ_2_3']['charge']=2
    if bad=='macro':(folder/'ionic_native_cycle_values.tex').write_text('wrong numbers')
    if bad=='zero_for_missing':data['cycles'][1]['delta_delta_G_by_cutoff_kJ_mol']={'50':0.}
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):checker().audit(folder)
