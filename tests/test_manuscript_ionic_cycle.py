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


def full_fixture(tmp_path):
    folder=fixture(tmp_path)
    data=json.loads((ROOT/'validation/ionic_cycle_completion_2026-11-02/complete_cycles_12.json').read_text())
    (folder/'ionic_native_cycles_v4.json').write_text(json.dumps(data))
    module=checker()
    (folder/'ionic_native_cycle_values_v4.tex').write_text(module.values(data))
    (folder/'ionic_native_cycle_table_v4.tex').write_text(module.table(data))
    return folder


def test_all_four_cycles_reconstruct_and_survive_compilation_with_their_signs(tmp_path):
    module=checker();report=module.audit(full_fixture(tmp_path))
    assert report['checked_components']==12 and report['complete_finite_cycles']==4 and report['pending_cycles']==0
    values=[r['G_kJ_mol'] for r in report['primary_cycle_values']]
    assert [v>0 for v in values]==[False,True,True,False]
    module.check_compiled_values('−3.94 22.58 8.28 −5.30',report)
    with pytest.raises(ValueError):module.check_compiled_values('−3.94 22.58 8.28 5.30',report)


@pytest.mark.parametrize('bad',['component','cycle','scope','frames','model'])
def test_the_recovered_fourth_component_cannot_be_promoted_with_incomplete_or_mixed_evidence(tmp_path,bad):
    folder=full_fixture(tmp_path);path=folder/'ionic_native_cycles_v4.json';data=json.loads(path.read_text())
    row=data['records']['alternating__occ_2_3']
    if bad=='component':row['G_by_cutoff_eV']['50']+=.1
    if bad=='cycle':data['cycles'][2]['delta_delta_G_by_cutoff_kJ_mol']['50']=0.
    if bad=='scope':data['conformer_ensemble_qualified']=True
    if bad=='frames':row['actual_SCF_checkpoint_frames_checked']=24
    if bad=='model':row['context']['solvent_model']={'name':'borrowed primary molecular model'}
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):checker().audit(folder)
