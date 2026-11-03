import importlib.util
import json
from pathlib import Path
import shutil

import pytest

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('ptz_checker',ROOT/'manuscript/check_phenothiazine_redox.py')
checker=importlib.util.module_from_spec(spec);spec.loader.exec_module(checker)


@pytest.fixture
def evidence(tmp_path):
    source=ROOT/'validation/phenothiazine_cds_2026-11-03/publication'
    for p in source.iterdir():shutil.copyfile(p,tmp_path/p.name)
    return tmp_path


def test_checked_actual_publication_replays(evidence):
    assert checker.audit(evidence)['primary_V']==pytest.approx(.6260868628669414)


@pytest.mark.parametrize('change',['potential','scope','review','tex','model'])
def test_changed_frozen_claims_are_rejected(evidence,change):
    path=evidence/'ptz_redox_cycle.json';data=json.loads(path.read_text())
    if change=='potential':data['cycles']['minus']['potential_by_cutoff_V']['50']+=.2
    if change=='scope':data['experimental_accuracy_validated']=True
    if change=='model':data['cycles']['minus']['components']['fc_ox']['context']['solvent']='water'
    if change=='review':
        p=evidence/'ptz_review_minus.json';v=json.loads(p.read_text());v['G_by_cutoff_eV']['50']+=.1
        p.write_text(json.dumps(v))
    if change=='tex':(evidence/'ptz_redox_values.tex').write_text('invented potential')
    path.write_text(json.dumps(data))
    with pytest.raises(ValueError):checker.audit(evidence)
