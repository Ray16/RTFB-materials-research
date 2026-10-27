import importlib.util
import json
import shutil
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('contact_slide_build',ROOT/'scripts/presentations/build_colleague_workflow.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_contact_slide_checks_all_four_source_bound_points():
    values,_=module.validate_local_contacts()
    assert values==pytest.approx([-24.7615245734,-26.2195879462])


def test_changed_separated_pose_cannot_be_presented(tmp_path,monkeypatch):
    paths=['results/reviews/2026-10-26/local_motif_complete/summary.json',
           'results/figures/polymer_models/local_motif_contacts.json',
           'results/figures/polymer_models/local_motif_contacts.pdf',
           'scripts/plotting/polymer_models/plot_local_motif_contacts.py']
    for relative in paths:
        target=tmp_path/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(ROOT/relative,target)
    metadata=tmp_path/paths[1];data=json.loads(metadata.read_text())
    next(p for p in data['points'] if p['parent_key']=='occ_2_3__pose1')['interaction_kJ_mol']+=1
    metadata.write_text(json.dumps(data));monkeypatch.setattr(module,'ROOT',tmp_path)
    with pytest.raises(ValueError,match='plotted point differs'):
        module.validate_local_contacts()
