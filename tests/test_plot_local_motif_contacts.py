import importlib.util
import json
from pathlib import Path

import pytest


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('local_motif_plot',ROOT/'scripts/plotting/polymer_models/plot_local_motif_contacts.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
SOURCE=ROOT/'results/reviews/2026-10-26/local_motif_complete/summary.json'


def test_checked_points_reconstruct_the_four_geometries():
    _,points=module.validated_points(SOURCE)
    assert len(points)==4
    by_key={p['parent_key']:p for p in points}
    assert by_key['occ_2_3__pose0']['interaction_kJ_mol']==pytest.approx(-24.7615245734)
    assert by_key['occ_2_4__pose1']['interaction_kJ_mol']==pytest.approx(-26.2195879462)
    assert all(p['cap_range_kJ_mol'][0]<=p['interaction_kJ_mol']<=p['cap_range_kJ_mol'][1] for p in points)


def test_changed_electronic_term_is_rejected(tmp_path):
    data=json.loads(SOURCE.read_text());data['states'][0]['exact_state_energies_eV']['AB']['acetonitrile']+=0.01
    changed=tmp_path/'changed.json';changed.write_text(json.dumps(data))
    with pytest.raises(ValueError,match='does not reconstruct'):
        module.validated_points(changed)


def test_free_energy_claim_is_rejected(tmp_path):
    data=json.loads(SOURCE.read_text());data['free_energies_computed']=True
    changed=tmp_path/'changed.json';changed.write_text(json.dumps(data))
    with pytest.raises(ValueError,match='scope differs'):
        module.validated_points(changed)
