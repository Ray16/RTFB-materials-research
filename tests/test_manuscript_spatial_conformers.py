"""Geometry and conformer-gap claims must remain bound to accepted parents."""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location('spatial_claims', ROOT/'manuscript/check_claims.py')
CLAIMS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CLAIMS)


def load(name):
    path = (ROOT/'results/figures/polymer_models/spacing_conformers.json' if name == 'spacing_conformers.json'
        else ROOT/'manuscript/evidence/2026-10-23_deployment_revision'/name)
    return json.loads(path.read_text())


def test_saved_contact_distances_and_free_energy_gaps_reconstruct():
    report = CLAIMS.audit_spatial_conformers(load)
    assert report['checked_conformers'] == 4
    assert 11.8 < report['within_graph_deltaG_kJ_mol']['occ_2_3__pose1'] < 12
    assert 7 < report['within_graph_deltaG_kJ_mol']['occ_2_4__pose0'] < 7.2
    assert not report['isolated_motif_binding_inferred']


@pytest.mark.parametrize('change', ['distance', 'energy_parent', 'graph', 'missing_pose', 'convergence', 'reflection'])
def test_unbound_or_promoted_spatial_claims_are_refused(change):
    def changed(name):
        data = load(name)
        if name == 'spacing_conformers.json':
            if change == 'distance': data['states'][0]['contact']['centroid_distance_A'] += 1
            elif change == 'energy_parent': data['states'][0]['G_eV'] += .01
            elif change == 'graph': data['states'][0]['graph'] = 'occ_2_4'
            elif change == 'missing_pose': data['states'].pop()
            elif change == 'reflection': data['states'][0]['coordinate_rotation'][0] = [-v for v in data['states'][0]['coordinate_rotation'][0]]
            else: data['conformer_convergence_established'] = True
        return data
    with pytest.raises(ValueError): CLAIMS.audit_spatial_conformers(changed)


def test_centroid_reconstruction_is_independent_of_a_proper_rigid_view():
    def transformed(name):
        data = load(name)
        if name == 'spacing_conformers.json':
            for row in data['states']:
                row['core_coordinates_A'] = [[[-y+1,x+2,z+3] for x,y,z in core] for core in row['core_coordinates_A']]
        return data
    assert CLAIMS.audit_spatial_conformers(transformed)['checked_conformers'] == 4
