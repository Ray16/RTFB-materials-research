"""The supported baseline cannot admit other viologen substitutions or ion cycles."""
import importlib.util
import json
from pathlib import Path

import pytest

from redox.build.grafting import grafting_reaction
from redox.screening.candidate_registry import publish


ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('tested_viologen_route',ROOT/'scripts/pipeline/screen_grafting_routes.py')
gate=importlib.util.module_from_spec(spec);spec.loader.exec_module(gate)


@pytest.fixture
def case():
    folder=ROOT/'validation/viologen_route_2026-10-30'
    load=lambda name:json.loads((folder/name).read_text())
    return load('candidate_rows.json'),load('reviewed_entries.json')['candidates']


def test_same_graph_aliases_have_exact_prospective_route_passes(case):
    rows,evidence=case
    for row in rows:
        result=gate.screen(row,evidence,{})
        assert result['route_status']=='proceed'
        assert result['exact_route_binding_verified']
        assert result['experimental_ease'].startswith('reviewed_workable_route')


def test_quaternization_releases_chloride_without_a_base_or_HCl(case):
    row=case[0][0]
    reaction=grafting_reaction(row['frag'],row['product_smiles'])
    assert reaction['atom_balanced'] and reaction['charge_balanced']
    assert (reaction['precursor_charge'],reaction['product_charge'],reaction['byproduct_charge'])==(1,2,-1)
    assert reaction['byproduct_smiles']=='[Cl-]'
    assert not reaction['base_correction_required']


def test_route_review_is_not_transferred_to_ethyl_cap_or_unreviewed_alias(case):
    rows,evidence=case
    other=dict(rows[0],id='unreviewed_viologen')
    assert gate.screen(other,evidence,{})['route_status']=='hold'
    other=dict(rows[0],frag='CC[n+]1ccc(-c2cc[n+]([*:1])cc2)cc1')
    other.pop('product_smiles');other.pop('precursor_smiles')
    result=gate.screen(other,evidence,{})
    assert result['route_status']=='hold'
    assert not result['exact_route_binding_verified']


@pytest.mark.parametrize('late_mutation',[False,True])
def test_registry_explicit_new_catalog_preserves_immutable_archive(tmp_path,case,monkeypatch,late_mutation):
    archive=tmp_path/'archive';archive.mkdir()
    row=case[0][0]
    (archive/'rows.json').write_text(json.dumps([row]))
    old=json.dumps({'candidates':{}})
    (archive/'routes.json').write_text(old)
    (archive/'assessment.json').write_text(json.dumps({'candidates':{}}))
    specification=dict(sources=[dict(name='enum',path='rows.json',role='enumerated')],
                       route_evidence='routes.json',assessment='assessment.json')
    plan=tmp_path/'spec.json';plan.write_text(json.dumps(specification))
    new=tmp_path/'reviewed.json';new.write_text(json.dumps({'candidates':case[1]}))
    before=publish(plan,archive,tmp_path/'before')
    if late_mutation:
        mkdir=Path.mkdir
        def begin_snapshot(path,*args,**kwargs):
            if path==tmp_path/'after':new.write_text(json.dumps({'candidates':{}}))
            return mkdir(path,*args,**kwargs)
        monkeypatch.setattr(Path,'mkdir',begin_snapshot)
    after=publish(plan,archive,tmp_path/'after',route_evidence=new)
    assert before['route_counts']=={'hold':1}
    assert after['route_counts']=={'proceed':1}
    assert (archive/'routes.json').read_text()==old
    frozen=json.loads((tmp_path/'after/inputs/route_evidence.json').read_text())
    assert frozen=={'candidates':case[1]}
    assert after['full_material_admissions']==0 and after['new_quantum_jobs']==0
