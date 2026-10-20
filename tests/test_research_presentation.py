"""The main narrative stays tied to its data and avoids stale scientific claims."""
import importlib.util
import json
from pathlib import Path
import pytest


@pytest.fixture
def talk(tmp_path):
    source=Path(__file__).resolve().parents[1]/'scripts/plotting/candidates/research_presentation.py'
    spec=importlib.util.spec_from_file_location('research_talk',source)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    figures=tmp_path/'results/candidates/figures'
    (figures/'metadata').mkdir(parents=True);(figures/'data').mkdir()
    manifest=dict(snapshot='audited-A',candidates=6,admitted_potential_candidates=3)
    (figures/'metadata/manifest.json').write_text(json.dumps(manifest))
    rows=[dict(id=gid,family=family,capacity_nominal_mAh_g=capacity) for gid,family,capacity in
          [('lead-a','quinone (n-type)',205.),('lead-b','quinone (n-type)',195.),('lead-c','quinone (n-type)',160.),
           ('ndi','imide (n-type)',135.),('pmdi','imide (n-type)',154.),('viologen','pyridine-multi-e',185.)]]
    highlights=dict(highlighted_ids=['lead-a','lead-b','lead-c'],matched_contrasts=[])
    (figures/'data/highlights.json').write_text(json.dumps(highlights))
    proposals=tmp_path/'results/discovery/vinyl_monomers';proposals.mkdir(parents=True)
    records=[dict(id='vm_'+gid,source_capped_id=source,capacity_nominal_mAh_g=capacity) for gid,source,capacity in
             [('design_ndi_ethyl','ndi',131.),('design_pmdi_ethyl','pmdi',149.),('ethylviologen','viologen',177.)]]
    (proposals/'proposals.json').write_text(json.dumps(dict(records=records)))
    generated=tmp_path/'slides/generated';generated.mkdir(parents=True)
    (generated/'matched_parent_validation_plot.tex').write_text(r'\begin{frame}{Validation}Plot\end{frame}')
    return module,tmp_path,dict(snapshot_path='audited-A',original39=rows,additional_candidates=[])


def test_main_story_has_background_workflow_method_and_current_results(talk):
    module,root,payload=talk
    result=module.research_frames(payload,root)
    assert result.count(r'\begin{frame}')==13
    assert result.index('The target is a redox-active polymer')<result.index('Workflow:')<result.index('Method:')<result.index('Quinone models lead')
    assert '205.0' in result and 'Preformed polymer' in result and 'at polymer sites' in result
    assert 'vinyl' not in result.lower()
    assert 'Degree of grafting' in result and 'preformed polymer' in result
    assert 'arithmetic checks' not in result and 'GPU' not in result
    payload['snapshot_path']='audited-B'
    with pytest.raises(RuntimeError,match='same audited snapshot'):module.research_frames(payload,root)


def test_a_new_leader_cannot_inherit_a_quinone_or_topology_claim(talk):
    module,root,payload=talk
    payload['original39'][0]['family']='imide (n-type)'
    result=module.research_frames(payload,root)
    assert 'Quinone models lead' not in result
    assert 'The two BQ models reach' not in result
    assert 'Substituted BQ models lead' not in result
    assert 'NDI core topology provides' not in result


def test_completed_base_route_replaces_stale_bq_hold(talk):
    module,root,payload=talk
    folder=root/'results/grafting';folder.mkdir(parents=True)
    decisions={key:dict(status='base_required',thermodynamic_status='favorable_base_assisted_model',
                       base_assisted_deltaG_kJ_mol={'dbu':{'free_ions_kJ_mol':value}})
               for key,value in [('exp_bq_oxy_methoxy',-45.1),('exp_bq_oxy_trimethyl',-46.7)]}
    (folder/'selection_assessment.json').write_text(json.dumps({'candidates':decisions}))
    result=module.research_frames(payload,root)
    assert 'Base required; BQ candidates retained for evaluation.' in result
    assert '-45.1/-46.7' in result
    assert 'selectivity and yield untested' in result
    assert 'on hold' not in result.lower().replace('unknown routes stay on hold','')
