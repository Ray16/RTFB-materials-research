"""Partial descriptor plots must not masquerade as a full kinetic/polymer ranking."""
import importlib.util
from pathlib import Path


def module():
    path=Path(__file__).resolve().parents[1]/'scripts/plotting/candidates/candidate_highlights.py'
    spec=importlib.util.spec_from_file_location('candidate_plot_highlights',path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result)
    return result


def row(gid,first=-1.,second=-2.,capacity=150.,**overrides):
    result=dict(id=gid,states=[['neu',0,1,0],['red1',-1,2,-1],['red2',-2,1,-2]],
        n_e_nominal=2,model_primary_admissible=True,solution_potentials_status='verified native calculation',
        E1_vs_Fc_V=first,E2_vs_Fc_V=second,capacity_nominal_mAh_g=capacity,
        accessible_capacity_mAh_g=capacity,SA_grafted=2.,lambda_i_eV=None)
    result.update(overrides);return result


def test_both_potentials_and_capacity_determine_the_partial_front():
    report=module().highlights([row('best'),row('dominated',first=-.5,second=-1.,capacity=100.),row('capacity',first=-.8,capacity=200.)])
    assert set(report['partial_front_ids'])=={'best','capacity'}
    assert not report['qualified_primary_Pareto_front_published']
    assert report['qualified_solution_lambda_count']==0


def test_catholytes_model_holds_and_one_electron_paths_cannot_enter_anolyte_front():
    rows=[row('ok'),row('cathode',states=[['neu',0,1,0],['ox1',1,2,1],['ox2',2,1,2]]),
          row('hold',model_primary_admissible=False),row('one_step',accessible_capacity_mAh_g=75.)]
    assert module().highlights(rows)['eligible_ids']==['ok']


def test_SA_proxy_cannot_create_dominance_or_a_polymer_synthesis_ranking():
    report=module().highlights([row('easy',SA_grafted=1.),row('hard',SA_grafted=6.)])
    assert set(report['partial_front_ids'])=={'easy','hard'}
    assert report['SA_is_annotation_not_objective'] and not report['polymer_synthesis_ranking']


def test_future_admitted_candidate_can_replace_the_highlighted_leader():
    report=module().highlights([row('old'),row('new',first=-1.2,second=-2.1,capacity=180.)])
    assert report['highlighted_ids']==['new']


def test_matched_control_must_be_admitted_and_both_events_lower():
    m=module();candidate=row('search_ndi_new',first=-1.1,second=-2.1,search_core='ndi');control=row('control')
    payload=dict(matched_comparisons=[dict(id=candidate['id'],baseline=control['id'])])
    assert not m.matched_highlights(payload,[candidate])
    assert len(m.matched_highlights(payload,[candidate,control]))==1
    candidate['E2_vs_Fc_V']=-1.9
    assert not m.matched_highlights(payload,[candidate,control])


def test_legend_collision_with_axis_text_blocks_figure_publication():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import pytest
    m=module();fig,ax=plt.subplots();ax.set_xlabel('Potential (V)')
    line,=ax.plot([0,1],[0,1],label='Candidate')
    fig.canvas.draw();box=ax.xaxis.label.get_window_extent(fig.canvas.get_renderer())
    center=((box.x0+box.x1)/(2*fig.bbox.width),(box.y0+box.y1)/(2*fig.bbox.height))
    legend=fig.legend(handles=[line],loc='center',bbox_to_anchor=center)
    with pytest.raises(ValueError,match='Legend overlaps axis text'):
        m.assert_clear_legend(fig,legend,[ax])
    plt.close(fig)
