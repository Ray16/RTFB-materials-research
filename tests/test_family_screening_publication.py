"""Preoptimization progress must remain separate from admitted solution states."""
import importlib.util
from pathlib import Path


def test_screening_frame_does_not_promote_uma_counts_to_native_counts():
    root=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('family_slide_publisher',root/'scripts/pipeline/update_verified_slides.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    row=dict(search_core='viologen',search_contrast='ring_CN',declared_states=3,
             uma_gas_verified_states=3,native_solution_ready_states=0)
    frame='\n'.join(module.family_screening_frame([row],'as of test'))
    assert '3/3 & 0/3' in frame and 'gas-phase screening' in frame
    assert 'Solution potentials and usable electrons require' in frame
    assert module.family_screening_frame([dict(id='historical')],'as of test')==[]


def test_publisher_tracks_partial_comparison_and_control_registry_aliases(tmp_path,monkeypatch):
    root=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('partial_scope_publisher',root/'scripts/pipeline/update_verified_slides.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    gids=['good','partial_registry','control_registry','unused','ferrocene']
    summary=dict(arithmetic_discrepancies=0,inputs_changed_during_snapshot=[],spin_input_gids=gids,
                 spin_input_paths_by_gid={gid:[] for gid in gids})
    payload=dict(summary=summary,original39=[dict(id='good',registry_id='good'),
        dict(id='partial_alias',registry_id='partial_registry')],additional_candidates=[],
        comparison_controls=[dict(id='control_alias',registry_id='control_registry')],
        matched_comparisons=[dict(id='partial_alias',baseline='control_alias',core='ndi',delta_E1_V=.1,delta_E2_V=None)])
    monkeypatch.setattr(module,'frames',lambda *args:('generated frames',['good']))
    observed=[]
    monkeypatch.setattr(module,'input_inventory_changed',lambda scoped:observed.append(set(scoped['spin_input_gids'])) or True)
    assert module.publish(payload,tmp_path)['status']=='waiting_for_current_input_inventory'
    assert observed==[{'good','partial_registry','control_registry','ferrocene'}]


def test_family_shift_frame_preserves_partial_events_and_omits_other_pilots():
    root=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('paired_slide_publisher',root/'scripts/pipeline/update_verified_slides.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    rows=[dict(core='pmdi',contrast='ring_F',delta_E1_V=.15,delta_E2_V=None,
               delta_wave_spacing_V=None),
          dict(core=None,contrast=None,delta_E1_V=.7,delta_E2_V=.8)]
    frame='\n'.join(module.paired_comparison_frame(rows,'as of test'))
    assert 'PMDI & Ring F & 0.150 & -- & --' in frame
    assert '0.700' not in frame and 'None' not in frame
    assert module.paired_comparison_frame([dict(core='ndi')],'as of test')==[]


def test_expanded_screening_and_matched_results_keep_every_row_across_pages():
    root=Path(__file__).resolve().parents[1]
    spec=importlib.util.spec_from_file_location('expanded_slide_publisher',root/'scripts/pipeline/update_verified_slides.py')
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    rows=[dict(search_core='ndi',search_contrast='ring_OMe',declared_states=3,
               uma_gas_verified_states=3,native_solution_ready_states=0) for _ in range(18)]
    frame='\n'.join(module.family_screening_frame(rows,'as of test'))
    assert frame.count(r'\begin{frame}')==2 and frame.count('NDI & Ring OMe & 3/3 & 0/3')==18
    pairs=[dict(core='ndi',contrast='ring_Me',delta_E1_V=.1,delta_E2_V=None,
                delta_wave_spacing_V=None) for _ in range(18)]
    frame='\n'.join(module.paired_comparison_frame(pairs,'as of test'))
    assert frame.count(r'\begin{frame}')==2 and frame.count('NDI & Ring Me & 0.100 & -- & --')==18
