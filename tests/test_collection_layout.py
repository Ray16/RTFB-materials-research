"""Generated outputs stay organized, including future formats and link rewrites."""
import importlib.util
from pathlib import Path
import pytest


def module():
    path=Path(__file__).resolve().parents[1]/'scripts/plotting/candidates/collection_layout.py'
    spec=importlib.util.spec_from_file_location('figure_collection_layout',path)
    result=importlib.util.module_from_spec(spec);spec.loader.exec_module(result);return result


def test_output_categories_are_explicit():
    m=module()
    assert m.artifact_path('redox_capacity.svg')=='plots/redox_capacity.svg'
    assert m.artifact_path('redox_states.pdf')=='structures/redox_states.pdf'
    assert m.artifact_path('coverage.png')=='supplementary/coverage.png'
    assert m.artifact_path('all_candidates.csv')=='data/all_candidates.csv'
    assert m.artifact_path('manifest.json')=='metadata/manifest.json'
    with pytest.raises(ValueError):m.artifact_path('unclassified.png')


def test_links_rewrite_only_known_local_generated_files():
    m=module();page='<a href="redox_capacity.pdf">Plot</a><img src="redox_states.svg"><a href="https://example.com">Source</a>'
    rewritten=m.route_html(page,['redox_capacity.pdf','redox_states.svg'])
    assert 'href="plots/redox_capacity.pdf"' in rewritten
    assert 'src="structures/redox_states.svg"' in rewritten
    assert 'href="https://example.com"' in rewritten


def test_staging_keeps_entry_point_and_routes_manifest_created_last(tmp_path):
    m=module()
    (tmp_path/'index.html').write_text('<a href="manifest.json">Audit</a><img src="redox_capacity.png">')
    (tmp_path/'redox_capacity.png').write_bytes(b'unchanged plot bytes')
    m.organize_stage(tmp_path)
    assert (tmp_path/'plots/redox_capacity.png').read_bytes()==b'unchanged plot bytes'
    assert not (tmp_path/'redox_capacity.png').exists()
    assert 'href="metadata/manifest.json"' in (tmp_path/'index.html').read_text()
    assert 'src="plots/redox_capacity.png"' in (tmp_path/'index.html').read_text()
