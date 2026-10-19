"""Slides cannot combine a current table with stale or changed comparison images."""
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest


@pytest.fixture
def collection(tmp_path):
    path=Path(__file__).resolve().parents[1]/'scripts/pipeline/update_verified_slides.py'
    spec=importlib.util.spec_from_file_location('comparison_slides',path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    folder=tmp_path/'results/candidates/figures';folder.mkdir(parents=True)
    slides=tmp_path/'slides';(slides/'generated').mkdir(parents=True)
    source=slides/'redox_workflow.tex';source.write_text(r'\input{generated/candidate_comparison_collection.tex}')
    hashes={}
    for name in ('redox_capacity.png','synthesis_tradeoff.png','new_candidate_shifts.png','candidate_structures.png','candidate_comparison_collection.tex'):
        data=name.encode();relative=module.artifact_path(name);target=folder/relative
        target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(data);hashes[relative]=hashlib.sha256(data).hexdigest()
    (slides/'generated/candidate_comparison_collection.tex').write_bytes(b'candidate_comparison_collection.tex')
    (folder/'metadata/manifest.json').write_text(json.dumps(dict(snapshot='snapshot-A',files=hashes)))
    return module,tmp_path,source,folder


def test_current_collection_binds_all_three_figures_and_the_include(collection):
    module,root,source,_=collection
    inputs=module.comparison_collection_inputs(dict(snapshot_path='snapshot-A'),root,source)
    assert len(inputs)==6


def test_stale_figure_snapshot_cannot_publish_with_new_table(collection):
    module,root,source,_=collection
    with pytest.raises(RuntimeError,match='displayed candidate snapshot'):
        module.comparison_collection_inputs(dict(snapshot_path='snapshot-B'),root,source)


def test_changed_figure_is_rejected_before_slide_build(collection):
    module,root,source,folder=collection
    (folder/'plots/redox_capacity.png').write_bytes(b'changed')
    with pytest.raises(RuntimeError,match='differs from its checked collection'):
        module.comparison_collection_inputs(dict(snapshot_path='snapshot-A'),root,source)
