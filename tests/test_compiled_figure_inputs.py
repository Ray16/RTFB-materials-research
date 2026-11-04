import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).resolve().parents[1]
definition=importlib.util.spec_from_file_location('checked_paper_figures',ROOT/'manuscript/build.py')
builder=importlib.util.module_from_spec(definition);definition.loader.exec_module(builder)


@pytest.mark.parametrize('change',[None,'changed_bytes','unbound_figure'])
def test_the_actual_compiled_figure_must_match_the_snapshot(tmp_path,change):
    paper=tmp_path/'manuscript';paper.mkdir();figures=tmp_path/'results/figures';figures.mkdir(parents=True)
    figure=figures/'actual.pdf';figure.write_bytes(b'checked figure')
    (paper/'manuscript.tex').write_text(r'\includegraphics[width=\linewidth]{../results/figures/actual.pdf}')
    (paper/'supporting_information.tex').write_text('no figures')
    snapshot={'results/figures/actual.pdf':builder.sha(figure)}
    if change=='changed_bytes':figure.write_bytes(b'changed after capture')
    if change=='unbound_figure':snapshot={}
    if change:
        with pytest.raises(ValueError):builder.compiled_figure_inputs(paper,tmp_path,snapshot)
    else:assert builder.compiled_figure_inputs(paper,tmp_path,snapshot)==snapshot
