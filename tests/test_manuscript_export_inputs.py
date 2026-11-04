import hashlib
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
specification = importlib.util.spec_from_file_location('paper_dependency_export', ROOT/'manuscript/export_overleaf.py')
exporter = importlib.util.module_from_spec(specification)
specification.loader.exec_module(exporter)


def setup(tmp_path, contents):
    paper = tmp_path/'paper'
    evidence = paper/'evidence/epoch'
    evidence.mkdir(parents=True)
    hashes = {}
    for name, text in contents.items():
        path = evidence/name
        path.write_text(text)
        hashes[str(path.relative_to(tmp_path))] = hashlib.sha256(path.read_bytes()).hexdigest()
    return paper, evidence, dict(snapshot_sha256=hashes), tmp_path/'export'


def test_nested_frozen_numeric_inputs_are_copied_and_rewritten(tmp_path):
    paper, evidence, snapshot, destination = setup(tmp_path, {
        'values.tex': r'\input{evidence/epoch/constants.tex}',
        'constants.tex': r'\newcommand{\Value}{-1.23}'})
    text = exporter.export_dependencies(r'\input{evidence/epoch/values.tex}', paper, tmp_path,
                                        evidence, '../figures/', snapshot, destination)
    assert text == r'\input{inputs/values.tex}'
    assert (destination/'inputs/values.tex').read_text() == r'\input{inputs/constants.tex}'
    assert (destination/'inputs/constants.tex').read_text() == r'\newcommand{\Value}{-1.23}'


def test_changed_numeric_fragment_cannot_enter_export(tmp_path):
    paper, evidence, snapshot, destination = setup(tmp_path, {'values.tex': 'checked numbers'})
    (evidence/'values.tex').write_text('changed numbers')
    with pytest.raises(ValueError, match='not bound'):
        exporter.export_dependencies(r'\input{evidence/epoch/values.tex}', paper, tmp_path,
                                     evidence, '../figures/', snapshot, destination)


def test_cyclic_inputs_fail_instead_of_recursing_indefinitely(tmp_path):
    paper, evidence, snapshot, destination = setup(tmp_path, {
        'a.tex': r'\input{evidence/epoch/b.tex}', 'b.tex': r'\input{evidence/epoch/a.tex}'})
    with pytest.raises(ValueError, match='Cyclic'):
        exporter.export_dependencies(r'\input{evidence/epoch/a.tex}', paper, tmp_path,
                                     evidence, '../figures/', snapshot, destination)


def test_include_keeps_tex_automatic_extension(tmp_path):
    paper, evidence, snapshot, destination = setup(tmp_path, {'chapter.tex': 'frozen chapter'})
    result = exporter.export_dependencies(r'\include{evidence/epoch/chapter}', paper, tmp_path,
                                          evidence, '../figures/', snapshot, destination)
    assert result == r'\include{inputs/chapter}'
    assert (destination/'inputs/chapter.tex').read_text() == 'frozen chapter'
