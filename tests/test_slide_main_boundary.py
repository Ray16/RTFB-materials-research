import importlib.util
from pathlib import Path
import pytest

PATH=Path(__file__).resolve().parents[1]/'scripts/pipeline/update_verified_slides.py'
spec=importlib.util.spec_from_file_location('slide_boundary_check',PATH);module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_preserved_extra_main_inputs_and_nested_benchmark_count(tmp_path):
    main=tmp_path/'main.tex';generated=tmp_path/'generated';generated.mkdir()
    (generated/'story.tex').write_text('\\begin{frame}{Story}\\end{frame}\n\\input{generated/benchmark.tex}\n')
    (generated/'benchmark.tex').write_text('\\begin{frame}{Benchmark}\\end{frame}\n')
    (generated/'extra.tex').write_text('\\begin{frame}{Extra}\\end{frame}\n')
    main.write_text('\\begin{document}\n\\input{generated/story.tex}\n\\input{generated/extra.tex}\n% \\begin{frame}{not a slide}\n\\appendix\n\\begin{frame}{Backup}\\end{frame}\n')
    count,dependencies=module.declared_main_story(main)
    assert count==3 and len(dependencies)==3


def test_cyclic_main_inputs_refused(tmp_path):
    main=tmp_path/'main.tex';main.write_text('\\begin{document}\n\\input{cycle.tex}\n\\appendix')
    (tmp_path/'cycle.tex').write_text('\\input{cycle.tex}')
    with pytest.raises(RuntimeError,match='Cyclic'):module.declared_main_story(main)
