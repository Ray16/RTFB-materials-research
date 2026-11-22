"""Build failures must not masquerade as a successful submission artifact."""
import importlib.util
from pathlib import Path
import subprocess

import pytest


ROOT = Path(__file__).resolve().parents[1]


def load_script(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts/manuscript" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("failed_pass", [0, 1, 2, 3])
def test_failed_compile_removes_stale_or_partial_pdf(monkeypatch, tmp_path, failed_pass):
    build = load_script("build_latex")
    monkeypatch.setattr(build, "OUT", tmp_path)
    pdf = tmp_path / "manuscript.pdf"
    pdf.write_bytes(b"stale PDF")
    commands = []

    def run(cmd, **kwargs):
        if not commands:
            assert not pdf.exists()
        commands.append(cmd)
        pdf.write_bytes(b"partial PDF")
        kwargs["stdout"].write("full diagnostic from compiler\n")
        if len(commands) - 1 == failed_pass:
            raise subprocess.CalledProcessError(1, cmd)

    monkeypatch.setattr(build.subprocess, "run", run)
    with pytest.raises(RuntimeError, match="LaTeX build failed"):
        build.compile_document()
    assert len(commands) == failed_pass + 1
    assert not pdf.exists()
    assert "full diagnostic" in (tmp_path / "build.log").read_text()


def test_missing_compiler_cannot_reuse_old_pdf(monkeypatch, tmp_path):
    build = load_script("build_latex")
    monkeypatch.setattr(build, "OUT", tmp_path)
    (tmp_path / "manuscript.pdf").write_bytes(b"stale PDF")

    def missing(*a, **kw):
        raise FileNotFoundError("xelatex")

    monkeypatch.setattr(build.subprocess, "run", missing)
    with pytest.raises(RuntimeError, match="LaTeX build failed at xelatex"):
        build.compile_document()
    assert not (tmp_path / "manuscript.pdf").exists()


def test_success_exit_without_pdf_is_failure(monkeypatch, tmp_path):
    build = load_script("build_latex")
    monkeypatch.setattr(build, "OUT", tmp_path)
    monkeypatch.setattr(build.subprocess, "run", lambda *a, **kw: None)
    with pytest.raises(RuntimeError, match="produced no PDF"):
        build.compile_document()


def test_si_reports_recorded_versions_without_executable_or_installed_packages(monkeypatch):
    import importlib.metadata

    build = load_script("build_si")
    monkeypatch.setenv("PATH", "")
    monkeypatch.setenv("XTB_BIN", "/nonexistent/xtb")
    monkeypatch.setattr(subprocess, "run", lambda *a, **kw: pytest.fail("ran an executable"))
    monkeypatch.setattr(importlib.metadata, "version", lambda *a: pytest.fail("read build environment"))
    section = build.s_software()
    recorded = (ROOT / "manuscript/supporting_information.md").read_text()
    assert section.strip() == ("## S9. Software" + recorded.split("## S9. Software", 1)[1]
                              .split("## S10.", 1)[0]).strip()


def test_raw_pareto_front_anolyte_and_catholyte():
    from redox.screening.pareto import raw_pareto_front

    anolytes = [
        {"id": "a_best_v", "E": -1.5, "cap": 140.0},
        {"id": "a_best_c", "E": -1.0, "cap": 190.0},
        {"id": "a_dom", "E": -1.2, "cap": 130.0},
        {"id": "a_mid", "E": -1.3, "cap": 160.0},
    ]
    assert [r["id"] for r in raw_pareto_front(anolytes, pool="anolyte")] == [
        "a_best_v",
        "a_best_c",
        "a_mid",
    ]

    catholytes = [
        {"id": "c_best_v", "E": 0.8, "cap": 140.0},
        {"id": "c_best_c", "E": 0.4, "cap": 170.0},
        {"id": "c_dom", "E": 0.3, "cap": 150.0},
    ]
    assert [r["id"] for r in raw_pareto_front(catholytes, pool="catholyte")] == [
        "c_best_v",
        "c_best_c",
    ]


def test_raw_pareto_front_rejects_invalid_pool():
    from redox.screening.pareto import raw_pareto_front

    with pytest.raises(ValueError, match="pool must be"):
        raw_pareto_front([], pool="ambipolar")


def test_scored_pool_has_no_placeholder_chemical_names():
    import csv

    rows = list(csv.DictReader((ROOT / "data/manuscript/full_pool_scored_legacyG_MeCN_TBAPF6.csv").open()))
    assert len(rows) >= 91
    placeholders = [r["id"] for r in rows if "canonical structure supplied" in r["name"]]
    assert placeholders == []


def test_latex_table_cell_allows_breaks_in_long_chemical_names():
    build = load_script("build_latex")
    tex = build.table(
        [
            "| Candidate | E (V) |",
            "|---|---|",
            "| N-(4-methylbenzyl)-3,7-dimethoxyphenothiazine | +0.78 |",
            "| 1-(4-methylbenzyloxy)anthraquinone | -1.48 |",
        ]
    )
    assert r"dimethoxy\-pheno\-thiazine" in tex
    assert r")\allowbreak{}anthra\-quinone" in tex
