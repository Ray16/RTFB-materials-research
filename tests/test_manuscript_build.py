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
