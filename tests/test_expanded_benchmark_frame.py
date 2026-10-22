import importlib.util
import hashlib
import json
from pathlib import Path

import pytest


def test_current_benchmark_frame_requires_unchanged_checked_artifact(tmp_path):
    script = Path(__file__).resolve().parents[1] / "scripts/plotting/candidates/research_presentation.py"
    spec = importlib.util.spec_from_file_location("expanded_benchmark_frame", script)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    folder = tmp_path / "slides/generated"
    folder.mkdir(parents=True)
    (folder / "matched_parent_validation_plot.tex").write_text("historical subset")
    assert module.benchmark_frame(tmp_path) == "historical subset"
    include = folder / "retrospective_redox_benchmark.tex"
    include.write_text("expanded checked panel")
    metadata = tmp_path / "results/figures/validation/retrospective_redox_benchmark.json"
    metadata.parent.mkdir(parents=True)
    metadata.write_text(json.dumps(dict(status="checked_retrospective_snapshot", output_sha256={
        str(include.relative_to(tmp_path)): hashlib.sha256(include.read_bytes()).hexdigest()})))
    assert module.benchmark_frame(tmp_path) == "expanded checked panel"
    include.write_text("changed artifact")
    with pytest.raises(ValueError, match="artifact changed"):
        module.benchmark_frame(tmp_path)
