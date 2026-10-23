"""Workspace recovery must preserve evidence and cannot hide physical failure."""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def case(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location("postpilot_prepare", ROOT / "scripts/pipeline/prepare_polymer_postpilot_continuation.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    parent = tmp_path / "parent"
    raw = tmp_path / "raw"
    key = "existing_neutral_pose0"
    for path in (parent / "controls", parent / "energy/pilot", parent / "energy/failures",
                 parent / "live_inputs" / key, raw / "energy" / (key + "_CDS.evaluations")):
        path.mkdir(parents=True, exist_ok=True)
    xyz = parent / "live_inputs" / key / "input.xyz"
    xyz.write_text("2\nexisting input\nH 0 0 0\nH 0 0 .74\n")
    task = dict(key=key, charge=0, mult=1)
    ready = dict(task, input_xyz=str(xyz), input_sha256={str(xyz): module.api.digest(xyz)})
    (xyz.parent / "task.json").write_text(json.dumps(ready))
    pilot = dict(status="memory_pilot_passed", observed_peak_GPU_MiB=1000,
                 output_sha256={}, pid=42, started=1234, finished=1235)
    (parent / "energy/pilot/status.json").write_text(json.dumps(pilot))
    (parent / "energy/pilot/energy_and_gradient.json").write_text(json.dumps(dict(scope="fixture point only")))
    (parent / "energy/failures" / (key + ".json")).write_text(json.dumps(dict(error="FileExistsError(17, 'File exists')")))
    plan = dict(tasks=[task], input_sha256={str(xyz): module.api.digest(xyz)},
                producer_sha256={}, raw_root=str(raw), maximum_memory_pilot_MiB=2000,
                quantum_settings=dict(protocol="unchanged_fixture"))
    (parent / "plan.json").write_text(json.dumps(plan))
    monkeypatch.setattr(module, "verify_stopped", lambda p: [dict(host="fixture", owners=[], matching_processes=[])])
    return module, parent, tmp_path / "successor", key, plan


def test_checked_pilot_is_explicitly_reused_without_moving_its_evidence(case):
    module, parent, out, key, original = case
    before = {str(p): p.read_bytes() for p in parent.rglob("*") if p.is_file()}
    module.prepare(parent, out)
    assert all(Path(p).read_bytes() == content for p, content in before.items())
    new = module.api.read_plan(out)
    assert new["tasks"] == original["tasks"] and new["quantum_settings"] == original["quantum_settings"]
    assert new["raw_root"] != original["raw_root"]
    copied = json.loads((out / "energy/pilot/status.json").read_text())
    assert copied["pid"] == 42 and copied["finished"] == 1235
    assert copied["resource_profile_reuse"]["new_calculation_performed"] is False
    ready = json.loads((out / "live_inputs" / key / "task.json").read_text())
    assert Path(ready["input_xyz"]).read_bytes() == (parent / "live_inputs" / key / "input.xyz").read_bytes()
    module.api.verify_sources(new["input_sha256"])


def test_geometry_or_SCF_failure_cannot_be_relabelled_as_workspace_collision(case):
    module, parent, out, key, _ = case
    p = parent / "energy/failures" / (key + ".json")
    p.write_text(json.dumps(dict(error="Independent solution gradient is not stationary")))
    with pytest.raises(ValueError, match="storage collision"):
        module.prepare(parent, out)
    assert not out.exists()
