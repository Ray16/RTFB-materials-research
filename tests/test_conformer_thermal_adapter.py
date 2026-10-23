import importlib.util
from pathlib import Path

import pytest


PATH = Path(__file__).resolve().parents[1] / "scripts/pipeline/polymer_conformer_thermals.py"
spec = importlib.util.spec_from_file_location("conformer_thermals", PATH)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def test_worker_dispatches_every_pose_and_preserves_native_input_identity(monkeypatch, tmp_path):
    records = {
        f"occ_2__pose{i}": {"task": {"key": f"occ_2__pose{i}", "model_key": "occ_2"}, "E_smd_eV": -100 - i}
        for i in range(2)
    }
    monkeypatch.setattr(module.api, "admitted_records", lambda folder, stage: (records if stage == "energy" else {}, []))
    backend = module.backend(tmp_path)
    tasks = backend.phase_tasks(tmp_path, {"tasks": [{}, {}]}, "thermal")
    # Both conformers are sent; the original one-leader-per-graph rule would
    # silently calculate only one. The actual native parent task is unmodified.
    assert len(tasks) == 2
    parents, invalid = backend.campaign_api.selected_graph_records(tmp_path)
    assert not invalid
    for task in tasks:
        assert parents[task["model_key"]]["task"]["key"] == task["key"]
        assert parents[task["model_key"]]["task"]["model_key"] == "occ_2"


def test_missing_or_invalid_energy_parent_blocks_all_native_dispatch(monkeypatch, tmp_path):
    monkeypatch.setattr(module.api, "admitted_records", lambda folder, stage: ({}, [{"error": "changed bytes"}]))
    with pytest.raises(ValueError, match="parents"):
        module.backend(tmp_path).phase_tasks(tmp_path, {"tasks": [{}]}, "thermal")


def test_pilot_cannot_silently_select_a_previously_completed_pose(monkeypatch, tmp_path):
    records = {
        "existing": {"task": {"key": "existing", "active_basis_functions": 3000}},
        "pending": {"task": {"key": "pending", "active_basis_functions": 2900}},
    }
    monkeypatch.setattr(module.api, "admitted_records", lambda folder, stage: (records if stage == "energy" else {"existing": records["existing"]}, []))
    tasks = module.backend(tmp_path).phase_tasks(tmp_path, {"tasks": [{}, {}]}, "thermal")
    assert max(tasks, key=lambda task: task["active_basis_functions"])["key"] == "pending"
