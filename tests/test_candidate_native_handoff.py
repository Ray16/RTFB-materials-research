"""The observer must wait for checked energy and preserve handoff identities."""
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def handoff(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location("candidate_handoff", ROOT / "scripts/pipeline/coordinate_candidate_baseline.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    base, source, target = (tmp_path / n for n in ("base", "source", "native"))
    base.mkdir()
    source.mkdir()
    key = "baseline_existing"
    plan = dict(tasks=[dict(key=key + "__pose0", model_key=key)], quantum_settings=dict(max_opt_steps=400),
                thermochemistry=dict(memory_cap_gib=24), solvent="acetonitrile")
    (source / "plan.json").write_text(json.dumps(plan))
    (source / "workspace_energy_status.json").write_text(json.dumps(dict(status="optimization_starting")))
    registry = dict(overrides={key: dict(campaign="original_failed_branch")}, external_continuations={key:
        dict(campaign=str(source), native_target=str(target), native_controls=str(tmp_path / "controls"),
             native_workspace_GiB=64, native_minimum_GPU_MiB=81920)})
    (base / "baseline_recoveries.json").write_text(json.dumps(registry))
    calls = []

    def prepare(s, out, controls, budget, minimum):
        calls.append((s, out))
        out.mkdir()
        new = dict(plan, source_campaign=str(source), inherited_energy_only=True,
                   execution_driver=str(m.WORKSPACE_DRIVER), producer_sha256={}, input_sha256={})
        (out / "plan.json").write_text(json.dumps(new))

    monkeypatch.setattr(m.api, "load_module", lambda *args: SimpleNamespace(prepare_native=prepare))
    monkeypatch.setattr(m.api, "admitted_records", lambda *args: ({key: dict(task=dict(model_key=key))}, []))
    return m, base, source, target, key, plan, calls


def test_observer_does_not_prepare_native_from_an_optimization_marker(handoff):
    m, base, source, target, _, _, calls = handoff
    before = (base / "baseline_recoveries.json").read_bytes()
    assert m.reconcile_external_energy(base) == []
    assert not calls and not target.exists()
    assert (base / "baseline_recoveries.json").read_bytes() == before


def test_checked_parent_handoff_preserves_source_and_registers_pending_native(handoff):
    m, base, source, target, key, _, calls = handoff
    marker = source / "workspace_energy_status.json"
    marker.write_text(json.dumps(dict(status="stationary_electronic_case_checked_native_G_pending")))
    before = {str(p): p.read_bytes() for p in source.iterdir()}
    events = m.reconcile_external_energy(base)
    assert len(calls) == len(events) == 1
    assert all(Path(p).read_bytes() == content for p, content in before.items())
    result = json.loads((base / "baseline_recoveries.json").read_text())
    assert result["overrides"][key]["campaign"] == str(target)
    assert key not in result["external_continuations"]
    assert events[0]["repeated_optimization"] is False and events[0]["native_G_computed_by_handoff"] is False
    assert m.reconcile_external_energy(base) == []


def test_existing_target_with_changed_quantum_model_is_preserved_and_refused(handoff):
    m, base, source, target, _, plan, calls = handoff
    (source / "workspace_energy_status.json").write_text(json.dumps(dict(
        status="stationary_electronic_case_checked_native_G_pending")))
    target.mkdir()
    wrong = dict(plan, quantum_settings=dict(max_opt_steps=1))
    (target / "plan.json").write_text(json.dumps(wrong))
    before = (base / "baseline_recoveries.json").read_bytes()
    with pytest.raises(ValueError, match="chemical/numerical model"):
        m.reconcile_external_energy(base)
    assert not calls and json.loads((target / "plan.json").read_text()) == wrong
    assert (base / "baseline_recoveries.json").read_bytes() == before
