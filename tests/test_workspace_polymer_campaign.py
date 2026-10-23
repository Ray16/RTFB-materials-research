"""Checked evidence can be reused across stages without inventing thermochemistry."""
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture
def case(tmp_path, monkeypatch):
    spec = importlib.util.spec_from_file_location("workspace_campaign", ROOT / "scripts/pipeline/workspace_polymer_campaign.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    source, raw, out = (tmp_path / s for s in ("source", "raw", "successor"))
    key = "existing_pose0"
    for p in (source / "live_inputs" / key, source / "controls", source / "energy/pilot", raw / "pilot_execution"):
        p.mkdir(parents=True, exist_ok=True)
    xyz = source / "live_inputs" / key / "input.xyz"
    xyz.write_text("2\nfixture only\nH 0 0 0\nH 0 0 .74\n")
    task = dict(key=key, charge=0, mult=1)
    (xyz.parent / "task.json").write_text(json.dumps(dict(task, input_xyz=str(xyz), input_sha256={str(xyz): m.api.digest(xyz)})))
    plan = dict(tasks=[task], producer_sha256={}, input_sha256={str(xyz): m.api.digest(xyz)},
                raw_root=str(raw), execution_workspace_GiB=64, maximum_memory_pilot_MiB=65536,
                minimum_large_GPU_MiB=81920, quantum_settings=dict(protocol="unchanged_fixture"),
                thermochemistry=dict(memory_cap_gib=24), solvent="acetonitrile")
    (source / "plan.json").write_text(json.dumps(plan))
    point = dict(status="fixed_geometry_workspace_qualified_larger_reference_unavailable", finished=1235,
                 free_energy_computed=False, optimized_geometry_qualified=False, native_Hessian_qualified=False,
                 observed_peak_GPU_MiB=39000, host_required_bytes=303720282499)
    (source / "footprint_qualification.json").write_text(json.dumps(point))
    energy = source / "energy/pilot/energy_and_gradient.json"
    energy.write_text(json.dumps(dict(scope="fixture point")))
    (source / "energy/pilot/status.json").write_text(json.dumps(dict(status="memory_pilot_passed", pid=42,
        finished=1235, output_sha256={str(energy): m.api.digest(energy)})))
    (raw / "pilot_execution/status.json").write_text(json.dumps(dict(gradient_solver_calls=1, budget_GiB=64,
        direct_LU_solver_retained=True, source_sha256={str(m.ADAPTER): m.api.digest(m.ADAPTER)})))
    monkeypatch.setattr(m, "controls", lambda plan: {})
    return m, source, out, task, xyz


def test_full_stage_reuses_point_evidence_in_fresh_force_storage(case):
    m, source, out, task, _ = case
    old = {str(p): p.read_bytes() for p in source.rglob("*") if p.is_file()}
    new = m.prepare_profile(source, out)
    assert new["tasks"] == [task]
    assert new["raw_root"] != m.api.read_plan(source)["raw_root"]
    assert all(Path(p).read_bytes() == content for p, content in old.items())
    copied = json.loads((out / "energy/pilot/status.json").read_text())
    assert copied["pid"] == 42 and copied["finished"] == 1235
    assert copied["resource_profile_reuse"]["new_calculation_performed"] is False
    assert not list((out / "thermal/done").iterdir())
    m.api.verify_sources(new["input_sha256"])


@pytest.mark.parametrize("change", [dict(status="energy_gradient_resource_pilot_running"),
                                     dict(native_Hessian_qualified=True),
                                     dict(observed_peak_GPU_MiB=70000)])
def test_incomplete_promoted_or_over_budget_point_cannot_seed_full_stage(case, change):
    m, source, out, _, _ = case
    p = source / "footprint_qualification.json"
    data = json.loads(p.read_text())
    data.update(change)
    p.write_text(json.dumps(data))
    with pytest.raises(ValueError, match="Resource point"):
        m.prepare_profile(source, out)
    assert not out.exists()


def test_native_handoff_keeps_energy_parent_and_leaves_G_missing(case, monkeypatch):
    m, source, out, task, xyz = case
    p = source / "energy/tasks" / task["key"] / "status.json"
    p.parent.mkdir(parents=True)
    report = dict(status="state_records_verified", task=task, physical_QC="passed",
                  E_smd_eV=-1., output_sha256={str(xyz): m.api.digest(xyz)})
    p.write_text(json.dumps(report))
    original = p.read_bytes()
    monkeypatch.setattr(m, "checked_energy", lambda source: (report, p))
    m.prepare_native(source, out, source / "controls", 64, 81920)
    admitted, invalid = m.api.admitted_records(out, "energy")
    assert not invalid and admitted[task["key"]]["E_smd_eV"] == -1.
    assert p.read_bytes() == original
    assert "G_eV" not in admitted[task["key"]] and not list((out / "thermal/done").iterdir())


@pytest.mark.parametrize("budget,device", [(64, 32768), (float("nan"), 81920), (-1, 81920)])
def test_native_handoff_requires_valid_workspace_and_device_headroom(case, budget, device):
    m, source, out, _, _ = case
    with pytest.raises(ValueError, match="device headroom"):
        m.prepare_native(source, out, source / "controls", budget, device)
    assert not out.exists()
