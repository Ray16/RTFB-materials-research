"""Real queue/receipt behavior with injected scientific tasks; no GPU execution."""

import json

import pytest

from redox.workflows import benchmark_native as queue
from redox.workflows.checkpoints import verify_state_checkpoint


def test_endpoint_selection_uses_reused_actual_names_and_deduplicates_waves():
    panel = dict(
        targets=[
            dict(id="tempo", event="ox->rad"),
            dict(id="mol", event="neu->red1"),
            dict(id="mol", event="red1->red2"),
        ]
    )
    reviewed = {"mol:neu": dict(status="ready")}
    assert queue.pending_states(panel, reviewed) == [
        ("mol", "red1"),
        ("mol", "red2"),
        ("tempo", "ox"),
        ("tempo", "rad"),
    ]


@pytest.fixture
def campaign(tmp_path, monkeypatch):
    import redox.workflows.resources as resources
    import redox.workflows.processes as processes

    monkeypatch.setattr(resources, "verify_gate_owner", lambda: None)
    monkeypatch.setattr(processes, "stop_owned_children", lambda: None)
    monkeypatch.setattr(queue.os, "sched_getaffinity", lambda _: {0, 1})
    monkeypatch.setattr(queue.signal, "signal", lambda *a: None)
    for name in queue.THREADS:
        monkeypatch.setenv(name, "2")
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "0")
    source = tmp_path / "source.py"
    source.write_text("frozen")
    for name in ("claims", "done", "tasks", "failures", "workers"):
        (tmp_path / name).mkdir()
    (tmp_path / "plan.json").write_text(
        json.dumps(
            dict(
                stop_new_work_at=queue.time.time() + 1000,
                producer_sha256={str(source): queue.digest(source)},
                jobs={
                    "mol": [dict(id="mol", state="neu"), dict(id="mol", state="red1")]
                },
            )
        )
    )

    def complete(item):
        output = tmp_path / (item["state"] + ".json")
        output.write_text(json.dumps(item))
        return dict(item, output_sha256={str(output): queue.digest(output)})

    monkeypatch.setattr(queue, "complete_state", complete)
    return tmp_path


def test_completed_job_receipt_binds_actual_outputs_and_resume_reuses(
    campaign, monkeypatch
):
    queue.worker(campaign)
    receipt = campaign / "done/mol.json"
    report = verify_state_checkpoint(receipt)
    assert len(report["states"]) == 2
    monkeypatch.setattr(
        queue, "complete_state", lambda _: pytest.fail("complete output rerun")
    )
    queue.worker(campaign)
    assert not list((campaign / "claims").iterdir())
    (campaign / "neu.json").write_text("changed")
    with pytest.raises(ValueError, match="Completed output changed"):
        queue.worker(campaign)


def test_failed_second_state_preserves_first_output_without_done(campaign, monkeypatch):
    original = queue.complete_state

    def failed(item):
        if item["state"] == "red1":
            raise RuntimeError("physical minimum failed")
        return original(item)

    monkeypatch.setattr(queue, "complete_state", failed)
    queue.worker(campaign)
    report = json.loads((campaign / "failures/mol.json").read_text())
    assert report["status"] == "failed_preserved" and len(report["states"]) == 1
    assert (campaign / "neu.json").exists()
    assert not (campaign / "done/mol.json").exists()
    assert not list((campaign / "claims").iterdir())
    monkeypatch.setattr(
        queue,
        "complete_state",
        lambda _: pytest.fail("failed input retried without review"),
    )
    queue.worker(campaign)


def test_interruption_releases_only_own_claim_and_preserves_failure(
    campaign, monkeypatch
):
    def interrupted(item):
        raise KeyboardInterrupt("interrupted")

    monkeypatch.setattr(queue, "complete_state", interrupted)
    with pytest.raises(KeyboardInterrupt):
        queue.worker(campaign)
    assert (campaign / "failures/mol.json").exists()
    assert not list((campaign / "claims").iterdir())


def test_changed_producer_is_rejected_before_claim(campaign):
    (campaign / "source.py").write_text("changed")
    with pytest.raises(ValueError, match="Changed benchmark evidence"):
        queue.worker(campaign)
    assert not list((campaign / "claims").iterdir())


def test_launcher_uses_all_suitable_free_slots_and_rejects_small_cards(
    campaign, monkeypatch
):
    import redox.workflows.resources as resources
    from types import SimpleNamespace

    (campaign / "launches").mkdir()
    plan = json.loads((campaign / "plan.json").read_text())
    plan["jobs"]["another"] = []
    (campaign / "plan.json").write_text(json.dumps(plan))
    monkeypatch.setattr(
        resources, "free_slots", lambda *a: [("node", 0), ("node", 1), ("node", 2)]
    )
    started = []

    def run(argv, **kwargs):
        if any("nvidia-smi" in str(arg) for arg in argv):
            return SimpleNamespace(
                stdout="0, 15360\n1, 32768\n2, 32768\n", returncode=0
            )
        if "--min-available-ram-mib" in argv:
            assert argv[argv.index("--min-available-ram-mib") + 1] == "16384"
            return SimpleNamespace(stdout='{"ok": true}', returncode=0)
        started.append(argv[-1])
        return SimpleNamespace(stdout="1234\n", returncode=0)

    monkeypatch.setattr(queue.subprocess, "run", run)
    queue.launch(campaign)
    assert len(started) == 2
    assert "--gpu 1" in started[0] and "--gpu 2" in started[1]
    assert all(
        "gpu_reserve" not in command for command in started
    )  # shared GPU wrapper owns the gate
