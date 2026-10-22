"""Bounded chunk packing, concurrent ownership, integrity and real crash recovery."""

import json
import multiprocessing as mp
import os
from pathlib import Path

import pytest

from redox.workflows import chunked as queue

STAGE = dict(
    name="control", protocol=dict(operation="square"), quality_scope="synthetic"
)


def records(count):
    for index in range(count):
        yield dict(
            identity=dict(index=index),
            payload=dict(index=index),
            provenance=dict(recipe=index),
        )


def square(task):
    index = task["payload"]["index"]
    return dict(index=index, square=index * index)


def claim_all(root, output):
    plan = queue.read_plan(root)
    selected = []
    with queue.owner(root) as owner:
        cursor = 0
        while True:
            chunk, cursor = queue.claim(root, owner, plan["chunks"], cursor)
            if chunk is None:
                break
            selected.append(chunk["key"])
    Path(output).write_text(json.dumps(selected))


def crash_after_checkpoint(root):
    plan = queue.read_plan(root)
    calls = 0

    def handler(task):
        nonlocal calls
        calls += 1
        if calls == 5:
            os._exit(9)
        return square(task)

    with queue.owner(root) as owner:
        chunk, _ = queue.claim(root, owner, plan["chunks"])
        queue.run_chunk(root, chunk, handler, checkpoint_every=1)


def test_stream_dedup_preserves_every_recipe_and_stable_chunk_order(tmp_path):
    data = [*records(11), *records(3)]
    left = queue.prepare(iter(data), tmp_path / "left", STAGE, chunk_size=4)
    right = queue.prepare(reversed(data), tmp_path / "right", STAGE, chunk_size=4)
    assert left["unique_tasks"] == 11 and left["duplicate_recipes"] == 3
    assert [c["sha256"] for c in left["chunks"]] == [
        c["sha256"] for c in right["chunks"]
    ]
    assert sum(c["count"] for c in left["chunks"]) == 11
    assert len((tmp_path / "left/lineage.jsonl").read_text().splitlines()) == 14
    assert queue.read_plan(tmp_path / "left") == left


def test_stage_charge_stereo_and_route_identity_are_not_collapsed():
    base = dict(
        identity=dict(smiles="C[C@H](O)F", charge=0, precursor="A"),
        payload=dict(index=1),
    )
    first = queue.task_record(base, STAGE)["task_id"]
    for changes in [dict(smiles="C[C@@H](O)F"), dict(charge=1), dict(precursor="B")]:
        changed = dict(base, identity=dict(base["identity"], **changes))
        assert queue.task_record(changed, STAGE)["task_id"] != first
    assert (
        queue.task_record(base, dict(STAGE, protocol=dict(operation="other")))[
            "task_id"
        ]
        != first
    )


def test_incomplete_or_changed_manifest_refuses_workers(tmp_path):
    root = tmp_path / "run"
    queue.prepare(records(3), root, STAGE)
    path = root / "plan.json"
    plan = json.loads(path.read_text())
    plan["unique_tasks"] += 1
    path.write_text(json.dumps(plan))
    with pytest.raises(ValueError, match="manifest or executor"):
        queue.read_plan(root)
    (root / "preparation.json").write_text('{"status":"preparing"}')
    with pytest.raises(ValueError, match="packing is not complete"):
        queue.read_plan(root)


def test_nonfinite_payload_preserves_unfinished_preparation(tmp_path):
    root = tmp_path / "run"
    with pytest.raises(ValueError):
        queue.prepare(
            [dict(identity=dict(id=1), payload=dict(v=float("nan")))], root, STAGE
        )
    assert (
        json.loads((root / "preparation.json").read_text())["status"]
        == "interrupted_preparation_preserved"
    )
    assert not (root / "plan.json").exists()


def test_chunk_memory_is_bounded_by_bytes_and_count(tmp_path):
    data = [dict(identity=dict(id=i), payload=dict(text="x" * 100)) for i in range(8)]
    plan = queue.prepare(
        data, tmp_path / "run", STAGE, chunk_size=1000, max_chunk_bytes=750
    )
    assert len(plan["chunks"]) > 1
    assert all(
        (tmp_path / "run" / c["path"]).stat().st_size <= 750 for c in plan["chunks"]
    )
    with pytest.raises(ValueError, match="Single task exceeds"):
        queue.prepare(data, tmp_path / "too_small", STAGE, max_chunk_bytes=1)


def test_historical_results_can_be_read_but_not_executed_with_a_new_producer(
    tmp_path, monkeypatch
):
    root = tmp_path / "run"
    plan = queue.prepare(records(3), root, STAGE)
    replacement = tmp_path / "new_producer.py"
    replacement.write_text("new executor version")
    monkeypatch.setattr(queue, "__file__", str(replacement))
    assert queue.read_plan(root) == plan
    with pytest.raises(ValueError, match="Worker executor differs"):
        queue.read_plan(root, require_current=True)
    (root / "producer_source.py").write_text("changed frozen source")
    with pytest.raises(ValueError, match="Historical executor source"):
        queue.read_plan(root)


def test_concurrent_processes_claim_each_chunk_once(tmp_path):
    root = tmp_path / "run"
    plan = queue.prepare(records(137), root, STAGE, chunk_size=7)
    ctx = mp.get_context("spawn")
    workers = [
        ctx.Process(
            target=claim_all, args=(str(root), str(tmp_path / f"worker{i}.json"))
        )
        for i in range(4)
    ]
    for worker in workers:
        worker.start()
    for worker in workers:
        worker.join(20)
        assert worker.exitcode == 0
    claimed = [
        key
        for path in tmp_path.glob("worker*.json")
        for key in json.loads(path.read_text())
    ]
    assert len(claimed) == len(set(claimed)) == len(plan["chunks"])
    assert not list((root / "pending").iterdir())
    for owner in (root / "workers").iterdir():
        queue.recover_local(root, owner)
    assert len(list((root / "pending").iterdir())) == len(plan["chunks"])


def test_live_claim_cannot_be_recovered_by_age_or_exit_marker(tmp_path):
    root = tmp_path / "run"
    plan = queue.prepare(records(3), root, STAGE)
    with queue.owner(root) as owner:
        queue.claim(root, owner, plan["chunks"])
        (owner / "exit.json").write_text('{"exited":0}')
        with pytest.raises(ValueError, match="still live"):
            queue.recover_local(root, owner)


def test_real_process_death_resumes_checked_tasks_without_repeating_them(tmp_path):
    root = tmp_path / "run"
    plan = queue.prepare(records(9), root, STAGE)
    child = mp.get_context("spawn").Process(
        target=crash_after_checkpoint, args=(str(root),)
    )
    child.start()
    child.join(20)
    assert child.exitcode == 9
    original = queue.journal(root, plan["chunks"][0])
    assert len(original) == 4
    old_owner = next((root / "workers").iterdir())
    assert queue.recover_local(root, old_owner) == [plan["chunks"][0]["key"]]
    called = []

    def resumed(task):
        called.append(task["payload"]["index"])
        return square(task)

    with queue.owner(root) as owner:
        chunk, _ = queue.claim(root, owner, plan["chunks"])
        report = queue.run_chunk(root, chunk, resumed, checkpoint_every=1)
    assert len(called) == 5 and report["completed"] == 9
    current = queue.journal(root, chunk)
    assert all(current[key] == value for key, value in original.items())
    queue.verify_done(root / "done" / (chunk["key"] + ".json"))


def test_failed_tasks_are_preserved_and_not_mistaken_for_scientific_qc(tmp_path):
    root = tmp_path / "run"
    plan = queue.prepare(records(6), root, STAGE)

    def broken(task):
        raise RuntimeError("dependency unavailable")

    with pytest.raises(RuntimeError, match="circuit breaker"):
        queue.run_chunk(root, plan["chunks"][0], broken, checkpoint_every=1)
    assert len(queue.journal(root, plan["chunks"][0])) == 3
    assert not list((root / "done").iterdir())
    report = queue.run_chunk(root, plan["chunks"][0], square)
    assert report["status"] == "executed_with_failures"
    assert report["completed"] == report["failed"] == 3
    assert report["scientific_QC_inferred_from_execution"] is False


def test_mutated_inputs_journals_and_completed_outputs_fail_closed(tmp_path):
    root = tmp_path / "run"
    plan = queue.prepare(records(4), root, STAGE)
    chunk = plan["chunks"][0]
    queue.run_chunk(root, chunk, square)
    done = root / "done" / (chunk["key"] + ".json")
    report = queue.verify_done(done)
    Path(report["results"]).write_text("changed")
    with pytest.raises(ValueError, match="results changed"):
        queue.verify_done(done)
    journal = root / "progress" / (chunk["key"] + ".json")
    value = json.loads(journal.read_text())
    value["rows"] = {}
    journal.write_text(json.dumps(value))
    with pytest.raises(ValueError, match="integrity mismatch"):
        queue.journal(root, chunk)
    (root / chunk["path"]).write_text("changed")
    with pytest.raises(ValueError, match="Task chunk changed"):
        queue.run_chunk(root, chunk, square)


def test_descriptor_adapter_matches_existing_serial_definitions_without_path_growth():
    import sys
    from redox.screening import streaming_search as search
    from redox.screening.discovery import describe
    from redox.core.common import load_config
    from itertools import islice

    examples = list(islice(search.family_records(), 10))
    expected = [
        describe(
            row["payload"]["product_smiles"], 2, load_config("electrolyte").FARADAY
        )
        for row in examples
    ]
    search.descriptor_tools()
    before = len(sys.path)
    for row, baseline in zip(examples, expected):
        actual = search.descriptor(row)
        assert all(actual[key] == value for key, value in baseline.items())
        assert actual["solution_properties_status"] == "not computed"
    assert len(sys.path) == before
