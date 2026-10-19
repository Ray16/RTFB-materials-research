"""Historical energy repair preserves evidence and refuses changed-source promotion."""
import json

import pytest

from redox.core.common import HARTREE_EV
from redox.core.io import atomic_write_json
from redox.core.protocol import ACTIVE_SP, geom_sha1, record_path, load_record
from redox.workflows import revalidate_sp as R


@pytest.fixture
def campaign(tmp_path):
    root = tmp_path / "dft"
    directory = root / "m" / "neu"
    directory.mkdir(parents=True)
    (directory / "opt.xyz").write_text("2\n\nH 0 0 0\nH 0 0 .74\n")
    atomic_write_json(directory / "result.json", dict(charge=0, mult=1))
    path = record_path(directory, 0, 1)
    old = dict(protocol=ACTIVE_SP, protocol_hash=path.stem, charge=0, mult=1,
               geom_sha1=geom_sha1(directory / "opt.xyz"), source="adopted:result.json",
               density_fit=True, conv_tol=1e-9, converged_smd=True, converged_gas=True,
               e_smd_Ha=-1., e_smd_eV=-HARTREE_EV, e_gas_Ha=-.9, e_gas_eV=-.9*HARTREE_EV)
    atomic_write_json(path, old)
    run = tmp_path / "run"
    queue = R.prepare(run, rows=[dict(id="m", state="neu")], root=root)
    return run, queue["tasks"][0], old


def finish(run, task, old, difference=0.):
    frozen = R.Path(task["frozen"])
    path = record_path(frozen, 0, 1)
    new = dict(old, source="computed:gpu:independent_revalidation",
               e_smd_eV=old["e_smd_eV"]+difference,
               e_smd_Ha=old["e_smd_Ha"]+difference/HARTREE_EV)
    atomic_write_json(path, new)
    atomic_write_json(run / "done" / f"{task['key']}.json", dict(record_sha256=R.sha(path)))
    return new


def test_freezing_does_not_change_any_canonical_evidence(campaign):
    run, task, _ = campaign
    assert R.inputs_current(task)
    assert load_record(R.Path(task["state_dir"]), 0, 1)["source"] == "adopted:result.json"
    assert not (R.Path(task["frozen"]) / "sp").exists()
    assert R.collect(run)["pending"] == 1


def test_agreeing_record_is_promoted_with_exact_original_archive_and_idempotence(campaign):
    run, task, old = campaign
    before = R.Path(task["old_record"]).read_bytes()
    finish(run, task, old, .0002)
    summary = R.collect(run, promote=True)
    assert summary["verified"] == summary["promoted"] == 1
    assert summary["unresolved"] == 0
    current = load_record(R.Path(task["state_dir"]), 0, 1)
    assert current["source"] == "computed:gpu:independent_revalidation"
    archived = list((R.Path(task["old_record"]).parent / "_superseded").glob("*.json"))
    assert len(archived) == 1 and archived[0].read_bytes() == before
    assert R.collect(run, promote=True)["promoted"] == 1


def test_changed_energy_is_a_review_item_and_does_not_overwrite(campaign):
    run, task, old = campaign
    before = R.Path(task["old_record"]).read_bytes()
    finish(run, task, old, .05)
    summary = R.collect(run, promote=True)
    assert summary["changed"] == summary["unresolved"] == 1
    assert summary["promoted"] == 0
    assert R.Path(task["old_record"]).read_bytes() == before


def test_changed_geometry_blocks_even_an_agreeing_record(campaign):
    run, task, old = campaign
    before = R.Path(task["old_record"]).read_bytes()
    finish(run, task, old)
    (R.Path(task["state_dir"]) / "opt.xyz").write_text("2\n\nH 0 0 0\nH 0 0 .80\n")
    summary = R.collect(run, promote=True)
    assert summary["promoted"] == 0 and summary["unresolved"] == 1
    assert R.Path(task["old_record"]).read_bytes() == before


@pytest.mark.parametrize("value", [float("nan"), 99.])
def test_invalid_hartree_units_do_not_establish_a_verified_candidate(campaign, value):
    run, task, old = campaign
    new = finish(run, task, old)
    path = record_path(R.Path(task["frozen"]), 0, 1)
    new["e_gas_Ha"] = value
    path.write_text(json.dumps(new))
    assert R.verified_candidate(task) is None
    assert not R.done_valid(run, task)


def test_recomputation_bypasses_the_adoption_cache(campaign, monkeypatch):
    run, task, old = campaign
    from redox.qm import dft
    output = dict(old, energy_protocol=ACTIVE_SP,
                  energy_geom_sha1=geom_sha1(R.Path(task["frozen"]) / "opt.xyz"))
    calls = []

    def compute(*args, **kwargs):
        calls.append((args, kwargs))
        return output

    monkeypatch.setattr(dft, "dft_smd", compute)
    monkeypatch.setattr(R, "software_versions", lambda: {})
    assert R.compute_task(task)["source"] == "computed:gpu:independent_revalidation"
    assert len(calls) == 1 and calls[0][1]["do_opt"] is False
    assert R.compute_task(task)["source"] == "computed:gpu:independent_revalidation"
    assert len(calls) == 1


def test_cpu_fallback_keeps_protocol_and_avoids_gpu_metadata_initialization(campaign, monkeypatch):
    run, task, old = campaign
    from redox.qm import dft
    from redox.core import protocol
    output = dict(old, energy_protocol=ACTIVE_SP,
                  energy_geom_sha1=geom_sha1(R.Path(task["frozen"]) / "opt.xyz"))
    original_cleanup, original_versions = dft._free_gpu_pool, protocol.software_versions

    def compute(*args, **kwargs):
        assert kwargs["backend"] == "cpu"
        assert protocol.software_versions()["gpu4pyscf"] is None
        dft._free_gpu_pool()
        return output

    monkeypatch.setattr(dft, "dft_smd", compute)
    monkeypatch.setattr(R, "software_versions", lambda: pytest.fail("CPU must not import GPU metadata"))
    result = R.compute_task(task, backend="cpu")
    assert result["source"] == "computed:cpu:independent_revalidation"
    assert result["protocol"] == ACTIVE_SP and R.verified_candidate(task)
    assert dft._free_gpu_pool is original_cleanup and protocol.software_versions is original_versions


def test_unqualified_host_storage_cannot_claim_a_scientific_task(campaign, monkeypatch):
    run, task, _ = campaign
    from redox.workflows import resources
    monkeypatch.setattr(resources, "verify_gate_owner", lambda: None)
    requests = []
    monkeypatch.setattr(resources, "host_capacity", lambda *a, **kw: requests.append((a, kw)))
    certificate = run / "unqualified.json"
    atomic_write_json(certificate, dict(status="running", control_differences_eV={"smd": 0.}))
    queue = json.loads((run / "queue.json").read_text())
    queue["execution_context"] = dict(density_fit_tensors="host", qualification_path=str(certificate),
                                       qualification_sha256=R.sha(certificate))
    atomic_write_json(run / "queue_host.json", queue)
    with pytest.raises(RuntimeError, match="has not passed its control"):
        R.worker(run, "cannot_start", "queue_host.json")
    assert list((run / "claims").iterdir()) == []
    assert requests[0][1] == {"minimum_ram_gib": 64}


def test_diagnostic_recheck_selects_computed_records_without_relabelling_sources(campaign):
    run, task, old = campaign
    canonical = R.Path(task["old_record"])
    atomic_write_json(canonical, dict(old, source="computed:gpu"))
    before = canonical.read_bytes()
    rows = [dict(id="m", state="neu")]
    queue = R.prepare(run.parent / "diagnostics", rows=rows,
                      root=R.Path(task["state_dir"]).parents[1], selection="missing-smd-diagnostics")
    assert len(queue["tasks"]) == 1 and queue["tasks"][0]["require_solution_diagnostics"]
    assert canonical.read_bytes() == before
    finish(run.parent / "diagnostics", queue["tasks"][0], old)
    assert not R.done_valid(run.parent / "diagnostics", queue["tasks"][0])
    path = record_path(R.Path(queue["tasks"][0]["frozen"]), 0, 1)
    record = json.loads(path.read_text())
    record.update(smd_homo_eV=-5., smd_s_squared=0., smd_spin_contam=0.)
    atomic_write_json(path, record)
    assert R.verified_candidate(queue["tasks"][0]) is not None
    atomic_write_json(canonical, dict(old, source="computed:gpu", **{k:record[k] for k in
                      ("smd_homo_eV", "smd_s_squared", "smd_spin_contam")}))
    assert not R.prepare(run.parent / "already_checked", rows=rows,
                         root=R.Path(task["state_dir"]).parents[1], selection="missing-smd-diagnostics")["tasks"]


@pytest.mark.parametrize("change", [{"smd_spin_contam": .2}, {"smd_homo_eV": float("nan")},
                                   {"smd_s_squared": True}, {"smd_spin_contam": None}])
def test_solution_diagnostic_admission_requires_finite_physical_values(change):
    record = dict(smd_homo_eV=-5., smd_s_squared=.75, smd_spin_contam=0.)
    assert R.solution_diagnostics_ok(record)
    record.update(change)
    assert not R.solution_diagnostics_ok(record)


@pytest.mark.parametrize("tamper", [None, "original", "reproduction", "alternative"])
def test_spin_refresh_requires_intact_independent_evidence(campaign, monkeypatch, tamper):
    from redox.validation import spins
    from redox.qm import dft
    run, task, old = campaign
    directory = R.Path(task["state_dir"])
    monkeypatch.setattr(R, "ROOT", run.parent)
    monkeypatch.setattr(spins, "ROOT", run.parent)
    monkeypatch.setattr(spins, "DFT", directory.parents[1])
    monkeypatch.setattr(spins, "declared", lambda *args: {"smiles": "[H][H]"})
    monkeypatch.setattr(dft, "_opt_protocol_mismatch", lambda *args: False)
    helper = run.parent / "helper"; helper.mkdir()
    (helper / "opt.xyz").write_bytes((directory / "opt.xyz").read_bytes())
    atomic_write_json(helper / "result.json", dict(charge=0, mult=3, spin_validation_status="validated"))
    alt_path = record_path(helper, 0, 3)
    alternative = dict(old, source="computed:gpu", mult=3, protocol_hash=alt_path.stem,
                       e_smd_eV=old["e_smd_eV"]+.7, smd_spin_contam=.01)
    atomic_write_json(alt_path, alternative)
    identity = spins.snapshot("m", "neu")
    parent = spins.comparison_path("m", "neu", identity)
    atomic_write_json(parent, dict(status="spin_gap_verified", chosen_identity=identity, chosen_mult=1,
                                  alternative_mult=3, electronic_gap_eV=.7,
                                  thermochemical_gap_validated=False, helper=str(helper)))
    before = parent.read_bytes()
    finish(run, task, old, .0002)
    assert R.collect(run, promote=True)["promoted"] == 1
    if tamper == "original":
        (R.Path(task["frozen"]) / "previous_record.json").write_text("{}")
    elif tamper == "reproduction":
        record_path(R.Path(task["frozen"]), 0, 1).write_text("{}")
    elif tamper == "alternative":
        atomic_write_json(alt_path, dict(alternative, e_smd_eV=alternative["e_smd_eV"]+.1))
    if tamper:
        with pytest.raises(ValueError):spins.refresh_reproduced_comparison("m", "neu")
        assert list(parent.parent.glob("*.json")) == [parent]
    else:
        result = spins.refresh_reproduced_comparison("m", "neu")
        assert result["status"] == "comparison_refreshed"
        derived = json.loads(R.Path(result["path"]).read_text())
        assert derived["electronic_gap_eV"] == pytest.approx(.6998)
        assert derived["refreshed_from_sha256"] == R.sha(parent)
        assert spins.refresh_reproduced_comparison("m", "neu")["status"] == "current_comparison_present"
    assert parent.read_bytes() == before
