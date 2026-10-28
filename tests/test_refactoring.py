"""Regression cases for shared evidence ownership and workflow isolation."""
import csv
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest

from redox.core.common import ROOT
from redox.core.io import atomic_text_writer, atomic_write_json
from redox.core.protocol import ACTIVE_SP, geom_sha1, load_record
from redox.qm import sp
from redox.workflows import tables


def test_scorecard_audit_respects_serialized_reorganization_precision():
    from test_independent_native_audit import audit_module
    audit = audit_module()
    raw = [dict(lam=.30004, lam_O=.10004, lam_R=.2),
           dict(lam=.42614, lam_O=.12614, lam_R=.3)]
    # CSV totals are .3000 and .4261; averaging the raw values would give .3631.
    assert audit.scorecard_inner_lambdas(raw) == (.363, .1131, .25)
    assert audit.scorecard_inner_lambdas([]) == (None, None, None)


def test_interrupted_record_write_keeps_previous_complete_file(tmp_path):
    path = tmp_path / "state.json"
    atomic_write_json(path, {"valid": True})
    before = path.read_bytes()
    with pytest.raises(RuntimeError):
        with atomic_text_writer(path) as stream:
            stream.write('{"partial":')
            assert path.read_bytes() == before
            raise RuntimeError("interrupted producer")
    assert path.read_bytes() == before
    assert list(tmp_path.iterdir()) == [path]
    with pytest.raises(ValueError):
        atomic_write_json(path, {"energy": float("nan")})
    assert path.read_bytes() == before


@pytest.fixture
def adoption(tmp_path):
    (tmp_path / "opt.xyz").write_text("2\n\nH 0 0 0\nH 0 0 .74\n")
    record = dict(energy_protocol=ACTIVE_SP, energy_geom_sha1=geom_sha1(tmp_path / "opt.xyz"),
                  optimized=True, charge=0, mult=1, density_fit=True, conv_tol=1e-9,
                  sp_xc=ACTIVE_SP["xc"], sp_basis=ACTIVE_SP["basis"],
                  sp_nlc=ACTIVE_SP["nlc"], sp_disp=ACTIVE_SP["disp"],
                  solvent=ACTIVE_SP["solvent"], converged_smd=True, converged_gas=True,
                  e_smd_eV=-1.0, e_gas_eV=-0.5)
    return tmp_path, record


@pytest.mark.parametrize("change", [
    {"energy_geom_sha1": "previous_geometry"}, {"charge": 1}, {"mult": 3},
    {"energy_protocol": None}, {"e_smd_eV": float("nan")},
    {"converged_gas": 1}, {"e_gas_eV": True},
])
def test_adoption_cannot_relabel_a_different_or_unverified_source(adoption, change):
    directory, record = adoption
    record.update(change)
    assert not sp._adoptable(record, charge=0, mult=1, geom_hash=geom_sha1(directory / "opt.xyz"))


def test_identity_proven_adoption_needs_no_quantum_calculation(adoption, monkeypatch):
    directory, record = adoption
    atomic_write_json(directory / "result.json", record)
    from redox.qm import dft
    monkeypatch.setattr(dft, "dft_smd", lambda *a, **kw: pytest.fail("unexpected computation"))
    status, result = sp.ensure("hydrogen", "neu", directory, 0, 1, backend="cpu")
    assert status == "adopted"
    assert len(result["adopted_result_sha256"]) == 64
    assert load_record(directory, 0, 1) == result


def test_legacy_source_without_geometry_proof_is_recomputed(adoption, monkeypatch):
    directory, record = adoption
    record.pop("energy_geom_sha1")
    atomic_write_json(directory / "result.json", record)
    from redox.qm import dft
    monkeypatch.setattr(dft, "dft_smd", lambda *a, **kw: record)
    monkeypatch.setattr(sp, "software_versions", lambda: {})
    status, _ = sp.ensure("hydrogen", "neu", directory, 0, 1, backend="cpu")
    assert status == "computed"
    assert load_record(directory, 0, 1) is not None


def test_changed_geometry_during_computation_cannot_publish_under_old_address(adoption, monkeypatch):
    directory, record = adoption
    record.pop("energy_geom_sha1")
    atomic_write_json(directory / "result.json", record)
    from redox.qm import dft

    def changed(*args, **kwargs):
        (directory / "opt.xyz").write_text("2\n\nH 0 0 0\nH 0 0 .80\n")
        return record

    monkeypatch.setattr(dft, "dft_smd", changed)
    monkeypatch.setattr(sp, "software_versions", lambda: {})
    with pytest.raises(RuntimeError, match="geometry changed"):
        sp.ensure("hydrogen", "neu", directory, 0, 1, backend="cpu")
    assert not list(directory.glob("sp/*.json"))


def test_unsupported_solver_protocol_does_not_touch_cached_evidence(adoption):
    directory, _ = adoption
    with pytest.raises(ValueError, match="solver protocol"):
        sp.ensure("hydrogen", "neu", directory, 0, 1, proto=dict(ACTIVE_SP, conv_tol=1e-5))
    assert not (directory / "sp").exists()


def test_standalone_and_validation_builders_serialize_the_shared_manifest(tmp_path):
    # Exercise both real main functions, substituting only structure generation.
    # Delaying the builders makes the original unlocked read/merge race visible.
    program = """
import sys,time
from pathlib import Path
from redox.build import standalone,validation
module=standalone if sys.argv[1]=='standalone' else validation
kind=sys.argv[1]
module.LIBRARY=Path(sys.argv[2]);module.MANIFEST=module.LIBRARY/'manifest.csv'
definition=dict(id=sys.argv[1],name=sys.argv[1],family='test',smiles='C',states=[('neu',0,1,0)])
def build(_):
    time.sleep(.3)
    return ('C',1) if kind=='standalone' else (['conf_00.xyz'],1)
if sys.argv[1]=='standalone':module._load=lambda:[definition];module.build_one=build
else:module._load_validation=lambda which:[definition];module.build_core=build
sys.argv=['builder'];module.main()
"""
    env = dict(os.environ, PYTHONPATH=str(ROOT / "src"), OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1")
    children = [subprocess.Popen([sys.executable, "-c", program, name, str(tmp_path)],
                                env=env, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                for name in ("standalone", "validation")]
    for child in children:
        output, error = child.communicate(timeout=30)
        assert child.returncode == 0, output + error
    with (tmp_path / "manifest.csv").open() as stream:
        rows = list(csv.DictReader(stream))
    assert {row["id"] for row in rows} == {"standalone", "validation"}


def test_failed_review_stops_and_never_mutates_parent_module_paths(tmp_path, monkeypatch):
    from redox.properties import potentials
    original = potentials.RESULTS
    calls = []
    monkeypatch.setattr(tables, "snapshot_inputs", lambda source: {})

    def failed(argv, **kwargs):
        calls.append((argv, kwargs))
        return SimpleNamespace(returncode=7)

    monkeypatch.setattr(tables.subprocess, "run", failed)
    destination = tmp_path / "review"
    report = tables.review_tables(destination, source=tmp_path / "missing_benchmarks")
    assert report["status"] == "failed" and report["failed_stage"] == "potentials"
    assert report["arithmetic_discrepancies"] is None
    assert report["copied_benchmark_inputs"] == []
    assert len(calls) == 1 and calls[0][0][1:3] == ["-m", "redox.workflows.tables"]
    assert potentials.RESULTS == original
    assert json.loads((destination / "review_status.json").read_text()) == report
    with pytest.raises(FileExistsError):
        tables.review_tables(destination)


def test_audit_failure_preserves_discrepancy_count(tmp_path, monkeypatch):
    monkeypatch.setattr(tables, "snapshot_inputs", lambda source: {})
    monkeypatch.setattr(tables, "STAGES", (tables.STAGES[-1],))

    def audit_failed(argv, **kwargs):
        destination = Path(argv[argv.index("--out") + 1])
        atomic_write_json(destination / "stages/audit_result.json", {"arithmetic_discrepancies": 3})
        return SimpleNamespace(returncode=1)

    monkeypatch.setattr(tables.subprocess, "run", audit_failed)
    report = tables.review_tables(tmp_path / "review", source=tmp_path)
    assert report["status"] == "failed" and report["arithmetic_discrepancies"] == 3
    assert report["scientific_accuracy_validated"] is False


def test_inputs_changing_during_review_prevent_a_stable_completion_claim(tmp_path, monkeypatch):
    artifact = tmp_path / "raw.json"; artifact.write_text("before")
    monkeypatch.setattr(tables, "STAGES", (tables.STAGES[0],))
    monkeypatch.setattr(tables, "snapshot_inputs", lambda source: {str(artifact): artifact.read_text()})

    def changed(*args, **kwargs):
        artifact.write_text("after")
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(tables.subprocess, "run", changed)
    report = tables.review_tables(tmp_path / "review", source=tmp_path)
    assert report["status"] == "inputs_changed"
    assert report["inputs_current"] is False and report["changed_inputs"] == [str(artifact)]


def test_failed_audit_also_records_input_changes_without_hiding_discrepancies(tmp_path, monkeypatch):
    artifact = tmp_path / "raw.json"; artifact.write_text("before")
    monkeypatch.setattr(tables, "STAGES", (tables.STAGES[-1],))
    monkeypatch.setattr(tables, "snapshot_inputs", lambda source: {str(artifact): artifact.read_text()})

    def changed_during_failed_audit(argv, **kwargs):
        artifact.write_text("after")
        destination = Path(argv[argv.index("--out") + 1])
        atomic_write_json(destination / "stages/audit_result.json", {"arithmetic_discrepancies": 1})
        return SimpleNamespace(returncode=1)

    monkeypatch.setattr(tables.subprocess, "run", changed_during_failed_audit)
    report = tables.review_tables(tmp_path / "review", source=tmp_path)
    assert report["status"] == "failed" and report["failed_stage"] == "audit"
    assert report["arithmetic_discrepancies"] == 1
    assert report["inputs_current"] is False and report["changed_inputs"] == [str(artifact)]
    assert json.loads((tmp_path / "review/review_status.json").read_text()) == report


def test_snapshot_tracks_consumed_outputs_without_optimizer_scratch(tmp_path, monkeypatch):
    monkeypatch.setattr(tables, 'ROOT', tmp_path)
    monkeypatch.setattr(tables, 'CALCS', tmp_path/'calcs')
    state=tmp_path/'calcs/dft/m/neu';state.mkdir(parents=True)
    for name in ('opt.xyz','result.json','_opt_00.xyz','_seed_00.xyz'):
        (state/name).write_text('{}')
    chunks=tmp_path/'results/reorg_anchors/our_chunks_v2';chunks.mkdir(parents=True)
    (chunks/'chunk0_results.csv').write_text('id,lambda\nm,.3\n')
    inputs=tables.snapshot_inputs(tmp_path/'results')
    assert str(state/'opt.xyz') in inputs and str(state/'result.json') in inputs
    assert str(state/'_opt_00.xyz') not in inputs and str(state/'_seed_00.xyz') not in inputs
    assert str(chunks/'chunk0_results.csv') in inputs


def test_finalizer_rejects_old_compute_switch_before_any_stage():
    proc = subprocess.run(["bash", str(ROOT / "scripts/pipeline/finalize_after_dft.sh")],
                          env=dict(os.environ, SKIP_REORG_CROSS="0"), capture_output=True, text=True)
    assert proc.returncode == 2
    assert "reserved run_motifs.py" in proc.stderr
    assert "data-contract tests" not in proc.stdout


def test_independent_mass_audit_reproduces_serialized_precision_at_a_decimal_tie():
    spec = importlib.util.spec_from_file_location("mass_rounding_audit", ROOT / "scripts/validation/audit/recompute_axes.py")
    audit = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(audit)
    import csv
    rows = list(csv.DictReader((ROOT / "tests/fixtures/library_manifest.csv").open()))
    smiles = next(row["smiles"] for row in rows if row["id"] == "exp_pmdi_pent2_r")
    mass, repeat = audit.scorecard_masses(smiles)
    assert mass == 390.439 and repeat == 402.45
    assert round(repeat, 1) == 402.4


def test_gpu_discovery_uses_gate_policy_and_keeps_all_reported_slots(monkeypatch):
    from redox.workflows import resources
    calls = []

    def listed(argv, **kwargs):
        calls.append(argv)
        return SimpleNamespace(stdout="lambda1 0\nlambda1 2\nlambda13 7\n")

    monkeypatch.setattr(resources.subprocess, "run", listed)
    assert resources.free_slots() == [("lambda1", 0), ("lambda1", 2), ("lambda13", 7)]
    assert calls == [[str(resources.GATE), "list"]]


@pytest.mark.parametrize("memory,valid", [(16384, False), (32510, True)])
def test_full_qm_suitability_requires_validated_gpu_memory(monkeypatch, memory, valid):
    from redox.workflows import resources
    monkeypatch.setattr(resources.subprocess, "run", lambda *a, **kw:
                        SimpleNamespace(stdout=f"2, GPU-uuid, V100, {memory}, 7.0\n"))
    if valid:
        assert resources.gpu_suitability(2)["memory_MiB"] == memory
    else:
        with pytest.raises(RuntimeError, match="32-GB-class"):
            resources.gpu_suitability(2)


def test_gpu_probe_cannot_import_cuda_outside_a_live_reservation(tmp_path, monkeypatch):
    from redox.workflows import gpu
    monkeypatch.delenv("CUDA_VISIBLE_DEVICES", raising=False)
    monkeypatch.setattr(gpu.importlib, "import_module", lambda *a: pytest.fail("unreserved CUDA import"))
    with pytest.raises(RuntimeError, match="live gpu_reserve"):
        gpu._probe_and_exec(tmp_path, ["unused"])


def test_launcher_plans_every_suitable_device_with_disjoint_shards(tmp_path, monkeypatch, capsys):
    from redox.workflows import launch
    monkeypatch.setattr(launch.socket, "gethostname", lambda: "lambda1")
    monkeypatch.setattr(launch, "free_slots", lambda hosts: [("lambda1", i) for i in range(8)])

    def suitable(index):
        if index == 5:
            raise RuntimeError("unvalidated small device")
        return {"index": index}

    monkeypatch.setattr(launch, "gpu_suitability", suitable)
    monkeypatch.setattr(launch, "host_capacity", lambda *a, **kw: {"verified": True})
    monkeypatch.setattr(launch, "cpu_allocation", lambda count: list(range(count)))
    monkeypatch.setattr(launch.subprocess, "Popen", lambda *a, **kw: pytest.fail("plan must not submit"))
    destination = tmp_path / "unsubmitted"
    assert launch.main(["--module", "dft", "--out", str(destination)]) == 0
    report = json.loads(capsys.readouterr().out)
    assert len(report["argv"]) == 7
    assert not destination.exists() and report["workers_submitted"] == 0
    assert report["cpu_sets"] == [[i, i + 1] for i in range(0, 14, 2)]
    shards = [command[command.index("--shard") + 1] for command in report["argv"]]
    assert shards == [f"7:{i}" for i in range(7)]


def test_global_shards_preserve_disjoint_cross_node_numbering(tmp_path):
    from redox.workflows.launch import commands
    jobs = commands("dft", [2, 4], tmp_path, shard_count=5, shard_start=3)
    assert [job[job.index("--shard") + 1] for job in jobs] == ["5:3", "5:4"]
    with pytest.raises(ValueError, match="shard range"):
        commands("dft", [2, 4], tmp_path, shard_count=4, shard_start=3)


def test_gpu_runner_preserves_evidence_and_refuses_occupied_devices(tmp_path, monkeypatch):
    from contextlib import nullcontext
    from redox.workflows import gpu
    monkeypatch.setattr(gpu, "allocation_guard", nullcontext)
    monkeypatch.setattr(gpu, "cpu_allocation", lambda *a: [2, 3])
    monkeypatch.setattr(gpu, "host_capacity", lambda *a: {})
    monkeypatch.setattr(gpu, "gpu_suitability", lambda *a: {"index": 0})
    monkeypatch.setattr(gpu, "free_slots", lambda *a: [])
    monkeypatch.setattr(gpu.subprocess, "Popen", lambda *a, **kw: pytest.fail("occupied GPU cannot run"))
    destination = tmp_path / "run"
    with pytest.raises(RuntimeError, match="occupied or reserved"):
        gpu.run_gpu(destination, 0, ["never"])
    status = json.loads((destination / "execution.json").read_text())
    assert status["status"] == "failed" and status["scientific_outputs_verified"] is False
    before = (destination / "execution.json").read_bytes()
    with pytest.raises(FileExistsError):
        gpu.run_gpu(destination, 0, ["never"])
    assert (destination / "execution.json").read_bytes() == before


def test_command_help_is_discoverable_and_forwards_to_target_without_gpu_use(monkeypatch, capsys):
    from redox import __main__ as cli
    assert cli.main(["--help"]) == 0
    assert "energy-audit" in capsys.readouterr().out
    previous = sys.argv
    observed = []
    monkeypatch.setattr(cli.runpy, "run_module", lambda name, **kw: observed.append((name, list(sys.argv))))
    assert cli.main(["energy-audit", "--only", "ferrocene"]) == 0
    assert observed == [("redox.qm.sp", ["redox.qm.sp", "--audit", "--only", "ferrocene"])]
    assert sys.argv is previous


def test_retired_orop_fanout_cannot_start_an_unreserved_process(monkeypatch):
    from redox.validation import orop as module
    with pytest.raises(RuntimeError, match="unreserved OROP fanout is retired"):
        module.fanout([1, 2])
