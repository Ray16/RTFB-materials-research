"""A reservation race must not retry an actual scientific failure."""
import importlib.util
import hashlib
from types import SimpleNamespace

from redox.core.common import ROOT


def pipeline():
    spec = importlib.util.spec_from_file_location("motif_pipeline", ROOT / "scripts/pipeline/run_motifs.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_gate_refusal_retries_without_losing_attempt_logs(tmp_path, monkeypatch):
    module = pipeline()
    calls = []

    def run(argv, **kwargs):
        calls.append(argv)
        kwargs["stdout"].write("gpu_reserve: gpu7 BUSY/reserved -> not running\n" if len(calls) == 1 else "converged\n")
        return SimpleNamespace(returncode=1 if len(calls) == 1 else 0)

    monkeypatch.setattr(module.subprocess, "run", run)
    monkeypatch.setattr(module.time, "sleep", lambda seconds: None)
    proc, attempts = module.run_stage([str(module.GATE), "run", "7", "--", "compute"],
                                      tmp_path / "sp.log", {}, retry_delays=(0,))
    assert proc.returncode == 0 and len(calls) == 2
    assert attempts[0]["reservation_refused"] and not attempts[1]["reservation_refused"]
    assert (tmp_path / "sp.log").read_text().startswith("gpu_reserve:")
    assert (tmp_path / "sp_gate_retry_1.log").read_text() == "converged\n"


def test_scientific_failure_is_never_retried(tmp_path, monkeypatch):
    module = pipeline()

    def run(argv, **kwargs):
        kwargs["stdout"].write("SCF convergence failed\n")
        return SimpleNamespace(returncode=1)

    monkeypatch.setattr(module.subprocess, "run", run)
    monkeypatch.setattr(module.time, "sleep", lambda seconds: (_ for _ in ()).throw(AssertionError("unexpected retry")))
    proc, attempts = module.run_stage([str(module.GATE), "run", "7", "--", "compute"], tmp_path / "sp.log", {})
    assert proc.returncode == 1 and len(attempts) == 1
    assert not attempts[0]["reservation_refused"]


def test_persistent_contention_has_a_bounded_exit(tmp_path, monkeypatch):
    module = pipeline()

    def run(argv, **kwargs):
        kwargs["stdout"].write("gpu_reserve: gpu7 BUSY/reserved -> not running\n")
        return SimpleNamespace(returncode=1)

    monkeypatch.setattr(module.subprocess, "run", run)
    monkeypatch.setattr(module.time, "sleep", lambda seconds: None)
    proc, attempts = module.run_stage([str(module.GATE), "run", "7", "--", "compute"],
                                      tmp_path / "sp.log", {}, retry_delays=(0, 0))
    assert proc.returncode == 1 and len(attempts) == 3
    assert all(a["reservation_refused"] for a in attempts)


def test_append_only_registry_growth_keeps_new_style_selected_inputs(tmp_path):
    module=pipeline();config=tmp_path/'config.py';config.write_text('# appended registry')
    groups=[dict(id='chosen',frag='[*:1]O',states=[('neu',0,1,0)])]
    previous=dict(selected_definitions_sha256=module.selection_hash('C[*:1]',groups))
    assert module.resume_config_matches(previous,config,'C[*:1]',groups)
    assert not module.resume_config_matches(previous,config,'C[*:1]',[dict(groups[0],frag='[*:1]N')])


def test_legacy_resume_requires_intact_archived_selected_definitions(tmp_path,monkeypatch):
    module=pipeline();monkeypatch.setattr(module,'ROOT',tmp_path)
    groups=[dict(id='chosen',frag='[*:1]O',states=[('neu',0,1,0)])]
    old=f'SCAFFOLD="C[*:1]"\nGROUPS={groups!r}\n'
    digest=hashlib.sha256(old.encode()).hexdigest()
    config=tmp_path/'config.py';config.write_text(old+'GROUPS += [dict(id="other",frag="[*:1]N")]\n')
    previous=dict(input_sha256=digest)
    assert not module.resume_config_matches(previous,config,'C[*:1]',groups)
    archive=tmp_path/'results/discovery/run/config_before'/f'{digest}.py';archive.parent.mkdir(parents=True);archive.write_text(old)
    assert module.resume_config_matches(previous,config,'C[*:1]',groups)
    assert not module.resume_config_matches(previous,config,'C[*:1]',[dict(groups[0],frag='[*:1]N')])
    archive.write_text(old+'# changed bytes')
    assert not module.resume_config_matches(previous,config,'C[*:1]',groups)
