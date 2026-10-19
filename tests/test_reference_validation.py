import importlib.util
import sys
from types import SimpleNamespace

import pytest

from redox.core.common import ROOT
from redox.properties import potentials
from redox.validation import readiness


@pytest.mark.parametrize("change", [dict(n_imag=3), dict(g_thermal_eV=None),
                                      dict(thermal_definition="different_model")])
def test_invalid_fc_reference_cannot_reference_clean_candidates(monkeypatch, change):
    ox = dict(e_smd_eV=5., g_thermal_eV=1., n_imag=0)
    ox.update(change)
    neu = dict(e_smd_eV=0., g_thermal_eV=1., n_imag=0)
    monkeypatch.setattr(potentials, "read_result", lambda gid, st: ox if st == "ox" else neu)
    assert potentials.fc_live() is None
    with pytest.raises(SystemExit, match="thermally valid"):
        potentials._fc_reference()


def test_strict_readiness_rejects_thermal_review(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["readiness", "--strict"])
    monkeypatch.setattr(readiness, "read_manifest", lambda: [])
    monkeypatch.setattr(readiness, "inspect_states", lambda rows: [dict(
        id="q", state="red2", status="REVIEW", issues="imaginary_thermal_modes")])
    with pytest.raises(SystemExit) as exc:
        readiness.main()
    assert exc.value.code == 1


def test_runner_avoids_existing_cpu_workers(monkeypatch):
    spec = importlib.util.spec_from_file_location("run_motifs_test", ROOT / "scripts/pipeline/run_motifs.py")
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    import psutil
    monkeypatch.setattr(runner.os, "sched_getaffinity", lambda _: set(range(8)))
    worker = SimpleNamespace(pid=-1, cmdline=lambda: ['python'], cpu_affinity=lambda: [0, 1])
    kernel = SimpleNamespace(pid=-2, cmdline=lambda: [], cpu_affinity=lambda: list(range(8)))
    monkeypatch.setattr(psutil, "process_iter", lambda: [worker, kernel])
    assert runner.cpu_allocation(2) == [2, 3]
    with pytest.raises(RuntimeError, match="overlaps"):
        runner.cpu_allocation(2, [0, 1])


def test_new_motif_plan_requires_matched_native_thermal_before_readiness(tmp_path):
    spec = importlib.util.spec_from_file_location("run_motifs_plan_test", ROOT / "scripts/pipeline/run_motifs.py")
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    plan = runner.commands(tmp_path / "config.py", ["new_q"], tmp_path, "gpu", 2, 3)
    labels = [entry[0] for entry in plan]
    assert labels.index("sp") < labels.index("native_thermal") < labels.index("coverage")
    thermal = next(argv for label, _, argv in plan if label == "native_thermal")
    assert "ferrocene,new_q" in thermal and "--select-policy" in thermal and "--repair-minimum" in thermal
    dft = next(argv for label, _, argv in plan if label == "dft")
    assert "--no-freq" in dft
