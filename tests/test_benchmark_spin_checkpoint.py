import importlib.util
import json
from pathlib import Path

import pytest


def load():
    path = Path(__file__).resolve().parents[1] / "scripts/validation/complete_benchmark_spin.py"
    spec = importlib.util.spec_from_file_location("benchmark_spin_checkpoint_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_checkpoint_preserves_existing_callback_settings_and_false_convergence(tmp_path, monkeypatch):
    from redox.qm import dft
    from redox.workflows import geometry_checkpoints

    module = load()
    module.OUT = tmp_path
    observed, forwarded = [], []

    def factory(folder):
        def callback(environment):
            observed.append((folder, environment))
        return callback

    def kernel(method, **kwargs):
        forwarded.append({k: v for k, v in kwargs.items() if k != "callback"})
        kwargs["callback"]({"point": "evaluated"})
        return False, "unqualified_molecule"

    monkeypatch.setattr(dft, "_geomopt_fn", lambda backend: kernel)
    previous = dft._geomopt_fn
    monkeypatch.setattr(geometry_checkpoints, "checkpoint_callback", factory)
    monkeypatch.setattr(geometry_checkpoints, "write_geometry", lambda molecule, path, comment: path.write_text(molecule))
    old_calls = []
    with module.checkpoint_optimizers():
        result = dft._geomopt_fn("gpu")("method", maxsteps=120, callback=old_calls.append, convergence_gmax=5e-5)
    assert result == (False, "unqualified_molecule")
    assert forwarded == [dict(maxsteps=120, convergence_gmax=5e-5)]
    assert old_calls == [{"point": "evaluated"}]
    assert len(observed) == 1
    assert json.loads((tmp_path / "geometry_checkpoints/optimizer_0/optimizer_returned.json").read_text())["optimizer_converged"] is False
    assert dft._geomopt_fn is previous


def test_checkpoint_restores_driver_after_failure(tmp_path, monkeypatch):
    from redox.qm import dft
    from redox.workflows import geometry_checkpoints

    module = load()
    module.OUT = tmp_path
    monkeypatch.setattr(geometry_checkpoints, "checkpoint_callback", lambda folder: lambda event: None)

    def fail(*args, **kwargs):
        raise RuntimeError("real optimizer failure")

    monkeypatch.setattr(dft, "_geomopt_fn", lambda backend: fail)
    previous = dft._geomopt_fn
    with pytest.raises(RuntimeError, match="real optimizer failure"):
        with module.checkpoint_optimizers():
            dft._geomopt_fn("gpu")("method")
    assert dft._geomopt_fn is previous
