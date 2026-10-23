"""Cavity-specific native admission must replace, not be maxed with, the blanket budget."""
import importlib.util
import inspect
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[1]


def load(name, relative):
    spec = importlib.util.spec_from_file_location(name, ROOT / relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_cavity_specific_budget_preserves_cache_limit_and_host_headroom():
    m = load("native_forecast", "scripts/pipeline/xc_workspace_campaign.py")
    assert m.native_host_forecast(24615) == pytest.approx(216.15862073004246)
    assert m.native_host_forecast(100000) == 256
    with pytest.raises(ValueError):
        m.native_host_forecast(0)


def test_native_launch_uses_the_actual_cavity_budget():
    m = load("native_forecast", "scripts/pipeline/xc_workspace_campaign.py")
    original = load("original_workspace", "scripts/pipeline/workspace_polymer_campaign.py")
    selected = []
    driver = SimpleNamespace()

    def slots(available, promises, budget):
        selected.append(budget)
        return int(.75*available/(budget*2**30))

    driver.additional_host_slots = slots
    driver.main = lambda:driver.additional_host_slots(321*2**30, [], 256)
    api = SimpleNamespace(load_module=lambda *args:driver,
        read_plan=lambda *args:dict(native_host_budget_GiB=216.15862073004246))
    ns = dict(api=api, ROOT=ROOT, SCRIPT=Path("driver.py"), Path=Path,
        sys=SimpleNamespace(argv=["driver.py", "launch", "--campaign", "/tmp/fixture", "--stage", "thermal"]),
        atomic_write_json=lambda *args:None)
    exec(m.host_execution_body(inspect.getsource(original.execution)), ns)
    ns["execution"]()
    assert selected == [pytest.approx(216.15862073004246)]
