"""An owned thermal subprocess must not make a repeated CPU review conflict."""

from contextlib import nullcontext
import importlib.util
import os
from pathlib import Path

import pytest


def test_repeated_validation_retains_its_known_cpu_allocation(monkeypatch):
    import redox.workflows.resources as resources

    path = (
        Path(__file__).resolve().parents[1]
        / "scripts/validation/validate_polymer_workflow.py"
    )
    spec = importlib.util.spec_from_file_location("validator_cpu_binding_test", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    affinity = {2, 3}
    allocations = []
    monkeypatch.setattr(resources, "allocation_guard", nullcontext)
    monkeypatch.setattr(resources, "host_capacity", lambda n: None)

    def allocate(n):
        assert not allocations, "Cached owned client must not trigger reallocation"
        allocations.append(n)
        return [3]

    monkeypatch.setattr(resources, "cpu_allocation", allocate)
    monkeypatch.setattr(os, "sched_getaffinity", lambda pid: affinity)

    def bind(pid, cpus):
        affinity.clear()
        affinity.update(cpus)

    monkeypatch.setattr(os, "sched_setaffinity", bind)
    module.bind_validation_cpu()
    module.bind_validation_cpu()
    assert allocations == [1] and affinity == {3}
    affinity.clear()
    affinity.add(2)
    with pytest.raises(RuntimeError, match="affinity changed"):
        module.bind_validation_cpu()
