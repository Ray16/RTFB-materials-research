from types import SimpleNamespace
import pytest
from redox.workflows import resources


def test_large_scientific_worker_is_not_reused(monkeypatch):
    import psutil
    monkeypatch.setattr(resources.os,'sched_getaffinity',lambda pid:set(range(64)))
    worker=SimpleNamespace(pid=-1,cmdline=lambda:['crest'],cpu_affinity=lambda:list(range(14)))
    ordinary=SimpleNamespace(pid=-2,cmdline=lambda:['bash'],cpu_affinity=lambda:list(range(64)))
    monkeypatch.setattr(psutil,'process_iter',lambda:[worker,ordinary])
    assert resources.cpu_allocation(14)==list(range(14,28))
    with pytest.raises(RuntimeError,match='overlaps'):resources.cpu_allocation(2,[0,1])
