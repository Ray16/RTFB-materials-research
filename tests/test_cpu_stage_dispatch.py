"""A shared CPU launcher must dispatch route tasks correctly and reject quantum stages."""
import pytest

from redox.workflows import cpu_chunks


def test_exact_route_stage_uses_route_worker_and_quantum_is_refused(tmp_path, monkeypatch):
    assert cpu_chunks.module_for_stage(dict(stage=dict(name='exact_routes'))) == 'redox.screening.route_stage'
    assert cpu_chunks.module_for_stage(dict(stage=dict(name='descriptors'))) == 'redox.screening.streaming_search'
    monkeypatch.setattr(cpu_chunks.chunked, 'read_plan', lambda *a, **kw: dict(stage=dict(name='quantum')))
    monkeypatch.setattr(cpu_chunks.subprocess, 'run', lambda *a, **kw: pytest.fail('No command may launch for a quantum stage'))
    destination = tmp_path / 'refused'
    with pytest.raises(ValueError, match='quantum jobs require'):
        cpu_chunks.launch(tmp_path / 'campaign', {'lambda1': 1}, destination)
    assert not destination.exists()
