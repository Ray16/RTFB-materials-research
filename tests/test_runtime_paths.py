from pathlib import Path
import pytest
from redox.core.runtime_paths import reservation_gate, reservation_lock_root


def test_explicit_gate_is_authoritative_and_does_not_fall_back(tmp_path,monkeypatch):
    tool=tmp_path/'gpu_reserve';tool.write_text('#!/bin/sh\nexit 0\n');tool.chmod(0o755)
    monkeypatch.setenv('REDOX_GPU_RESERVE',str(tool))
    assert reservation_gate()==tool
    tool.chmod(0o644)
    with pytest.raises(RuntimeError,match='Authoritative'):reservation_gate()


def test_path_discovery_and_shared_lock_override(tmp_path,monkeypatch):
    monkeypatch.setattr('redox.core.runtime_paths.site_defaults',lambda:{})
    tool=tmp_path/'gpu_reserve';tool.write_text('#!/bin/sh\nexit 0\n');tool.chmod(0o755)
    monkeypatch.delenv('REDOX_GPU_RESERVE',raising=False)
    monkeypatch.setattr('redox.core.runtime_paths.shutil.which',lambda _:str(tool))
    assert reservation_gate()==tool
    monkeypatch.setenv('GPU_LOCK_ROOT',str(tmp_path/'leases'))
    assert reservation_lock_root()==tmp_path/'leases'


def test_cpu_only_import_can_declare_unavailable_gate(tmp_path,monkeypatch):
    monkeypatch.setattr('redox.core.runtime_paths.site_defaults',lambda:{})
    monkeypatch.delenv('REDOX_GPU_RESERVE',raising=False)
    monkeypatch.setattr('redox.core.runtime_paths.shutil.which',lambda _:None)
    monkeypatch.setattr(Path,'home',classmethod(lambda cls:tmp_path))
    assert reservation_gate(required=False)==tmp_path/'bin/gpu_reserve'
    with pytest.raises(RuntimeError):reservation_gate()


def test_shared_site_gate_works_with_a_different_node_home(tmp_path,monkeypatch):
    tool=tmp_path/'shared_gpu_reserve';tool.write_text('#!/bin/sh\nexit 0\n');tool.chmod(0o755)
    monkeypatch.delenv('REDOX_GPU_RESERVE',raising=False);monkeypatch.delenv('GPU_LOCK_ROOT',raising=False)
    monkeypatch.setattr('redox.core.runtime_paths.site_defaults',lambda:{'gpu_reserve':str(tool),'gpu_lock_root':str(tmp_path/'shared_leases')})
    monkeypatch.setattr('redox.core.runtime_paths.shutil.which',lambda _:None)
    monkeypatch.setattr(Path,'home',classmethod(lambda cls:tmp_path/'different_home'))
    assert reservation_gate()==tool
    assert reservation_lock_root()==tmp_path/'shared_leases'
