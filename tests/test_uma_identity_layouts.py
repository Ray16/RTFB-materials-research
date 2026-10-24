"""Fingerprint both wheel and editable registries without importing a GPU stack."""
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

from redox.qm.uma_model_identity import registry_file


@pytest.fixture(autouse=True)
def isolated_registry_module(monkeypatch):
    monkeypatch.delitem(sys.modules, 'fairchem.core.calculate', raising=False)


def distribution(root,editable=None):
    return SimpleNamespace(locate_file=lambda name:root/name,
        read_text=lambda name:json.dumps(dict(url=editable.as_uri(),dir_info=dict(editable=True))) if editable else None)


def test_wheel_registry_is_read_directly(tmp_path):
    registry=tmp_path/'fairchem/core/calculate/pretrained_models.json';registry.parent.mkdir(parents=True)
    registry.write_text('{"model":{}}')
    assert registry_file(distribution(tmp_path))==registry


def test_editable_source_is_discovered_without_import(tmp_path):
    project=tmp_path/'editable project';registry=project/'src/fairchem/core/calculate/pretrained_models.json'
    registry.parent.mkdir(parents=True);registry.write_text('{"model":{}}')
    assert registry_file(distribution(tmp_path/'site-packages',project))==registry


def test_missing_registry_fails_instead_of_initializing_model(tmp_path):
    with pytest.raises(RuntimeError,match='No local'):
        registry_file(distribution(tmp_path))


def test_conflicting_source_and_installed_registry_is_refused(tmp_path):
    project=tmp_path/'editable';other=project/'src/fairchem/core/calculate/pretrained_models.json'
    wheel=tmp_path/'site/fairchem/core/calculate/pretrained_models.json'
    for p in (other,wheel):p.parent.mkdir(parents=True)
    other.write_text('{"model":1}');wheel.write_text('{"model":2}')
    with pytest.raises(RuntimeError,match='Conflicting'):
        registry_file(distribution(tmp_path/'site',project))
