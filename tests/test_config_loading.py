"""Trusted configuration loading preserves execution and useful diagnostics."""
from pathlib import Path

import pytest

from redox.core.common import load_python_config


def test_config_executes_with_its_original_module_name(tmp_path):
    path = tmp_path / 'custom.py'
    path.write_text("SCAFFOLD = 'site'\nGROUPS = []\n")
    module = load_python_config(path)
    assert (module.SCAFFOLD, module.GROUPS, module.__name__) == ('site', [], 'custom')


def test_missing_config_has_actionable_path(tmp_path):
    with pytest.raises(FileNotFoundError, match='--config path'):
        load_python_config(tmp_path / 'missing.py')


def test_configuration_exception_retains_its_original_diagnostic(tmp_path):
    path = tmp_path / 'broken.py'
    path.write_text("raise RuntimeError('configuration detail')\n")
    with pytest.raises(RuntimeError, match='configuration detail'):
        load_python_config(path)
