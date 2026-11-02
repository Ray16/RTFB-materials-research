"""The larger electronic profile must be explicit; native limits stay intact."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest
from redox.core.common import ROOT


def guard():
    path = ROOT / 'scripts/fleet/bounded_peak_host_python_worker.py'
    spec = importlib.util.spec_from_file_location('test_hybrid_peak_guard', path)
    module = importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize('bound', [160, 240])
def test_larger_bound_requires_explicit_energy_profile(bound):
    with pytest.raises(ValueError): guard().validate_resources(bound, 32, .1, 5)
    guard().validate_resources(bound, 32, .1, 5, 'corrected_hybrid_energy')


@pytest.mark.parametrize('args', [(80, 32, .1, 5), (241, 32, .1, 5), (160, 16, .1, 5), (160, 32, 2, 5)])
def test_hybrid_profile_cannot_bypass_its_exact_limits(args):
    with pytest.raises(ValueError): guard().validate_resources(*args, 'corrected_hybrid_energy')


def test_actual_same_process_guard_runs_the_declared_hybrid_profile(tmp_path):
    script = tmp_path / 'program.py'; script.write_text('import os; print(os.getpid())\n')
    output = tmp_path / 'guard'
    run = subprocess.run([sys.executable, str(ROOT / 'scripts/fleet/bounded_peak_host_python_worker.py'),
        '--out', str(output), '--maximum-rss-gib', '160', '--resource-profile', 'corrected_hybrid_energy',
        str(script)], env=os.environ, capture_output=True, text=True, check=True, timeout=20)
    result = json.loads((output / 'status.json').read_text())
    assert int(run.stdout) == result['pid']
    assert result['maximum_RSS_GiB'] == 160 and result['resource_profile'] == 'corrected_hybrid_energy'
    assert result['status'] == 'effective_memory_guard_finished_profile_review_pending'
    assert result['peak_tree_RSS_bytes'] > 0 and result['self_VmHWM_bytes'] > 0
