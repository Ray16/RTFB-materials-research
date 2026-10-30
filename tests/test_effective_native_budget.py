import copy
import importlib.util
from pathlib import Path

import pytest

from redox.workflows.effective_native_namespace import validate_class

ROOT = Path(__file__).resolve().parents[1]


def load(relative, name):
    spec = importlib.util.spec_from_file_location(name, ROOT/relative)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


guard = load('scripts/fleet/bounded_peak_host_python_worker.py', 'effective_peak_guard')
worker = load('scripts/pipeline/run_native_effective_budget.py', 'effective_capacity_adapter')


def fixture():
    task = dict(smiles='CCO', atoms=105, active_basis_functions=2537, charge=0, mult=1)
    plan = dict(tasks=[task], maximum_host_RSS_GiB=80, experimental_method_id='model',
                experimental_native_settings=dict(gradient_limit=.00045), experimental_model=dict(solvent='acetonitrile'),
                quantum_settings=dict(convergence=1e-12), runtime_profile='/original/runtime.json')
    report = dict(status='alternate_native_state_independently_reconstructed', task=copy.deepcopy(task),
                  primary_candidate_admission=False, experimental_method_id='model', peak_host_RSS_bytes=27.5*2**30)
    return plan, copy.deepcopy(plan), report


def test_smaller_guard_retains_the_full_original_numerical_protocol():
    target, reference, report = fixture()
    before = copy.deepcopy(target)
    result = validate_class(target, reference, report, 40)
    assert target == before and result['profile_qualified'] is False
    assert result['effective_bound_GiB'] == 40


@pytest.mark.parametrize('change', ['headroom', 'constitution', 'native_settings', 'runtime'])
def test_unmatched_or_underprovisioned_class_cannot_transfer_budget(change):
    target, reference, report = fixture()
    bound = 40
    if change == 'headroom': bound = 32
    if change == 'constitution': report['task']['smiles'] = 'CCN'
    if change == 'native_settings': target['experimental_native_settings']['gradient_limit'] = .001
    if change == 'runtime': target['runtime_profile'] = '/other/runtime.json'
    with pytest.raises(ValueError):
        validate_class(target, reference, report, bound)


def test_admission_adapter_changes_only_the_exact_native_resource_floor():
    seen, calls = [], []

    def original(threads, **kwargs):
        seen.append((threads, kwargs))

    adapted = worker.capacity_adapter(original, 40, calls)
    adapted(2, minimum_ram_gib=112)
    adapted(1, minimum_ram_gib=8)
    assert seen == [(2, dict(minimum_ram_gib=72)), (1, dict(minimum_ram_gib=8))]
    assert calls == [dict(original_floor_GiB=112, effective_floor_GiB=72)]


def test_linux_peak_units_are_explicit():
    assert guard.high_water_bytes('Name:\tpython\nVmHWM:\t12345 kB\n') == 12345*1024
    with pytest.raises(ValueError):
        guard.high_water_bytes('VmHWM: 12345 MB')
    with pytest.raises(ValueError):
        guard.high_water_bytes('Name: python')
