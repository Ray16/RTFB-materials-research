from pathlib import Path

import pytest

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


@pytest.mark.parametrize('spin,expression', [('rhf', 'natm*3*nao*nao'), ('uhf', 'natm*3*2*nao**2')])
def test_transform_changes_only_storage_and_records_allocation(spin, expression):
    module = api.load_module('DF_workspace_transform_test', ROOT / 'scripts/pipeline/df_derivative_workspace.py')
    original = (f'    mem_sufficient = {expression} < word_avail * .6\n'
        '    batch_size = int(word_avail*.8 / _unit)\n'
        "    assert blksize > 0, 'Insufficient GPU memory'\n"
        '    exact_integral_contraction()\n')
    changed = module.transformed_source(original, spin)
    assert 'mem_sufficient = False' in changed
    assert 'batch_size = int(word_avail*.25 / _unit)' in changed
    assert '_workspace_trace' in changed
    assert changed.endswith('    exact_integral_contraction()\n')
    with pytest.raises(ValueError, match='template changed'):
        module.transformed_source(original.replace('word_avail*.8', 'word_avail*.7'), spin)


def test_unsupported_kernel_fails_closed():
    module = api.load_module('DF_workspace_spin_guard', ROOT / 'scripts/pipeline/df_derivative_workspace.py')
    with pytest.raises(ValueError, match='restricted/unrestricted'):
        module.transformed_source('', 'other')
