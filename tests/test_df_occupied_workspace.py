import pytest

from redox.core.common import ROOT
from redox.workflows import polymer_campaign as api


def test_cap_preserves_all_orbitals_and_records_each_original_loop_block():
    module = api.load_module('occupied_cap_source_test', ROOT / 'scripts/pipeline/df_occupied_workspace.py')
    body = ('    nocc_per_batch = min(nocc, (nocc_per_batch + 7) & 0xffffff8)\n'
            '            nocc_in_batch = nocc1 - nocc0\n'
            '            evaluate_all_integrals()\n')
    got = module.cap_source(body)
    assert 'min(nocc, _occupied_cap,' in got
    assert '_occupied_trace' in got
    assert got.endswith('            evaluate_all_integrals()\n')
    with pytest.raises(ValueError, match='template changed'):
        module.cap_source(body.replace('nocc1 - nocc0', 'different'))


def test_unqualified_cap_is_rejected_before_gpu_import():
    module = api.load_module('occupied_cap_bound_test', ROOT / 'scripts/pipeline/df_occupied_workspace.py')
    with pytest.raises(ValueError, match='two replayed occupied caps'):
        with module.compact_smd_execution('unused', 28, occupied_cap=128):
            pass
