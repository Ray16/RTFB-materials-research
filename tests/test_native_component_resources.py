import copy
import math

import pytest

from redox.workflows.native_component_resources import (
    check_peak, larger_profile, peak_limit, validate_plan,
)


def plan():
    return dict(charged_polymer_native_component=True, atoms=137, response_nao=1319,
                response_nocc=249, charge=4, mult=1, prepared_case={'occupied_sites':[2,4]},
                maximum_host_RSS_GiB=80, maximum_wall_seconds=28800,
                minimum_initial_storage_GiB=38, minimum_live_storage_GiB=32,
                **larger_profile(137,1319,249,4,1,[2,4]))


def test_original_32GiB_limit_and_explicit_larger_pilot_are_separate():
    assert peak_limit({})==31744
    with pytest.raises(ValueError):check_peak({},32000)
    p=plan();assert validate_plan(p)['maximum_durable_output_bytes']==6*2**30
    assert p['estimated_retained_bytes']*1.2<p['maximum_durable_output_bytes']
    assert peak_limit(p)==88000
    check_peak(p,32000)


@pytest.mark.parametrize('field,value',[
    ('charge',2),('response_nao',1320),('mult',3),('maximum_peak_GPU_MiB',96000),
    ('minimum_GPU_memory_MiB',31744),('maximum_durable_output_bytes',4*2**30),
    ('minimum_initial_storage_GiB',36),('minimum_live_storage_GiB',0),
    ('native_resource_class','unreviewed'),('matched_empty_polymer_reference',True),
])
def test_larger_profile_cannot_silently_relax_state_dimensions_or_resource_gates(field,value):
    p=copy.deepcopy(plan());p[field]=value
    with pytest.raises(ValueError):validate_plan(p)


@pytest.mark.parametrize('value',[math.nan,math.inf,-1,0,88001])
def test_nonfinite_or_excess_peaks_fail(value):
    with pytest.raises(ValueError):check_peak(plan(),value)


def test_expanded_RHF_budget_is_explicit_and_does_not_expand_other_stages(monkeypatch):
    from contextlib import contextmanager
    from types import SimpleNamespace,ModuleType
    import sys
    from redox.workflows.native_component_resources import EXPANDED_RHF_CLASS,execution_budget
    p=plan();p.update(larger_profile(137,1319,249,4,1,[2,4],EXPANDED_RHF_CLASS))
    validate_plan(p)
    original=lambda:64*2**30
    rhf=SimpleNamespace(get_avail_mem=original)
    for name in ('gpu4pyscf','gpu4pyscf.df','gpu4pyscf.df.hessian'):
        module=ModuleType(name);module.__path__=[];monkeypatch.setitem(sys.modules,name,module)
    sys.modules['gpu4pyscf.df.hessian'].rhf=rhf
    @contextmanager
    def native_budget():
        rhf.get_avail_mem=lambda:6*2**30
        try:yield {'cap_GiB':24,'protocol_default_cap_GiB':6}
        finally:rhf.get_avail_mem=original
    native=SimpleNamespace(execution_budget=native_budget)
    def _get_veff():return rhf.get_avail_mem()
    def _partial_hess_elec():return rhf.get_avail_mem()
    with execution_budget(native,p) as policy:
        assert _get_veff()==16*2**30 and _partial_hess_elec()==6*2**30
        assert policy['separate_larger_dimension_qualification_required']
        assert not policy['equations_changed'] and not policy['precision_changed']
    assert rhf.get_avail_mem is original
