import copy
import importlib.util
from pathlib import Path

import pytest

ROOT=Path(__file__).parents[1]
spec=importlib.util.spec_from_file_location('review_xtb_path_binding',ROOT/'scripts/validation/review_xtb_path_binding.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def point():
    return {'status':'qualified_xtb_relocated_kernel_probe_checked','loaded_library_path_verified':True,
        'quantum_calculations_run':False,'native_state_qualified':False,'full_platform_portability_inferred':False,
        'masses_amu':[1.0,12.0],'thermal':{'G_Ha':.1,'H_Ha':.2}}


def test_identical_kernel_controls_match_without_promoting_native_qualification():
    assert module.compare_probes([point(),point(),point()])==0


@pytest.mark.parametrize('change',['thermal','mass','scope','missing_path_verification'])
def test_changed_kernel_or_promoted_scope_is_rejected(change):
    a,b,c=point(),point(),point()
    if change=='thermal':c['thermal']['G_Ha']+=1e-5
    elif change=='mass':c['masses_amu'][0]+=.01
    elif change=='scope':c['native_state_qualified']=True
    else:c['loaded_library_path_verified']=False
    with pytest.raises(ValueError):module.compare_probes([a,b,c])
